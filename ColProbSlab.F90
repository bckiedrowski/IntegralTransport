module ColProbSlab
  implicit none

  real(8), parameter :: pi = acos(-1.d0)

  private :: collision_probability

  type :: slab_type
    integer              :: n
    real(8)              :: width, dx
    real(8), allocatable :: x(:)
    CONTAINS
      procedure :: collision_probability
  end type

  type :: block_type
    integer              :: nx, ny, nw
    real(8)              :: xmax, ymax, dx, dy, dh
    real(8), allocatable :: x(:), y(:), w(:)
  end type

  type :: ray_type
    real(8) :: x0, y0, w
    CONTAINS
      procedure :: integrate_ray
  end type ray_type
  type :: ray_vector_type
    integer                     :: n
    type(ray_type), allocatable :: ray(:)
    CONTAINS
      procedure :: create_rays
  end type ray_vector_type

CONTAINS

subroutine test2d( G )
  use Utility, only : linspace
  implicit none

  real(8), allocatable, intent(inout) :: G(:,:)

  type(block_type) :: block

  type(ray_vector_type) :: rvec


  real(8) :: dx, dy, dh, dw
  integer :: i, j, n, M

  block%nx   = 40
  block%ny   = 40
  block%xmax = 5.d0
  block%ymax = 5.d0

  block%x = linspace( 0.d0, block%xmax, block%nx+1 )
  block%y = linspace( 0.d0, block%ymax, block%ny+1 )

  block%nw = 512 !256
  block%dh = 0.001d0
  block%w  = linspace( pi/(block%nw+1), pi, block%nw + 1 )

  dx = block%xmax / block%nx
  dy = block%ymax / block%ny
  dh = block%dh
  dw = block%w(2) - block%w(1)

  M = block%nx*block%ny
  if ( allocated( G ) )  deallocate( G )
  allocate( G(1:M,1:M) ) ; G = 0.d0

  n = 0
  do i=1,block%nw
    call rvec%create_rays( block, i )
    n = n + rvec%n

    ! loop over each ray and perform integration
    do j=1,rvec%n
      call rvec%ray(j)%integrate_ray( block, G )
    enddo

    write(*,*) block%w(i) / pi * 180.d0, rvec%n
  enddo

  !G = dh*dw*G/(pi*dx*dy)

  G = G*dh*dw/(2.d0*pi*dx*dy)
  do j=1,M
    G(j,j) = 1.d0 - 2.d0*G(j,j) !2.d0*dx*dy*G(j,j)
!!    do i=j+1,M
!!      G(j,i) = G(i,j)
!!    enddo
  enddo

!  do j=1,M
!    !write(*,'(25es8.1)')  G(:,j)
!    !write(*,'(es12.4)')  G(j,1)
!  enddo
!  write(*,*) sum(G(:,1))
!  write(*,*) n

end subroutine test2d

real(8) function collision_probability( slab, i, j )  result(p)
  use ExpIntegral, only : En => ExpIntN
  implicit none

  class(slab_type), intent(in) :: slab
  integer,          intent(in) :: i, j

  real(8) :: dx
  real(8) :: t  ! optical depth between elements j and i

  dx = slab%dx
  if ( i == j ) then
    p = 1.d0 - ( 1.d0 - 2.d0*En(3,dx) )/(2.d0*dx)
  else
    t = (abs(j-i)-1)*dx
    p = 0.5d0/dx*( En(3,t) - 2.d0*En(3,t+dx) + En(3,t+2.d0*dx) )
  endif

end function collision_probability


subroutine create_rays( rvec, geom, n )
  implicit none

  class(ray_vector_type), intent(out) :: rvec
  class(block_type),      intent(in)  :: geom
  integer,                intent(in)  :: n

  real(8), parameter :: half_pi = 0.5d0*acos(-1.d0)

  real(8) :: du, dv, u, v

  integer :: k

  if ( sin( geom%w(n) ) == 0.d0 .or. cos( geom%w(n) ) == 0.d0 ) then
    write(*,*) 'error: ray may not be orthogonal to x-y axis.'
    stop
  endif

  du = abs( geom%dh / sin( geom%w(n) ) )
  dv = abs( geom%dh / cos( geom%w(n) ) )

  ! compute number of rays and allocate structure
  rvec%n = 0
  u = 0.0d0 ; v = 0.0d0
  !u = 0.5d0*du ; v = 0.5d0*dv
  do while( u < geom%xmax ) 
    rvec%n = rvec%n + 1 
    u = u + du
  enddo
  do while( v < geom%ymax ) 
    rvec%n = rvec%n + 1 
    v = v + dv
  enddo

  if ( allocated( rvec%ray ) )  deallocate( rvec%ray )
  allocate( rvec%ray(1:rvec%n) )

  !u = 0.0d0 ; k = 1
  u = 0.5d0*du; k = 1
  do while ( u < geom%xmax )
    rvec%ray(k)%x0 = u
    rvec%ray(k)%y0 = 0.d0
    u = u + du ; k = k + 1
  enddo

  !v = 0.0d0
  v = 0.5d0*dv
  do while ( v < geom%ymax )
    rvec%ray(k)%x0 = merge( 0.d0, geom%xmax, geom%w(n) < half_pi )
    rvec%ray(k)%y0 = v
    v = v + dv ; k = k + 1
  enddo

!  do k=1,rvec%n
!    write(*,'(i4,2f12.5)') k, rvec%ray(k)%x0, rvec%ray(k)%y0
!  enddo
!  write(*,*)

  ! all rays have same direction
  rvec%ray(1:rvec%n)%w = geom%w(n)

end subroutine create_rays


subroutine integrate_ray( ray, geom, A )
  use Bickley, only : Kin
  implicit none

  class(ray_type),  intent(in)    :: ray
  type(block_type), intent(in)    :: geom
  real(8),          intent(inout) :: A(:,:)

  type :: grid_list_type
    integer :: k
    real(8) :: s
  end type
  type(grid_list_type) :: g(1:geom%nx*geom%ny)

  real(8) :: x, y, x0, y0, xe, ye, dx, dy, wx, wy, r, s, t, p

  integer :: i, j, k, l, m, n, sgn 

  wx  = cos(ray%w) ; wy = sin(ray%w) 
  sgn = int( sign( 1.d0, wx ) )
  k   = int( merge( 0.d0, 1.d0, sgn == -1 ) )
  xe  = merge( 0.d0, geom%xmax, sgn == -1 ) 
  ye  = geom%ymax

  x = min( ray%x0, geom%xmax - 1.d-10 )
  y = ray%y0

  i = 1 + floor( x / geom%xmax * geom%nx )
  j = 1 + floor( y / geom%ymax * geom%ny )
!  k = i + (j-1)*geom%nx

  ! accumulate the rays
  n = 0
  do while( abs(x - xe) > 1.d-10 .and. abs(y - ye) > 1.d-10 )

    dx = ( geom%x(i+k) - x )/wx
    dy = ( geom%y(j+1) - y )/wy 

    s = min( merge( dx, 1.d10, dx > 0.d0), dy )

    n = n + 1
    g(n)%k = i + (j-1)*geom%nx
    g(n)%s = s

    x = x + wx*s
    y = y + wy*s
    if ( s == dx ) then
      i = i+sgn
    else
      j = j+1
    endif
  enddo

  do l=1,n
    ! handle within element collision probability
    j = g(l)%k
    s = g(l)%s
    t = 0.d0
    A(j,j) = A(j,j) + ( 0.25d0*pi - Kin(3,s) )
    do m=l+1,n
      i = g(m)%k
      r = g(m)%s
      p = Kin(3,t) - Kin(3,t+s) - Kin(3,t+r) + Kin(3,t+r+s)
      A(i,j) = A(i,j) + p
      A(j,i) = A(j,i) + p
      t = t + r
    enddo
  enddo
  !write(*,*) 

end subroutine integrate_ray

end module ColProbSlab
