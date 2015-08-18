module ColProbSlab
  implicit none

  real(8), parameter :: pi = acos(-1.d0)

  private :: collision_probability_slab, collision_probability_block

  type, abstract :: geom_type
    real(8), allocatable :: pScatter(:)
    real(8), allocatable :: SigmaF(:)
    real(8), allocatable :: nubar(:)
    real(8), allocatable :: nuSigmaF(:)
    CONTAINS
      procedure                  :: collision_probability
      procedure, non_overridable :: fission_matrix
  end type geom_type

  type, extends(geom_type) :: slab_type
    integer              :: n
    real(8)              :: width, dx
    real(8), allocatable :: x(:)
    CONTAINS
      procedure :: collision_probability => collision_probability_slab
  end type slab_type

  type, extends(geom_type) :: block_type
    integer              :: nx, ny, nw
    real(8)              :: xmax, ymax, dx, dy, dh
    real(8), allocatable :: x(:), y(:), w(:)
    CONTAINS
      procedure :: collision_probability => collision_probability_block !collision_probability_slab
  end type block_type

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

!------------------------------------------------------------------------------
subroutine fission_matrix( geom, G )
  use MatrixInverse, only : matinv
  implicit none

  class(geom_type), intent(in)    :: geom
  real(8),          intent(inout) :: G(:,:)     ! collision prob. matrix -> fission matrix

  real(8), allocatable :: ICG(:,:)     ! identity matrix - diagonal scattering matrix * transfer matrix
  real(8), allocatable :: ICGinv(:,:)  ! inverse of I - CG

  integer :: i, n

  ! >>>>> first compute scattering matrix if needed
  n = size( G, dim=1 )
  if ( any( geom%pScatter > 0.d0 ) ) then
    if ( allocated( ICG ) )     deallocate( ICG )
    if ( allocated( ICGinv ) )  deallocate( ICGinv )

    allocate( ICG(1:n,1:n), ICGinv(1:n,1:n) )
    do i=1,n
      ICG(i,:) = -geom%pScatter(i) * G(i,:)
      ICG(i,i) = 1.d0 + ICG(i,i)
    enddo
    call matinv( ICG, ICGinv, n )     ! returns ICGinv, inverse of ICG; destroys input matrix ICG
    G = matmul( G, ICGinv )           ! multiply by initial matrix G to include scattering

    deallocate( ICG, ICGinv )
  endif

  ! >>>>> then compute fission matrix
  ! multiply by nu-sigmaf diagonal matrix to convert to fission matrix
  do i=1,n
    G(i,:) = geom%nuSigmaF(i) * G(i,:)
  enddo

end subroutine fission_matrix

!------------------------------------------------------------------------------
subroutine collision_probability( geom, G ) 
  implicit none

  class(geom_type),     intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  write(*,*) 'geometry has not been initialized.'
  stop

end subroutine collision_probability

!------------------------------------------------------------------------------
subroutine collision_probability_slab( geom, G )
  use ExpIntegral, only : En => ExpIntN
  implicit none

  class(slab_type),     intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  real(8) :: t, p, dx
  integer :: i, j, n

  n = geom%n ; dx = geom%dx
  if ( allocated( G ) )  deallocate( G )
  allocate( G(1:n,1:n) )

  G = 0.d0
  do j=1,n
    G(j,j) = 1.d0 - ( 1.d0 - 2.d0*En(3,dx) )/(2.d0*dx)
    do i=j+1,n
      t = (abs(j-i)-1)*dx
      p = 0.5d0/dx*( En(3,t) - 2.d0*En(3,t+dx) + En(3,t+2.d0*dx) )
      G(i,j) = p ; G(j,i) = p
    enddo
  enddo

end subroutine collision_probability_slab

!------------------------------------------------------------------------------
subroutine collision_probability_block( geom, G )
  use Utility, only : linspace
  implicit none

  class(block_type),    intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  type(ray_vector_type) :: rvec

  real(8) :: dx, dy, dh, dw
  integer :: i, j, n, M

  dx = geom%xmax / geom%nx
  dy = geom%ymax / geom%ny
  dh = geom%dh
  dw = geom%w(2) - geom%w(1)

  M = geom%nx * geom%ny
  if ( allocated( G ) )  deallocate( G )
  allocate( G(1:M,1:M) ) ; G = 0.d0

  n = 0
  do i=1,geom%nw
    ! for each direction, get a vector of rays intersecting geometry
    call rvec%create_rays( geom, i )
    n = n + rvec%n

    ! loop over each ray and perform integration of collision kernel
    do j=1,rvec%n
      call rvec%ray(j)%integrate_ray( geom, G )
    enddo
    !write(*,*) geom%w(i) / pi * 180.d0, rvec%n
  enddo

  ! finalize collision probability matrix (normalize and compute diagnonal)
  G = G*dh*dw/(2.d0*pi*dx*dy)
  do j=1,M
    G(j,j) = 1.d0 - 2.d0*G(j,j) 
  enddo

end subroutine collision_probability_block

!------------------------------------------------------------------------------
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

  u = 0.0d0 ; k = 1
  !u = 0.5d0*du; k = 1
  do while ( u < geom%xmax )
    rvec%ray(k)%x0 = u
    rvec%ray(k)%y0 = 0.d0
    u = u + du ; k = k + 1
  enddo

  v = 0.0d0
  !v = 0.5d0*dv
  do while ( v < geom%ymax )
    rvec%ray(k)%x0 = merge( 0.d0, geom%xmax, geom%w(n) < half_pi )
    rvec%ray(k)%y0 = v
    v = v + dv ; k = k + 1
  enddo

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

  ! construct vector of n elements the current ray intersects 
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

  ! loop over intersecting elements and compute collision probabilities for this ray
  do l=1,n
    ! "within element" collision probability
    j = g(l)%k
    s = g(l)%s
    t = 0.d0
    A(j,j) = A(j,j) + ( 0.25d0*pi - Kin(3,s) ) 
    do m=l+1,n
      ! collision probabilities for all other elements (use reciprocity/optic symmetry)
      i = g(m)%k
      r = g(m)%s
      p = Kin(3,t) - Kin(3,t+s) - Kin(3,t+r) + Kin(3,t+r+s)
      A(i,j) = A(i,j) + p
      A(j,i) = A(j,i) + p
      t = t + r
    enddo
  enddo

end subroutine integrate_ray

end module ColProbSlab
