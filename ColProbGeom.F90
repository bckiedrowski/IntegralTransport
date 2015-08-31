module ColProbGeom
  implicit none

  private
  public  :: initialize_geom

  real(8), parameter :: pi = acos(-1.d0)

  type, public, abstract :: geom_type
    real(8), allocatable :: SigmaT(:)
    real(8), allocatable :: pScatter(:)
    real(8), allocatable :: SigmaF(:)
    real(8), allocatable :: nubar(:)
    real(8), allocatable :: nuSigmaF(:)
    CONTAINS
      procedure,              non_overridable :: fission_matrix
      procedure(meshsize),    deferred        :: mesh_size
      procedure(colprob),     deferred        :: collision_probability
      procedure(reflgeom),    deferred        :: reflect
      procedure(coarsengeom), deferred        :: coarsen_mesh
      procedure (assigngeom), deferred        :: assign_geom
      generic :: assignment(=) => assign_geom
  end type geom_type

  abstract interface
    subroutine colprob( geom, G )
      import geom_type
      class(geom_type),     intent(in)    :: geom
      real(8), allocatable, intent(inout) :: G(:,:)
    end subroutine colprob
    function meshsize( geom )
      import geom_type
      class(geom_type), intent(in) :: geom
      integer                      :: meshsize
    end function meshsize
    subroutine reflgeom( geom, G )
      import geom_type
      class(geom_type),     intent(inout) :: geom
      real(8), allocatable, intent(inout) :: G(:,:)
    end subroutine reflgeom
    subroutine coarsengeom( geom, G, success )
      import geom_type
      class(geom_type),     intent(inout) :: geom
      real(8), allocatable, intent(inout) :: G(:,:)
      logical,              intent(out)   :: success
    end subroutine coarsengeom
    subroutine assigngeom( gnew, geom )
      import geom_type
      class(geom_type), intent(out) :: gnew
      class(geom_type), intent(in)  :: geom
    end subroutine assigngeom
  end interface

  type, extends(geom_type) :: slab_type
    integer              :: n
    integer              :: fold
    real(8)              :: width, dx
    real(8), allocatable :: x(:)
    logical              :: coarsen
    CONTAINS
      procedure :: assign_geom           => assign_slab
      procedure :: collision_probability => collision_probability_slab
      procedure :: mesh_size             => mesh_size_slab
      procedure :: reflect               => reflect_slab
      procedure :: coarsen_mesh          => coarsen_slab
  end type slab_type

  type, extends(geom_type) :: block_type
    integer              :: nx, ny, nw
    integer              :: fold
    real(8)              :: xmax, ymax, dx, dy, dh
    real(8), allocatable :: x(:), y(:), w(:)
    logical              :: coarsen
    CONTAINS
      procedure :: assign_geom           => assign_block
      procedure :: collision_probability => collision_probability_block !collision_probability_slab
      procedure :: mesh_size             => mesh_size_block
      procedure :: reflect               => reflect_block
      procedure :: coarsen_mesh          => coarsen_block
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
subroutine initialize_geom( geom )
  use Input
  use Utility, only : vector, linspace
  implicit none

  class(geom_type), allocatable, intent(inout) :: geom

  type(slab_type)  :: slab
  type(block_type) :: block

  integer :: i

  if ( ProbDim == 1 ) then
    allocate( geom, source = slab )
    write(*,'(" 1d slab ")')
  elseif ( ProbDim == 2 ) then
    allocate( geom, source = block )
    write(*,'(" 2d block ")')
  else
    write(*,*) 'program supports slabs in 1d or 2d only.' ; stop
  endif

  ! >>>>> initialization
  select type ( geom )
  type is ( slab_type ) 
    geom%n       = NMesh_1D
    geom%width   = SlabWidth_1D
    geom%dx      = geom%width / geom%n
    geom%fold    = merge( Reflect_1D, 0, Reflect_1D == 1 )
    geom%coarsen = coarsen_1D

    call linspace( geom%x, 0.d0, geom%width, geom%n+1 )
  type is ( block_type )
    geom%nx      = NMesh_X_2D 
    geom%ny      = NMesh_Y_2D
    geom%xmax    = SlabWidth_X_2D
    geom%ymax    = SlabWidth_Y_2D
    geom%nw      = NAngles_2D
    geom%dh      = RaySpacing_2D
    geom%fold    = merge( Reflect_2D, 0, Reflect_2D == 1 .or. Reflect_2D == 2 )
    geom%coarsen = coarsen_2D

    call linspace( geom%x, 0.d0, geom%xmax, geom%nx+1 )
    call linspace( geom%y, 0.d0, geom%ymax, geom%ny+1 )
    call linspace( geom%w, pi/(geom%nw+1), pi, geom%nw + 1 )
  end select

  ! set cross sections
  call vector( geom%SigmaT, base_SigmaT, geom%mesh_size() )
  call vector( geom%pScatter, base_pScatter, geom%mesh_size() )
  call vector( geom%SigmaF, base_SigmaF, geom%mesh_size() )
  call vector( geom%nubar, base_nubar, geom%mesh_size() )
  call vector( geom%nuSigmaF, base_nubar*base_SigmaF, geom%mesh_size() )

  ! special modifications to cross sections
  select type ( geom )
  type is ( slab_type )
    ! split slab, set zone between xleft and xright to be pure capture
    if ( split_1D ) then
      do i=1,geom%n
        if ( geom%x(i) >= Split_xleft_1D .and. geom%x(i+1) <= Split_xright_1D ) then
          geom%pScatter(i) = 0.d0
          geom%nuSigmaF(i) = 0.d0
        endif
      enddo
    endif
  type is ( block_type ) 
    if ( nonuniform_2D ) then
      do i=1,geom%ny/2
        geom%SigmaT( 1 + (i-1)*geom%nx : geom%nx/2 + (i-1)*geom%nx ) = 2.d0
        geom%SigmaF( 1 + (i-1)*geom%nx : geom%nx/2 + (i-1)*geom%nx ) = 2.d0*(1.d0 - base_pScatter)
      enddo
      geom%nuSigmaF = geom%nubar * geom%SigmaF
    endif
  end select

end subroutine initialize_geom

!------------------------------------------------------------------------------
subroutine assign_slab( gnew, geom )
  use Utility, only : copy
  implicit none

  class(slab_type), intent(out) :: gnew
  class(geom_type), intent(in)  :: geom

  select type(geom) 
  type is (slab_type) 
    gnew%n       = geom%n      
    gnew%width   = geom%width  
    gnew%dx      = geom%dx     
    gnew%fold    = geom%fold   
    gnew%coarsen = geom%coarsen

    call copy(gnew%x,        geom%x)      
    call copy(gnew%SigmaT,   geom%SigmaT) 
    call copy(gnew%pScatter, geom%pScatter) 
    call copy(gnew%SigmaF,   geom%SigmaF) 
    call copy(gnew%nubar,    geom%nubar) 
    call copy(gnew%nuSigmaF, geom%nuSigmaF) 
  end select

end subroutine assign_slab

!------------------------------------------------------------------------------
subroutine assign_block( gnew, geom )
  use Utility, only : copy
  implicit none

  class(block_type), intent(out) :: gnew
  class(geom_type),  intent(in)  :: geom

  select type(geom) 
  type is (block_type) 
    gnew%nx   = geom%nx   
    gnew%ny   = geom%ny   
    gnew%xmax = geom%xmax 
    gnew%ymax = geom%ymax 
    gnew%nw   = geom%nw   
    gnew%dh   = geom%dh   
    gnew%fold = geom%fold 

    call copy( gnew%x,  geom%x) 
    call copy( gnew%y,  geom%y) 
    call copy( gnew%w,  geom%w) 

    call copy( gnew%x,        geom%x)      
    call copy( gnew%SigmaT,   geom%SigmaT) 
    call copy( gnew%pScatter, geom%pScatter) 
    call copy( gnew%SigmaF,   geom%SigmaF) 
    call copy( gnew%nubar,    geom%nubar) 
    call copy( gnew%nuSigmaF, geom%nuSigmaF) 
  end select

end subroutine assign_block

!------------------------------------------------------------------------------
integer function mesh_size_slab( geom )  result(n)
  implicit none

  class(slab_type), intent(in) :: geom

  n = geom%n

end function mesh_size_slab

!------------------------------------------------------------------------------
integer function mesh_size_block( geom )  result(n)
  implicit none

  class(block_type), intent(in) :: geom

  n = geom%nx * geom%ny

end function mesh_size_block

!------------------------------------------------------------------------------
subroutine fission_matrix( geom, G )
  use Utility,       only : matrix
  use MatrixInverse, only : matinv
  implicit none

  class(geom_type),     intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:) 

  real(8), allocatable :: ICG(:,:)     ! identity matrix - diagonal scattering matrix * col. prob. matrix
  real(8), allocatable :: ICGinv(:,:)  ! inverse of I - CG

  integer :: i, n

  ! >>>> compute first collision probability matrix
  call geom%collision_probability( G )

  ! >>>>> compute scattering matrix if needed
  n = size( G, dim=1 )
  if ( any( geom%pScatter > 0.d0 ) ) then
    call matrix( ICG, n, n )
    call matrix( ICGinv, n, n )
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
    G(i,:) = geom%nuSigmaF(i)/geom%SigmaT(i) * G(i,:)
  enddo

end subroutine fission_matrix

!------------------------------------------------------------------------------
subroutine reflect_slab( geom, G )
  use Utility, only : copy, matrix, linspace
  implicit none

  class(slab_type),     intent(inout)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  real(8), allocatable :: F(:,:)

  integer :: N

  if ( geom%fold == 1 ) then
    write(*,'(" reflected at midplane ")') 
    N = geom%mesh_size() / 2
    geom%n = N
    call linspace( geom%x, 0.d0, 0.5d0*geom%width, N+1 )
    call copy( F, G )

    call matrix( G, N, N )
    G(1:N,1:N) = F(1:N,1:N) + F(2*N:N+1:-1,1:N)
    deallocate( F )
  endif

end subroutine reflect_slab

!------------------------------------------------------------------------------
subroutine reflect_block( geom, G )
  implicit none

  class(block_type),     intent(inout)   :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

end subroutine reflect_block

!------------------------------------------------------------------------------
subroutine coarsen_slab( geom, G, success )
  use Utility, only : copy, matrix, linspace
  implicit none

  class(slab_type),     intent(inout) :: geom
  real(8), allocatable, intent(inout) :: G(:,:)
  logical,              intent(out)   :: success

  real(8), allocatable :: tmp(:,:)

  integer :: i, j, n

  n = geom%mesh_size()
  if ( geom%coarsen .and. mod(n,2) == 0 ) then
    ! coarsen by a factor of two
    call matrix( tmp, n/2, n/2 )
    tmp = 0.d0
    do j=1,n,2
      do i=1,n,2
        tmp(1+i/2,1+j/2) = sum( G(i:i+1,j:j+1) )
      enddo
    enddo  
    call copy( G, tmp )
    deallocate( tmp )

    geom%n  = n/2
    call linspace( geom%x, 0.d0, geom%width, geom%n+1 )
    success = .true.  
  else
    success = .false.
  endif

end subroutine coarsen_slab

!------------------------------------------------------------------------------
subroutine coarsen_block( geom, G, success )
  use Utility, only : copy, matrix, vector, linspace
  implicit none

  class(block_type),    intent(inout) :: geom
  real(8), allocatable, intent(inout) :: G(:,:)
  logical,              intent(out)   :: success

  integer :: i, j, k, l, m, n, p, mi, mj, q

  real(8), allocatable :: tmp(:,:) !, ti(:), tj(:)
  integer, allocatable :: ti(:), tj(:)

  m = geom%nx
  n = geom%ny
  if ( geom%coarsen .and. mod(m,2) == 0 .and. mod(n,2) == 0 ) then
    p = m*n/4
    call matrix( tmp, p, p ) ; tmp = 0.d0
    call vector( ti, 4 ) ; ti = 0.d0
    call vector( tj, 4 ) ; tj = 0.d0

    do i=1,p
      mi = 1 + (i-1)/(m/2) !ceiling( i / ( m/2 ) )
      do j=1,p
        mj = 1 + (j-1)/(m/2) !ceiling( j / ( m/2 ) )
        ! construct index vectors ti and tj
        do k=1,2
          do l=1,2
            q = l + 2*(k-1)
            ti(q) = l + 2*(i-1) + m*( mi + k - 2 )
            tj(q) = l + 2*(j-1) + m*( mj + k - 2 )
          enddo
        enddo

        ! loop over index vectors and compute corasened matrix element i,j
        do k=1,4
          do l=1,4
            tmp(i,j) = tmp(i,j) + G( ti(k), tj(l) )
          enddo
        enddo
      enddo
    enddo

    call copy( G, tmp )
    deallocate( tmp )

    geom%nx = geom%nx / 2
    geom%ny = geom%ny / 2
    call linspace( geom%x, 0.d0, geom%xmax, geom%nx+1 )
    call linspace( geom%y, 0.d0, geom%ymax, geom%ny+1 )

    success = .true.  
  else
    success = .false.  ! not implemented yet
  endif

end subroutine coarsen_block

!------------------------------------------------------------------------------
subroutine collision_probability_slab( geom, G )
  use Utility,     only : matrix
  use ExpIntegral, only : En => ExpIntN
  implicit none

  class(slab_type),     intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  real(8) :: r, s, t, p, dx
  integer :: i, j, n

  n = geom%n ; dx = geom%dx
  call matrix( G, n, n )
  G = 0.d0

  do j=1,n
    s = geom%SigmaT(j) * dx
    G(j,j) = 1.d0 - ( 1.d0 - 2.d0*En(3,s) )/(2.d0*s)
    t = 0.d0
    do i=j+1,n
      r = geom%SigmaT(i) * dx
      p = 0.5d0*( En(3,t) - En(3,t+r) - En(3,t+s) + En(3,t+r+s) )
      G(i,j) = p / s ; G(j,i) = p / r
      t = t + r
    enddo
  enddo

end subroutine collision_probability_slab

!------------------------------------------------------------------------------
subroutine collision_probability_block( geom, G )
  use Utility, only : matrix
  implicit none

  class(block_type),    intent(in)    :: geom
  real(8), allocatable, intent(inout) :: G(:,:)

  type(ray_vector_type) :: rvec

  real(8) :: dx, dy, dh, dw, c
  integer :: i, j, n, M

  dx = geom%xmax / geom%nx
  dy = geom%ymax / geom%ny
  dh = geom%dh
  dw = geom%w(2) - geom%w(1)

  M = geom%nx * geom%ny
  call matrix( G, M, M ) ; G = 0.d0

  n = 0
  do i=1,geom%nw
    ! for each direction, get a vector of rays intersecting geometry
    call rvec%create_rays( geom, i )
    n = n + rvec%n

    ! loop over each ray and perform integration of collision kernel
    do j=1,rvec%n
      call rvec%ray(j)%integrate_ray( geom, G )
    enddo
  enddo

  ! finalize collision probability matrix (normalize and compute diagnonal)
  c = dh*dw/(2.d0*pi*dx*dy)
  do j=1,M
    G(:,j) = G(:,j) * c / geom%SigmaT(j)
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
  do while ( u < geom%xmax )
    rvec%ray(k)%x0 = u
    rvec%ray(k)%y0 = 0.d0
    u = u + du ; k = k + 1
  enddo

  v = 0.0d0
  do while ( v < geom%ymax )
    rvec%ray(k)%x0 = merge( 0.d0, geom%xmax, geom%w(n) < half_pi )
    rvec%ray(k)%y0 = v
    v = v + dv ; k = k + 1
  enddo

  ! all rays have same direction
  rvec%ray(1:rvec%n)%w = geom%w(n)

end subroutine create_rays


subroutine integrate_ray( ray, geom, G )
  use Bickley, only : Kin
  implicit none

  class(ray_type),  intent(in)    :: ray
  type(block_type), intent(in)    :: geom
  real(8),          intent(inout) :: G(:,:)

  type :: elem_list_type
    integer :: k
    real(8) :: s
  end type
  type(elem_list_type) :: e(1:geom%nx*geom%ny)

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
    e(n)%k = i + (j-1)*geom%nx
    e(n)%s = s

    x = x + wx*s
    y = y + wy*s
    if ( s == dx ) then
      i = i+sgn
    else
      j = j+1
    endif
  enddo

  ! loop over intersecting elements and compute collision probabilities for this ray
  ! for efficiency multiplying factors common to all rays intersecting i,j must be applied afterward
  do l=1,n
    ! "within element" collision probability
    j = e(l)%k
    s = e(l)%s * geom%SigmaT(j)
    t = 0.d0
    G(j,j) = G(j,j) + 0.25d0*pi - Kin(3,s) 
    do m=l+1,n
      ! collision probabilities for all other elements (use reciprocity/optic symmetry)
      i = e(m)%k
      r = e(m)%s * geom%SigmaT(i)
      p = Kin(3,t) - Kin(3,t+s) - Kin(3,t+r) + Kin(3,t+r+s) 
      G(i,j) = G(i,j) + p 
      G(j,i) = G(j,i) + p 
      t = t + r
    enddo
  enddo

end subroutine integrate_ray

end module ColProbGeom
