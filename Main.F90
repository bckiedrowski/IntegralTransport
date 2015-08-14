program MutualInfoIntegralTransport
  use Utility,       only : linspace
  use ColProbSlab,   only : slab_type
  use Eigen,         only : eigen_type
  implicit none

  real(8), parameter :: SlabWidth = 20.d0
  integer, parameter :: Nmax      = 1000 !5120 !20480 !5120
  integer, parameter :: nIter     = 1
  real(8), parameter :: CorrelTol = 1.d-6

  logical, parameter :: coarsen = .false.
  logical, parameter :: reflect = .false.
  logical, parameter :: split   = .true.

  ! SigmaT = 1.d0
  real(8), parameter :: base_pScatter = 0.50d0
  real(8), parameter :: base_SigmaF   = 1.d0 * ( 1.d0 - base_pScatter )
  real(8), parameter :: base_nubar    = 1.d0

  real(8), parameter :: xleft = 9.d0, xright = 11.1d0

  !>>>>

  real(8) :: pScatter(1:Nmax) = base_pScatter
  real(8) :: SigmaF(1:Nmax)   = 1.d0 * ( 1.d0 - base_pScatter )
  real(8) :: nubar(1:Nmax)    = base_nubar
  real(8) :: nuSigmaF(1:Nmax) = base_nubar * base_SigmaF

  type(eigen_type) :: Eig
  type(slab_type)  :: slab

  real(8) :: dx

  real(8), allocatable :: x(:), p(:), q(:), r(:), F(:,:), G(:,:), B(:,:)

  real(8), allocatable :: MutualInfo(:), Correl(:)
  real(8), allocatable :: tmp(:,:)
  real(8) :: Entropy, Hq, DomRatio

  integer :: nMesh
  integer :: i, j, m, N

  ! >>>>> initialization
  N = Nmax
  x = linspace(0.d0,SlabWidth,N+1)
  dx = SlabWidth / N

  slab % n     = Nmax
  slab % width = SlabWidth
  slab % x     = linspace(0.d0,SlabWidth,N+1)
  slab % dx    = slab%width / slab%n

  ! split slab, set zone between xleft and xright to be pure capture
  if ( split ) then
    do i=1,N
      if ( slab%x(i) >= xleft .and. slab%x(i+1) <= xright ) then
        pScatter(i) = 0.d0
        nuSigmaF(i) = 0.d0
      endif
    enddo
  endif

  if ( allocated( F ) )  deallocate( F )
  if ( allocated( G ) )  deallocate( G )
  if ( allocated( B ) )  deallocate( B )
  if ( allocated( p ) )  deallocate( p )
  if ( allocated( q ) )  deallocate( q )
  allocate( F(1:N,1:N), G(1:N,1:N), B(1:N,1:N), p(1:N), q(1:N) )

  ! >>>>> calculate fission matrix
  ! first calculate first collision probability matrix and then get fission matrix
  ! loop over source elements
  F = 0.d0
  do j=1,N
    ! loop over destination elements
    do i=j,N
      F(i,j) = slab%collision_probability(i,j)
      F(j,i) = F(i,j)
    enddo
  enddo
  ! apply scattering to first collision probability matrix
  if ( any(pScatter > 0.d0) ) then
    call scattering( F )
  endif

  ! multiply by nu-sigmaf diagonal matrix to convert to fission matrix
  do i=1,n
    F(i,:) = nuSigmaF(i) * F(i,:)
  enddo

  ! >>>>> reflect slab if desired
  if ( reflect ) then
    if ( allocated( G ) )  deallocate( G )
    allocate( G(1:N,1:N) )

    N = N/2
    x = linspace(0.d0,5.d-1*SlabWidth,N+1)
    G = F

    deallocate( F ) ; allocate( F(1:N,1:N) )
    F(1:N,1:N) = G(1:N,1:N) + G(2*N:N+1:-1,1:N)
  endif

  ! >>>>> eigenvalues and eigenvectors
  ! find first two eigenvalues and eigenvectors of fission matrix
  call Eig%eig( F, 2 )
  DomRatio = Eig%val(2) / Eig%val(1)
  
  ! scale fission matrix by 1/k and find pmf vector
  F = F / Eig%val(1)
  p = Eig%vec(1)%v / sum( Eig%vec(1)%v )

  write(*,'(" mesh            = ", i12)')   N
  write(*,'(" keff            = ", f12.5)') Eig%val(1)
  write(*,'(" dominance ratio = ", f12.5)')  DomRatio
  write(*,'(" mutual info     = ")') 

  ! >>>>> mutual information calculation
  if ( allocated( MutualInfo ) )  deallocate( MutualInfo )
  if ( allocated( Correl ) )      deallocate( Correl )
  allocate( MutualInfo(1:nIter), Correl(1:nIter) )
  G = F
  
  if ( allocated( r ) )  deallocate( r )
  allocate( r(1:N) )
  r = p

  do m=1,nIter

    if ( allocated( B ) )  deallocate( B )
    allocate( B(1:N,1:N) )
    do i=1,N
      B(i,:) = F(i,:) * r(:)
    enddo  
 
    nMesh = N
    do
 
      if ( allocated( p ) )  deallocate( p )
      if ( allocated( q ) )  deallocate( q )
      allocate( p(1:nMesh), q(1:nMesh) )
      do i=1,nMesh
        p(i) = sum( B(i,:) )
        q(i) = sum( B(:,i) )
      enddo
    
      Entropy = 0.d0
      do i=1,nMesh
        if ( p(i) > 0.d0 ) then
          Entropy = Entropy - p(i) * log(p(i))
        endif
      enddo

      MutualInfo(m) = 0.d0
      do i=1,nMesh
        do j=1,nMesh
           if ( B(j,i) /= 0.d0 ) then
             MutualInfo = MutualInfo + B(j,i) * log( B(j,i) / ( q(i)*p(j) ) )
           endif
        enddo
      enddo
  
      Correl(m) = sqrt( 1.d0 - exp(-2.d0*MutualInfo(m)) )
  
      write(*,'(2i6,3es14.4)') m, nMesh, Entropy, MutualInfo(m), Correl(m) !/ m

      if ( coarsen .and. mod(nMesh,2) == 0 ) then
        ! coarsen by a factor of two
        allocate( tmp(1:nMesh/2,1:nMesh/2) )
        tmp = 0.d0
        do j=1,nMesh,2
          do i=1,nMesh,2
            tmp(1+i/2,1+j/2) = sum( B(i:i+1,j:j+1) )
          enddo
        enddo  
        nMesh = nMesh / 2

        deallocate( B )
        allocate( B(1:nMesh,1:nMesh) )
        B = tmp
        deallocate( tmp )  
      else  
        exit
      endif

    enddo

    if ( m /= nIter ) then
      F = matmul(F,G)
    endif
  enddo

  write(*,'("correction factor   = " f12.5)') sqrt( 1.d0 + 2.d0*Correl(1)/(1.d0 - DomRatio) )

  ! >>>>> all done
CONTAINS

! convert first collision matrix G into collision source matrix = G*inv( I - C*G )
subroutine scattering( G )
  use MatrixInverse, only : matinv
  implicit none

  real(8), intent(inout) :: G(:,:)     ! transfer matrix -> scattered transfer matrix

  real(8), allocatable :: ICG(:,:)     ! identity matrix - diagonal scattering matrix * transfer matrix
  real(8), allocatable :: ICGinv(:,:)  ! inverse of I - CG

  integer :: i, n

  n = size( G, dim=1 )
  if ( allocated( ICG ) )     deallocate( ICG )
  if ( allocated( ICGinv ) )  deallocate( ICGinv )

  allocate( ICG(1:n,1:n), ICGinv(1:n,1:n) )
  do i=1,n
    ICG(i,:) = -pScatter(i) * G(i,:)
    ICG(i,i) = 1.d0 + ICG(i,i)
  enddo
  call matinv( ICG, ICGinv, n )     ! returns ICGinv, inverse of ICG; destroys input matrix ICG
  G = matmul( G, ICGinv )           ! multiply by initial matrix G to include scattering

end subroutine scattering

end program MutualInfoIntegralTransport
