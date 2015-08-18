program MutualInfoIntegralTransport
  use Utility,       only : linspace, vec, copy_matrix
  use ColProbSlab,   only : geom_type, slab_type, block_type
  use Eigen,         only : eigen_type
  implicit none

  real(8), parameter :: SlabWidth = 20.d0
  integer, parameter :: Nmax      = 1000 !5120 !20480 !5120
  integer, parameter :: nIter     = 1
  real(8), parameter :: CorrelTol = 1.d-6

  integer, parameter :: ProbDim = 2

  logical, parameter :: coarsen = .false.
  logical, parameter :: reflect = .false.
  logical, parameter :: split   = .false.

  real(8), parameter :: base_SigmaT   = 1.0d0
  real(8), parameter :: base_pScatter = 0.9d0
  real(8), parameter :: base_SigmaF   = base_SigmaT * ( 1.d0 - base_pScatter )
  real(8), parameter :: base_nubar    = 1.d0

  real(8), parameter :: xleft = 9.d0, xright = 11.1d0

  real(8), parameter :: pi = acos(-1.d0)

  !>>>>

  type(eigen_type) :: Eig
  type(slab_type)  :: slab
  type(block_type) :: block

  class(geom_type), allocatable :: geom

  real(8) :: dx

  real(8), allocatable :: x(:), p(:), q(:), r(:), F(:,:), G(:,:), B(:,:)

  real(8), allocatable :: MutualInfo(:), Correl(:)
  real(8), allocatable :: tmp(:,:)
  real(8) :: Entropy, Hq, DomRatio

  integer :: nMesh
  integer :: i, j, m, N

  if ( ProbDim == 1 ) then
    allocate( geom, source = slab )
  elseif ( ProbDim == 2 ) then
    allocate( geom, source = block )
  else
    write(*,*) 'program supports slabs in 1d or 2d only.' ; stop
  endif

  select type ( geom )
  type is ( slab_type ) 
    ! >>>>> initialization
    N = Nmax
    x = linspace(0.d0,SlabWidth,N+1)
    dx = SlabWidth / N

    geom%n     = Nmax
    geom%width = SlabWidth
    geom%x     = linspace(0.d0,SlabWidth,N+1)
    geom%dx    = geom%width / geom%n

    geom%SigmaT   = vec( base_SigmaT, N )
    geom%pScatter = vec( base_pScatter, N )
    geom%SigmaF   = vec( base_SigmaF, N )
    geom%nubar    = vec( base_nubar, N )
    geom%nuSigmaF = vec( base_nubar*base_SigmaF, N )

    ! split slab, set zone between xleft and xright to be pure capture
    if ( split ) then
      do i=1,N
        if ( geom%x(i) >= xleft .and. geom%x(i+1) <= xright ) then
          geom%pScatter(i) = 0.d0
          geom%nuSigmaF(i) = 0.d0
        endif
      enddo
    endif
  type is ( block_type )
    geom%nx   = 20
    geom%ny   = 20
    geom%xmax = 5.d0
    geom%ymax = 5.d0
    geom%nw   = 256
    geom%dh   = 0.01d0

    geom%x = linspace( 0.d0, geom%xmax, geom%nx+1 )
    geom%y = linspace( 0.d0, geom%ymax, geom%ny+1 )
    geom%w = linspace( pi/(geom%nw+1), pi, geom%nw + 1 )
 
    N = geom%nx * geom%ny

    geom%SigmaT   = vec( base_SigmaT, N )
    geom%pScatter = vec( base_pScatter, N )
    geom%SigmaF   = vec( base_SigmaF, N )
    geom%nubar    = vec( base_nubar, N )
    geom%nuSigmaF = vec( base_nubar*base_SigmaF, N )
  end select

  call geom%collision_probability( F )
  call geom%fission_matrix( F )

  ! >>>>> reflect slab if desired
  if ( reflect .and. same_type_as(geom,slab) ) then
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
  G = copy_matrix( F )
  
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

      if ( coarsen .and. mod(nMesh,2) == 0 .and. same_type_as(geom,slab) ) then
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

end program MutualInfoIntegralTransport
