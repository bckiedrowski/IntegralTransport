program MutualInfoIntegralTransport
  use Utility,       only : linspace, vec, matrix, copy 
  use Input,         only : coarsen, reflect, nIter
  use ColProbSlab,   only : geom_type, slab_type, block_type, initialize_geom
  use Eigen,         only : eigen_type
  implicit none

  type(eigen_type) :: Eig
  type(slab_type)  :: slab
  type(block_type) :: block

  class(geom_type), allocatable :: geom

  real(8), allocatable :: p(:), q(:), r(:), F(:,:), G(:,:), B(:,:)

  real(8), allocatable :: MutualInfo(:), Correl(:)
  real(8), allocatable :: tmp(:,:)
  real(8) :: Entropy, Hq, DomRatio

  integer :: nMesh
  integer :: i, j, m, N

  ! >>>>> initialize, calculate fission matrix, apply symmetry if desired
  call initialize_geom( geom )
  call geom%fission_matrix( F )
  call geom%reflect( F )

  N = geom%mesh_size()

  ! >>>>> eigenvalues and eigenvectors
  ! find first two eigenvalues and eigenvectors of fission matrix
  call Eig%eig( F, 2 )
  DomRatio = Eig%val(2) / Eig%val(1)
  
  ! scale fission matrix by 1/k and find pmf vector
  F = F / Eig%val(1)
  r = copy( Eig%vec(1)%v ) / sum( Eig%vec(1)%v )

  write(*,'(" mesh            = ", i12)')   N
  write(*,'(" keff            = ", f12.5)') Eig%val(1)
  write(*,'(" dominance ratio = ", f12.5)')  DomRatio
  write(*,'(" mutual info     = ")') 

  ! >>>>> mutual information calculation
  MutualInfo = vec( nIter )
  Correl     = vec( nIter )
  G = copy( F )

  do m=1,nIter

    B = matrix( N, N )
    do i=1,N
      B(i,:) = F(i,:) * r(:)
    enddo  
 
    nMesh = N
    do

      p = vec( nMesh ) ; q = vec( nMesh ) 
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
        tmp = matrix( nMesh/2, nMesh/2 )
        tmp = 0.d0
        do j=1,nMesh,2
          do i=1,nMesh,2
            tmp(1+i/2,1+j/2) = sum( B(i:i+1,j:j+1) )
          enddo
        enddo  
        nMesh = nMesh / 2

        B = copy( tmp )
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
