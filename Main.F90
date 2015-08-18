program MutualInfoIntegralTransport
  use Utility,       only : vec, matrix, copy 
  use Input,         only : nIter
  use ColProbGeom,   only : geom_type, initialize_geom
  use Eigen,         only : eigen_type
  implicit none

  type(eigen_type)              :: Eig
  class(geom_type), allocatable :: geom, geom_copy

  real(8), allocatable :: p(:), q(:), r(:), F(:,:), G(:,:), B(:,:)
  real(8), allocatable :: MutualInfo(:), Correl(:)
  real(8) :: Entropy, Hq, DomRatio
  logical :: again
  integer :: i, j, m, nMesh 

  ! >>>>> initialize, calculate fission matrix, apply symmetry if desired
  call initialize_geom( geom )
  call geom%fission_matrix( F )
  call geom%reflect( F )

  ! >>>>> eigenvalues and eigenvectors
  ! find first two eigenvalues and eigenvectors of fission matrix
  call Eig%eig( F, 2 )
  DomRatio = Eig%val(2) / Eig%val(1)
  
  ! scale fission matrix by 1/k and find pmf vector
  F = F / Eig%val(1)
  r = copy( Eig%vec(1)%v ) / sum( Eig%vec(1)%v )

  write(*,'(" mesh            = ", i12)')   geom%mesh_size()
  write(*,'(" keff            = ", f12.5)') Eig%val(1)
  write(*,'(" dominance ratio = ", f12.5)') DomRatio
  write(*,'(" mutual info     = ")') 

  ! >>>>> mutual information calculation
  MutualInfo = vec( nIter )
  Correl     = vec( nIter )
  G = copy( F )
  allocate( geom_copy, source = geom )

  do m=1,nIter
    
    nMesh = geom%mesh_size()
    B = matrix( nMesh, nMesh )
    do i=1,nMesh
      B(i,:) = F(i,:) * r(:)
    enddo  
 
    do   ! until done coarsening
      nMesh = geom%mesh_size()

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

      call geom%coarsen_mesh( B, again )
      if ( .not. again ) exit
    enddo

    if ( m /= nIter ) then
      F = matmul(F,G)
    endif
    geom = geom_copy   ! restore geometry that may have been coarsened
  enddo

  write(*,'("correction factor   = " f12.5)') sqrt( 1.d0 + 2.d0*Correl(1)/(1.d0 - DomRatio) )

  ! >>>>> all done

end program MutualInfoIntegralTransport
