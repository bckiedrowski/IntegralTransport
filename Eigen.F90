module Eigen
  implicit none

  private

  type, private :: eigenvector
    real(8), allocatable :: v(:)
  end type eigenvector

  type, public :: eigen_type
    real(8),           allocatable :: val(:)
    type(eigenvector), allocatable :: vec(:)
    CONTAINS
      procedure :: eig  => DeflationPowerIteration
      procedure :: init => InitEigen
  end type eigen_type

CONTAINS

subroutine PowerIteration( E, A )
  implicit none

  class(eigen_type), intent(inout) :: E
  real(8),           intent(in)    :: A(:,:)

  real(8), allocatable :: Av(:), v(:)

  real(8) :: d, err
  integer :: i, N

  N = size( A, dim=1 )
  if ( allocated(v)  )  deallocate(v)
  if ( allocated(Av) )  deallocate(Av)
  allocate( v(1:N), Av(1:N) )
  Av = 1.d0

  err = 1.d0 ; E%val(1) = 1.d0
  do while ( err > 1.d-9 )
    v = Av
    Av = matmul( A, v )

    E%val(1) = dot_product( Av, v ) / dot_product( v, v )

    Av = Av / maxval( Av )
    err   = sum( abs( v - Av ) )
  enddo
  E%vec(1)%v = Av

end subroutine PowerIteration

!-------------------------------------------------------------------------------
recursive subroutine DeflationPowerIteration( E, A, k, neig )
  use Utility, only : outer_product
  implicit none

  class(eigen_type), intent(out) :: E
  real(8),           intent(in)  :: A(:,:)
  integer,           intent(in)  :: k
  integer, optional, intent(in)  :: neig

  real(8), allocatable :: B(:,:), x(:)
  type(eigen_type)     :: F 

  integer :: i, m, N

  if ( present( neig ) ) then
    m = neig
  else
    m = 1
  endif

  if ( m > 2 ) then
    write(*,*) ' eigenvalue method only computes first two eigenpairs currently.'
    stop
  endif

  !call InitEigen( E, A, k )
  call E%init( A, k )
  N = size( A, dim=1 )

  call PowerIteration( E, A )

  if ( m < k ) then
    if ( allocated( x ) )  deallocate( x )
    if ( allocated( B ) )  deallocate( B )
    allocate( x(1:N), B(1:N,1:N) )

    x = A(1,:) / ( E%val(1) * E%vec(1)%v(1) )
    B = A - E%val(1) * outer_product( x, E%vec(1)%v )

    call DeflationPowerIteration( F, B, k, m+1 )

    E%val(2)        = F%val(1)
    E%vec(2)%v = F%vec(1)%v
  endif

end subroutine DeflationPowerIteration


!-------------------------------------------------------------------------------
subroutine InitEigen( E, A, k )
  implicit none

  class(eigen_type), intent(out) :: E
  real(8),           intent(in)  :: A(:,:)
  integer,           intent(in)  :: k

  integer :: i, N

  if ( size(A,dim=1) /= size(A,dim=2) ) then
    write(*,*) 'error in eig: non-square matrix'
    stop
  endif

  N = size( A, dim=1 )

  if ( allocated(E%vec) )  deallocate(E%vec)
  allocate( E%vec(1:k) )

  do i=1,k
    if ( allocated(E%vec(i)%v) )  deallocate(E%vec(i)%v)
    allocate( E%vec(i)%v(1:N) )
    E%vec(i)%v = 1.d0 
  enddo

  if ( allocated(E%val) )  deallocate(E%val)
  allocate( E%val(1:k) )
  E%val(:) = 0.d0 

end subroutine InitEigen

end module Eigen
