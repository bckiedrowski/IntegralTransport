module Utility
  implicit none

CONTAINS

!-------------------------------------------------------------------------------
function vec( a, N )  result(x)
  implicit none

  real(8), intent(in)  :: a
  integer, intent(in)  :: N
  real(8), allocatable :: x(:)

  if ( allocated(x) )  deallocate( x )
  allocate( x(1:max(1,N)) )
  x(1:N) = a

end function vec

!-------------------------------------------------------------------------------
function linspace( a, b, N )  result(x)
  implicit none

  real(8), intent(in)  :: a, b
  integer, intent(in)  :: N
  real(8), allocatable :: x(:)

  integer :: i

  if ( allocated(x) )  deallocate( x )
  allocate( x(1:max(1,N)) )

  if ( N >= 2 ) then
    do i=1,N
      x(i) = a + real(i-1,8)*(b - a)/(N-1)
    enddo
  else
    x(1) = 0.5d0 * ( b - a )
  endif

end function linspace

!-------------------------------------------------------------------------------
function outer_product( u, v )  result(A)
  implicit none

  real(8), intent(in)  :: u(:), v(:)
  real(8), allocatable :: A(:,:)

  integer :: i, N

  N = size(u)
  if ( N == size(v) ) then
    if ( allocated( A ) )  deallocate( A )
    allocate( A(1:N,1:N) )

    do i=1,N
      A(:,i) = u(:) * v(i)
    enddo

  else
    write(*,*) 'error in outer_product, sizes not equal'
    stop
  endif

end function outer_product

!-------------------------------------------------------------------------------
real(8) function norm( v )   result( z )
  implicit none

  real(8), intent(in) :: v(:)

  integer :: i

  z = 0.d0
  do i=1,size(v)
    z = z + v(i)*v(i)
  enddo
  z = sqrt(z)

end function norm

end module Utility
