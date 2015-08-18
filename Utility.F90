module Utility
  implicit none

  interface copy
    module procedure copy_vector, copy_matrix
  end interface
  interface vector
    module procedure new_vector, uniform_vector
  end interface
  interface matrix
    module procedure new_matrix
  end interface

CONTAINS

!-------------------------------------------------------------------------------
function new_vector( N )  result(x)
  implicit none

  integer, intent(in)  :: N
  real(8), allocatable :: x(:)

  if ( allocated(x) )  deallocate( x )
  allocate( x(1:max(1,N)) )

end function new_vector

!-------------------------------------------------------------------------------
function uniform_vector( a, N )  result(x)
  implicit none

  real(8), intent(in)  :: a
  integer, intent(in)  :: N
  real(8), allocatable :: x(:)

  x = new_vector( N )
  x(1:N) = a

end function uniform_vector

!-------------------------------------------------------------------------------
function new_matrix( n, m )  result(A)
  implicit none

  integer, intent(in)  :: n, m
  real(8), allocatable :: A(:,:)

  if ( allocated(A) )  deallocate( A )
  allocate( A(1:max(1,n), 1:max(1,m)) )

end function new_matrix

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
function copy_vector( u )  result(v)
  implicit none

  real(8), intent(in) :: u(:)

  real(8), allocatable :: v(:)

  integer :: n

  n = size( u, dim=1 )
  if ( allocated( v ) )  deallocate( v )
  allocate( v(1:n) )
  v(:) = u(:)

end function copy_vector

!-------------------------------------------------------------------------------
function copy_matrix( A )  result(B)
  implicit none

  real(8), intent(in) :: A(:,:)

  real(8), allocatable :: B(:,:)

  integer :: n, m

  n = size( A, dim=1 )
  m = size( A, dim=2 )
  if ( allocated( B ) )  deallocate( B )
  allocate( B(1:n,1:m) )
  B(:,:) = A(:,:)

end function copy_matrix

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
