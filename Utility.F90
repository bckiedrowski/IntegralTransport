module Utility
  implicit none

  public

  interface copy
    module procedure copy_vector, copy_matrix
  end interface
  interface vector
    module procedure new_vector, new_vector_int, uniform_vector
  end interface
  interface matrix
    module procedure new_matrix
  end interface

CONTAINS

!-------------------------------------------------------------------------------
subroutine new_vector( x, N )  
  implicit none

  real(8), allocatable, intent(out) :: x(:)
  integer, intent(in)  :: N

  if ( allocated(x) )  deallocate( x )
  allocate( x(1:max(1,N)) )

end subroutine new_vector

!-------------------------------------------------------------------------------
subroutine new_vector_int( x, N )  
  implicit none

  integer, allocatable, intent(out) :: x(:)
  integer, intent(in)  :: N

  if ( allocated(x) )  deallocate( x )
  allocate( x(1:max(1,N)) )

end subroutine new_vector_int

!-------------------------------------------------------------------------------
subroutine uniform_vector( x, a, N )
  implicit none

  real(8), allocatable, intent(out) :: x(:)
  real(8), intent(in)  :: a
  integer, intent(in)  :: N

  call new_vector( x, N )
  x(1:N) = a

end subroutine uniform_vector

!-------------------------------------------------------------------------------
subroutine new_matrix( A, n, m )  
  implicit none

  integer, intent(in)  :: n, m
  real(8), allocatable :: A(:,:)

  if ( allocated(A) )  deallocate( A )
  allocate( A(1:max(1,n), 1:max(1,m)) )

end subroutine new_matrix

!-------------------------------------------------------------------------------
subroutine linspace( x, a, b, N )
  implicit none

  real(8), allocatable, intent(out) :: x(:)
  real(8), intent(in)  :: a, b
  integer, intent(in)  :: N

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

end subroutine linspace

!-------------------------------------------------------------------------------
subroutine copy_vector( v, u ) 
  implicit none

  real(8), allocatable, intent(out) :: v(:)
  real(8), intent(in)               :: u(:)


  integer :: n

  n = size( u, dim=1 )
  if ( allocated( v ) )  deallocate( v )
  allocate( v(1:n) )
  v(:) = u(:)

end subroutine copy_vector

!-------------------------------------------------------------------------------
subroutine copy_matrix( B, A ) 
  implicit none

  real(8), allocatable, intent(out) :: B(:,:)
  real(8), intent(in)               :: A(:,:)

  integer :: n, m

  n = size( A, dim=1 )
  m = size( A, dim=2 )
  if ( allocated( B ) )  deallocate( B )
  allocate( B(1:n,1:m) )
  B(:,:) = A(:,:)

end subroutine copy_matrix

!-------------------------------------------------------------------------------
function outer_product( u, v )   result(A)
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
