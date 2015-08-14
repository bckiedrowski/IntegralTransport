module ColProbSlab
  implicit none

  private :: collision_probability

  type :: slab_type
    integer              :: n
    real(8)              :: width, dx
    real(8), allocatable :: x(:)
    CONTAINS
      procedure :: collision_probability
  end type

CONTAINS

real(8) function collision_probability( slab, i, j )  result(p)
  use ExpIntegral, only : En => ExpIntN
  implicit none

  class(slab_type), intent(in) :: slab
  integer,          intent(in) :: i, j

  real(8) :: dx
  real(8) :: t  ! optical depth between elements j and i

  dx = slab%dx
  if ( i == j ) then
    p = 1.d0 - ( 1.d0 - 2.d0*En(3,dx) )/(2.d0*dx)
  else
    t = (abs(j-i)-1)*dx
    p = 0.5d0/dx*( En(3,t) - 2.d0*En(3,t+dx) + En(3,t+2.d0*dx) )
  endif

end function collision_probability

end module ColProbSlab
