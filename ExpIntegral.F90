module ExpIntegral
  implicit none

  private
  public  :: ExpIntN

CONTAINS

real(8) function ExpIntN ( n, x )  result(en)
  implicit none

  integer, intent(in) :: n
  real(8), intent(in) :: x

  integer :: i

  if ( x /= 0.d0 ) then
    en = ExpInt1(x);
    do i=2,n
      en = ( exp(-x) - x*en )/(real(i,8)-1.d0);
    enddo
  else
     en = 1.d0/(real(n,8)-1.d0);
  endif

end function ExpIntN


real(8) function ExpInt1 ( x )  result(e1)

!*****************************************************************************80
!
!! E1XA computes the exponential integral E1(x).
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    06 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Output, real ( kind = 8 ) E1, the function value.
!
  implicit none

  real ( kind = 8 ) es1
  real ( kind = 8 ) es2
  real ( kind = 8 ) x

  if ( x == 0.0D+00 ) then

    e1 = 1.0D+300

  else if ( x <= 1.0D+00 ) then

    e1 = - log ( x ) + (((( &
        1.07857D-03 * x &
      - 9.76004D-03 ) * x &
      + 5.519968D-02 ) * x &
      - 0.24991055D+00 ) * x &
      + 0.99999193D+00 ) * x &
      - 0.57721566D+00

  else

    es1 = ((( x &
      + 8.5733287401D+00 ) * x &
      +18.059016973D+00  ) * x &
      + 8.6347608925D+00 ) * x &
      + 0.2677737343D+00

    es2 = ((( x &
      +  9.5733223454D+00 ) * x &
      + 25.6329561486D+00 ) * x &
      + 21.0996530827D+00 ) * x &
      +  3.9584969228D+00

    e1 = exp ( - x ) / x * es1 / es2

  end if

  return
end function ExpInt1

end module ExpIntegral
