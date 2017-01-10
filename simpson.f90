!-----------------------------------------------------------------------
! SUBROUTINE: SIMPSON.F90
!-----------------------------------------------------------------------
!
! This subroutine performs the numerical integration of a function f(x) 
! in the  interval [a,b] by using the Simpson rule. 
! The arguments of the subroutine are the following:
!
!     > a: Lower limit of the integration region
!     > b: Upper limit of the integration region
!     > N: Number of points in the integration grid
!     > f: Function to be integrated, it is an N-element array containing 
!          the value of the function in the grid points
!
!------------------------------------------------------------------------
!
! Adrián Macía (2013)
!
!------------------------------------------------------------------------

real (kind=8) function simpson(a,b,N,f)

implicit none

real (kind=8)                 :: a,b,dx
integer (kind=4)              :: j,N
real (kind=8), dimension(0:N) :: f

dx = (b-a)/real(N)

simpson = 0.d0

do j=1,N/2
   simpson = simpson+4.d0*f(2*j-1)
end do

do j=1,N/2-1
   simpson = simpson+2.d0*f(2*j)
end do

simpson = dx*(simpson+f(0)+f(N))/3.d0

return
end function simpson

!-----------------------------------------------------------------------
