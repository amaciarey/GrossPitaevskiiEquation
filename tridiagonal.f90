!-----------------------------------------------------------------------
! SUBROUTINE: TRIDIAGONAL.F90
!-----------------------------------------------------------------------
!
! This rourine solves a linear system of equations given by:
!
!                        A*y = f
!
! by performing the LU decomposition of the matrix A. The matrix A must be 
! tridiagonal, which allows us to use the Thomas algorithm.
!
! The arguments of the routine are the following:
!
!    > N : size of the matrix (A(N,N))
!    > a : Vector containing the elements of the SUB-diagonal of the matrix
!    > b : Vector containing the elements of the diagonal of the matrix
!    > c : Vector containing the elements of the SUPER-diagonal of the matrix
!    > f : Vector containing the independent terms f
!    > y : Vector containing the N solutions of the system
!
!-------------------------------------------------------------------------------
!
! Adrián Macía (2013)
!
!-------------------------------------------------------------------------------

subroutine tridiagonal(N,a,b,c,f,y)

implicit none

!Definición de variables

integer (kind=4)            :: i,j,N
complex (kind=8), dimension(N) :: a,b,c,f,y,alpha,beta,z

beta(1)  = c(1)/b(1)
alpha(1) = b(1)
z(1)     = f(1)/alpha(1)

do i=2,N
   alpha(i) = b(i)-a(i)*beta(i-1)
   beta(i)  = c(i)/alpha(i)
   z(i)     = (f(i)-a(i)*z(i-1))/alpha(i)
end do

y(N) = z(N)

do i=1,N-1
   j    = N-i
   y(j) = z(j)-beta(j)*y(j+1)
end do

return
end subroutine tridiagonal

!-----------------------------------------------------------------------

