program gpe3d

use omp_lib
use, intrinsic :: iso_c_binding

implicit none

!include 'fftw3.f'

real (kind=8)    :: xmax,ymax,zmax
real (kind=8)    :: pi,g,ScatteringLength
real (kind=8)    :: dx,dy,dz,dt
real (kind=8)    :: alpha,beta,gamma,W0,WR
real (kind=8)    :: Norm,NormSq
real (kind=8)    :: tmp
real (kind=8)    :: diff,tol
real (kind=8)    :: simpson
real (kind=8)    :: ax1,ax2,ax3
real (kind=8)    :: ay1,ay2,ay3
real (kind=8)    :: az1,az2,az3
real (kind=8)    :: bx1,bx3
real (kind=8)    :: by1,by3
real (kind=8)    :: begin,end
real (kind=8)    :: Chem,En
complex (kind=8) :: iii
integer (kind=4) :: Np
integer (kind=4) :: Nx,Ny,Nz
integer (kind=4) :: i,j,k
integer (kind=4) :: is,ie,js,je,ks,ke
integer (kind=4) :: nthread,ithread,nrem,nchunk
integer (kind=4) :: it,Nit

real (kind=8), dimension (:,:,:), allocatable :: Pot
real (kind=8), dimension (:,:), allocatable   :: Ayz,Axy,Lx,Ly,Lz,L1x,L1y
real (kind=8), dimension (:), allocatable     :: Az
real (kind=8), dimension (:), allocatable     :: x,y,z

complex (kind=8), dimension (:,:,:), allocatable :: Psi,Psi_prev,Psi_aux
complex (kind=8), dimension (:), allocatable     :: bx,by,bz

!open (unit=1,file='3dgpe.in')

! Reading input parameters

read (*,*) xmax,ymax,zmax
read (*,*) dt,Nit
read (*,*) Nx,Ny,Nz
read (*,*) ScatteringLength
read (*,*) Np
read (*,*) W0,WR
read (*,*) alpha,beta,gamma

!close (unit=1)

pi  = acos(-1.d0)
tol = 1.d-6

dx = 2.d0*xmax/real(Nx)
dy = 2.d0*ymax/real(Ny)
dz = 2.d0*zmax/real(Nz)

g   = 4.d0*pi*Np*ScatteringLength
iii = (0.d0,1.d0) 

allocate (Psi(0:Nx,0:Ny,0:Nz),Psi_prev(0:Nx,0:Ny,0:Nz),Psi_aux(0:Nx,0:Ny,0:Nz))
allocate (Pot(0:Nx,0:Ny,0:Nz))
allocate (x(0:Nx),y(0:Ny),z(0:Nz))

! Initial guess for the condensate wave function and initialization of 
! the array potential

do k=0,Nz
   do j=0,Ny
      do i=0,Nx
         
         x(i) = -xmax+i*dx
         y(j) = -ymax+j*dy
         z(k) = -zmax+k*dz
         
         tmp        = (alpha*x(i))**2+(beta*y(j))**2+(gamma*z(k))**2
         Pot(i,j,k) = 0.5d0*(tmp)
         if (i==0 .or. i==Nx) then
            Psi(i,j,k) = (0.d0,0.d0)
         else if (j==0 .or. j==Ny) then
            Psi(i,j,k) = (0.d0,0.d0)
         else if (k==0 .or. k==Nz) then
            Psi(i,j,k) = (0.d0,0.d0)
         else
            Psi(i,j,k) = (1.d0-WR)*exp(-0.5d0*tmp)+WR*(x(i)+iii*y(j))*exp(-0.5d0*tmp)
            !Psi(i,j,k) = exp(-0.5d0*tmp)
         end if
         Psi_prev(i,j,k) = Psi(i,j,k)

      end do
   end do
end do

! Normalization of the wave function

call Integrate3D(Nx,Ny,Nz,xmax,ymax,zmax,abs(Psi)**2,Norm)

Psi      = Psi/sqrt(Norm)
Psi_prev = Psi

do i=0,Nx
   write (94,*) x(i),real(Psi(i,Ny/2,Nz/2)),aimag(Psi(i,Ny/2,Nz/2))
end do

do j=0,Ny
   write (95,*) y(j),real(Psi(Nx/2,j,Nz/2)),aimag(Psi(Nx/2,j,Nz/2))
end do

do k=0,Nz
   write (96,*) z(k),real(Psi(Nx/2,Ny/2,k)),aimag(Psi(Nx/2,Ny/2,k))
end do

allocate (Axy(0:Nx,0:Ny))

do j=0,Ny
   do i=0,Nx
      Axy(i,j) = simpson(-zmax,zmax,Nz,abs(Psi(i,j,:))**2)
      write (101,*) i,j,Axy(i,j)
   end do
   write (101,*)
end do

! Definition of the matrices that corresponds to the differential 
! operators

allocate (Lx(Nx-1,-1:1),Ly(Ny-1,-1:1),Lz(Nz-1,-1:1))
allocate (L1x(Nx-1,-1:1),L1y(Ny-1,-1:1))

do i=1,Nx-1
   
   if (i==1) then
      
      Lx(i,-1) = 0.d0
      Lx(i,0)  = -1.d0/dt-1.d0/(dx*dx)
      Lx(i,1)  = 0.5d0/(dx*dx)

      L1x(i,-1) = 0.d0
      L1x(i,0)  = 0.d0
      L1x(i,1)  = -0.5d0*WR/dx

   else if (i==Nx-1) then

      Lx(i,-1) = 0.5d0/(dx*dx)
      Lx(i,0)  = -1.d0/dt-1.d0/(dx*dx)
      Lx(i,1)  = 0.d0

      L1x(i,-1) = 0.5d0*WR/dx
      L1x(i,0)  = 0.d0
      L1x(i,1)  = 0.d0

   else

      Lx(i,-1) = 0.5d0/(dx*dx)
      Lx(i,0)  = -1.d0/dt-1.d0/(dx*dx)
      Lx(i,1)  = 0.5d0/(dx*dx)

      L1x(i,-1) = 0.5d0*WR/dx
      L1x(i,0)  = 0.d0
      L1x(i,1)  = -0.5d0*WR/dx
      
   end if

end do

do j=1,Ny-1
   
   if (j==1) then
      
      Ly(j,-1) = 0.d0
      Ly(j,0)  = -1.d0/dt-1.d0/(dy*dy)
      Ly(j,1)  = 0.5d0/(dy*dy)

      L1y(j,-1) = 0.d0
      L1y(j,0)  = 0.d0
      L1y(j,1)  = 0.5d0*WR/dy

   else if (j==Ny-1) then

      Ly(j,-1) = 0.5d0/(dy*dy)
      Ly(j,0)  = -1.d0/dt-1.d0/(dy*dy)
      Ly(j,1)  = 0.d0
      
      L1y(j,-1) = -0.5d0*WR/dy
      L1y(j,0)  = 0.d0
      L1y(j,1)  = 0.d0

   else

      Ly(j,-1) = 0.5d0/(dy*dy)
      Ly(j,0)  = -1.d0/dt-1.d0/(dy*dy)
      Ly(j,1)  = 0.5d0/(dy*dy)

      L1y(j,-1) = -0.5d0*WR/dy
      L1y(j,0)  = 0.d0
      L1y(j,1)  = 0.5d0*WR/dy

   end if
   
end do

do k=1,Nz-1
   
   if (k==1) then
      
      Lz(k,-1) = 0.d0
      Lz(k,0)  = -1.d0/dt-1.d0/(dz*dz)
      Lz(k,1)  = 0.5d0/(dz*dz)

   else if (k==Nx-1) then

      Lz(k,-1) = 0.5d0/(dz*dz)
      Lz(k,0)  = -1.d0/dt-1.d0/(dz*dz)
      Lz(k,1)  = 0.d0

   else

      Lz(k,-1) = 0.5d0/(dz*dz)
      Lz(k,0)  = -1.d0/dt-1.d0/(dz*dz)
      Lz(k,1)  = 0.5d0/(dz*dz)

   end if

end do

! Begin of the time evolution loop.
! At each time step we apply the following decomposition of the time 
! evolution operator:
!
! U(t,0) = e^(-0.5*Veff*dt) e^(-D*dt) e^(-0.5*Veff*dt)
!
! where 
!
! Veff = V(x,y,z)+g|Psi(x,y,z)|^2 
!
! and
! 
! D = -0.5*(Laplacian)
! 
! This decomposition of the time evolution operator has O(dt^3)
! accuracy 

!allocate (bx(1:Nx-1),by(1:Ny-1),bz(1:Nz-1))

ax1 = -0.25d0/(dx*dx)
ax2 = -1.d0/dt+0.5d0/(dx*dx)
ax3 = ax1

ay1 = -0.25d0/(dy*dy)
ay2 = -1.d0/dt+0.5d0/(dy*dy)
ay3 = ay1

az1 = -0.25d0/(dz*dz)
az2 = -1.d0/dt+0.5d0/(dz*dz)
az3 = az1

bx1 = -0.25d0*WR/dx
bx3 = -bx1

by1 = 0.25d0*WR/dy
by3 = -by1

allocate (Az(0:Nz),Ayz(0:Ny,0:Nz))

diff = 1.d0
it   = 0

!do while (diff/dt >= tol) 
do it=1,Nit

   call cpu_time(begin)
   diff = 0.d0

   !it = it+1

   !$OMP PARALLEL&
   !$OMP& PRIVATE(ithread,nrem,nchunk,ks,ke,js,je,is,ie,tmp,bx,by,bz)&
   !$OMP& REDUCTION(max:diff)

   allocate (bx(1:Nx-1),by(1:Ny-1),bz(1:Nz-1))

   nthread = omp_get_num_threads()   
   ithread = omp_get_thread_num()
   nrem    = mod(Nz-1,nthread)
   nchunk  = (Nz-1-nrem)/nthread 

   if (ithread < nrem) then
      ks = 1+ithread*(nchunk+1)
      ke = ks+nchunk
   else
      ks = 1+ithread*nchunk+nrem
      ke = ks+nchunk
   end if

   ! First step of the propagation

   do k=ks,ke
      do j=1,Ny-1
         do i=1,Nx-1

            tmp            = Pot(i,j,k)+g*abs(Psi_prev(i,j,k))**2
            Psi(i,j,k)     = exp(-tmp*dt)*Psi_prev(i,j,k)
            Psi_aux(i,j,k) = Psi(i,j,k) 
            
         end do
      end do
   end do

   !$OMP BARRIER

   ! Start the propagation of the linear term
   ! Define the independent term of the linear system to be solved
   ! Solving the propagation associated with the x derivatives

   do k=ks,ke
      do j=1,Ny-1
         do i=1,Nx-1
            bx(i) = -Psi_aux(i,j,k)/dt
         end do
         call tridiagonal(Nx-1,Lx(:,-1)+iii*y(j)*L1x(:,-1),&
                         &Lx(:,0)+iii*0.d0,Lx(:,1)+iii*y(j)*L1x(:,1),&
                         &bx,Psi(1:Nx-1,j,k))

         Psi_aux(1:Nx-1,j,k) = Psi(1:Nx-1,j,k)
      end do
   end do

   !$OMP BARRIER
   
   ! Solving the propagation associated with the y derivatives

   do k=ks,ke
      do i=1,Nx-1
         do j=1,Ny-1
            by(j) = -Psi_aux(i,j,k)/dt
         end do
         call tridiagonal(Ny-1,Ly(:,-1)+iii*x(i)*L1y(:,-1),&
                         &Ly(:,0)+iii*0.d0,Ly(:,1)+iii*x(i)*L1y(:,1),&
                         &by,Psi(i,1:Ny-1,k))
         Psi_aux(i,1:Ny-1,k) = Psi(i,1:Ny-1,k)
      end do
   end do
   
   !$OMP BARRIER

   ! Solving the propagation associated with the z derivatives

   nrem   = mod(Ny-1,nthread)
   nchunk = (Ny-1-nrem)/nthread 
   
   if (ithread < nrem) then
      js = 1+ithread*(nchunk+1)
      je = js+nchunk
   else
      js = 1+ithread*nchunk+nrem
      je = js+nchunk
   end if

   do j=js,je
      do i=1,Nx-1
         do k=1,Nz-1
            bz(k) = -Psi_aux(i,j,k)/dt
         end do
         call tridiagonal(Nz-1,Lz(:,-1)+iii*0.d0,Lz(:,0)+iii*0.d0,Lz(:,1)+iii*0.d0,bz,Psi(i,j,1:Nz-1))
         Psi_aux(i,j,1:Nz-1) = Psi(i,j,1:Nz-1)
      end do
   end do

   !$OMP BARRIER

   nrem   = mod(Nz,nthread)
   nchunk = (Nz-nrem)/nthread 
   
   if (ithread < nrem) then
      ks = 1+ithread*(nchunk+1)
      ke = ks+nchunk
   else
      ks = 1+ithread*nchunk+nrem
      ke = ks+nchunk-1
   end if

   do k=ks,ke
      do j=0,Ny
         Ayz(j,k) = simpson(-xmax,xmax,Nx,abs(Psi(:,j,k))**2)
      end do
      Az(k) = simpson(-ymax,ymax,Ny,Ayz(:,k))
   end do
   !$OMP BARRIER

   !$OMP SINGLE
   Norm   = simpson(-zmax,zmax,Nz,Az)
   NormSq = sqrt(Norm)
   !$OMP END SINGLE

   do k=ks,ke
      do j=0,Ny
         do i=0,Nx
            Psi(i,j,k)      = Psi(i,j,k)/NormSq
            Psi_aux(i,j,k)  = Psi(i,j,k)
            diff            = max(diff,abs(Psi(i,j,k)-Psi_prev(i,j,k))/dt)
            Psi_prev(i,j,k) = Psi(i,j,k)
         end do
      end do
   end do

   deallocate (bx,by,bz)

   !$OMP END PARALLEL

   call cpu_time(end)

   if (mod(it,1)==0) then
      print *, "Time elapsed in iteration",it,": ",end-begin,diff
   end if

end do

deallocate (Psi_prev,Psi_aux,Lx,Ly,Lz,L1x,L1y,Ayz,Az)

! Evaluate chemical potential and Energy

call Chem_and_En(Nx,Ny,Nz,xmax,ymax,zmax,g,Pot,Psi,Chem,En)

print *, '=============================================================='
print *, '                     FINAL RESULTS:                           '
print *, '=============================================================='
print *, ''
print *, ' # Time steps to achieve convergence :',it
print *, ' # Chemical potential                :',Chem
print *, ' # Total energy per particle         :',En
print *, ' '

do i=0,Nx
   write (97,*) x(i),abs(Psi(i,Ny/2,Nz/2))
end do

do j=0,Ny
   write (98,*) y(j),abs(Psi(Nx/2,j,Nz/2))
end do

do k=0,Nz
   write (99,*) z(k),abs(Psi(Nx/2,Ny/2,k))
end do

do j=0,Ny
   do i=0,Nx
      Axy(i,j) = simpson(-zmax,zmax,Nz,abs(Psi(i,j,:))**2)
      write (100,*) i,j,Axy(i,j)
   end do
   write (100,*)
end do

deallocate (x,y,z,Pot,Psi,Axy)

end program gpe3d

!=======================================================================

subroutine Chem_and_En(Nx,Ny,Nz,xmax,ymax,zmax,g,Pot,Psi,Chem,En)

implicit none

real (kind=8)    :: g
real (kind=8)    :: Chem,En
real (kind=8)    :: dx,dy,dz
real (kind=8)    :: xmax,ymax,zmax
integer (kind=4) :: Nx,Ny,Nz
integer (kind=4) :: i,j,k

real (kind=8), dimension (0:Nx,0:Ny,0:Nz)    :: Pot
complex (kind=8), dimension (0:Nx,0:Ny,0:Nz) :: DPsiX,DPsiY,DPsiZ
complex (kind=8), dimension (0:Nx,0:Ny,0:Nz) :: Psi


dx = 2.d0*xmax/real(Nx)
dy = 2.d0*ymax/real(Ny)
dz = 2.d0*zmax/real(Nz)

! Initialize the boundaries

DPsiX(0,:,:) = 0.d0
DPsiX(:,0,:) = 0.d0
DPsiX(:,:,0) = 0.d0

DPsiX(Nx,:,:) = 0.d0
DPsiX(:,Ny,:) = 0.d0
DPsiX(:,:,Nz) = 0.d0

DPsiY(0,:,:) = 0.d0
DPsiY(:,0,:) = 0.d0
DPsiY(:,:,0) = 0.d0

DPsiY(Nx,:,:) = 0.d0
DPsiY(:,Ny,:) = 0.d0
DPsiY(:,:,Nz) = 0.d0

DPsiZ(0,:,:) = 0.d0
DPsiZ(:,0,:) = 0.d0
DPsiZ(:,:,0) = 0.d0

DPsiZ(Nx,:,:) = 0.d0
DPsiZ(:,Ny,:) = 0.d0
DPsiZ(:,:,Nz) = 0.d0

do k=1,Nz-1
   do j=1,Ny-1
      do i=1,Nx-1
         DPsiX(i,j,k) = 0.5d0*(Psi(i+1,j,k)-Psi(i-1,j,k))/dx
         DPsiY(i,j,k) = 0.5d0*(Psi(i,j+1,k)-Psi(i,j-1,k))/dy
         DPsiZ(i,j,k) = 0.5d0*(Psi(i,j,k+1)-Psi(i,j,k-1))/dz
      end do
   end do
end do

call Integrate3D(Nx,Ny,Nz,xmax,ymax,zmax,0.5d0*(abs(DPsiX)**2+&
                &abs(DPsiY)**2+abs(DPsiZ)**2)+Pot*abs(Psi)**2+&
                &0.5d0*g*abs(Psi)**4,En)
call Integrate3D(Nx,Ny,Nz,xmax,ymax,zmax,0.5d0*(abs(DPsiX)**2+&
                &abs(DPsiY)**2+abs(DPsiZ)**2)+Pot*abs(Psi)**2+&
                &g*abs(Psi)**4,Chem)


return
end subroutine Chem_and_En

!-----------------------------------------------------------------------

subroutine Integrate3D(Nx,Ny,Nz,xmax,ymax,zmax,Psi,Atot)

implicit none

real (kind=8)    :: xmax,ymax,zmax
real (kind=8)    :: Atot
real (kind=8)    :: simpson
integer (kind=4) :: Nx,Ny,Nz
integer (kind=4) :: j,k

real (kind=8), dimension (0:Nx,0:Ny,0:Nz) :: Psi
real (kind=8), dimension (0:Ny,0:Nz)      :: Ayz
real (kind=8), dimension (0:Nz)           :: Az

do k=0,Nz
   do j=0,Ny
      Ayz(j,k) = simpson(-xmax,xmax,Nx,Psi(:,j,k))
   end do
   Az(k) = simpson(-ymax,ymax,Ny,Ayz(:,k))
end do

Atot = simpson(-zmax,zmax,Nz,Az)

return
end subroutine Integrate3D

         
