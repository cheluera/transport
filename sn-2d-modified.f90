Program twoDSn
!---------------------------------------------------------!
! Program to solve one-group, isotropic, 2D S_n equations !
!
!
!
!---------------------------------------------------------!


implicit none

double precision, allocatable, dimension(:,:,:) :: psi, sSource
double precision, allocatable, dimension(:,:)   :: phi_new, phi_old
double precision, allocatable, dimension(:,:,:) :: left, right, top, bottom
double precision :: error, conv, ss, pi, wsum, s0

!Quadrature (x,y,z,w) - (xi,eta,mu,w)
double precision, allocatable, dimension(:) :: xi, eta, mu, w 

double precision sigmaA, sigmaS, sigmaT, c

!Mesh parameters
double precision :: Lx, Ly
double precision :: dx, dy
integer :: nx, ny

!Integers
integer :: i, j, k, l, m, n, iter, nOrds
integer :: n1, n2, n3, n4, n5, n6, n7, n8

!Maximum number of iterations
integer :: maxIter

!===========================================================!
!===========================================================!

error   = 10.0
conv    = 1.0e-6
iter    = 0
maxIter = 500

pi = acos(-1.0)

!Setup mesh
Lx = 40.0
Ly = 40.0

nx = 800
ny = 800

dx = Lx/(1.0*nx)
dy = Ly/(1.0*ny)

!Cross section information -- at this time only isotropic
sigmaA = 0.5
sigmaS = 0.5
sigmaT = sigmaA + sigmaS
c      = sigmaS/sigmaT

write(*,*) 'c = ', c
write(*,*) 'sigmaT*mesh size: ', sigmaT * dx, sigmaT * dy

!Read quadrature set
open(unit=6,file='s10-ords.txt',status='old')
read(6,*) n1, n2, n3, n4
close(6)

!Number of discrete ordinates
nOrds = n1+n2+n3+n4
write(*,*) 'nOrds=', nOrds



!Allocate memory
allocate(    psi(nx,ny,nOrds))
allocate(sSource(nx,ny,nOrds))
allocate(phi_new(nx,ny))
allocate(phi_old(nx,ny))
allocate( left(0:nx+1,0:ny+1,nOrds), right(0:nx+1,0:ny+1,nOrds), &
           top(0:nx+1,0:ny+1,nOrds),bottom(0:nx+1,0:ny+1,nOrds))

allocate( xi(nOrds))
allocate(eta(nOrds))
allocate( mu(nOrds))
allocate(  w(nOrds))

wsum = 0.0
open(unit=6,file='s10.txt',status='old')
do i=1,nOrds
  read(6,*) xi(i), eta(i), mu(i), w(i)
  wsum = wsum + w(i)
end do
close(6)
w = w*(4.0*pi/wsum)
wsum = 0.0
wsum = SUM(w)
write(*,*) 'Sum of weights: ', wsum


!Set boundary angular flux values
do j=1,ny  
!
!Left face
  do n=1,n1
    left(0,j,n) = 0.0  ! 1.0
  end do
!
  do n=n1+n2+n3+1, nOrds
    left(0,j,n) = 0.0  ! 1.0
  end do
!
!Right face
  do n=n1+1, n1+n2+n3
    right(nx+1,j,n) = 0.0
  end do
!
end do

do i=1,nx  
!
!bottom face
  do n=1,n1+n2
    bottom(i,0,n) = 0.0
  end do
!
!Top face
  do n=n1+n2+1, n1+n2+n3
    top(i,ny+1,n) = 0.0
  end do
!
end do


!Initialize psi and phi
psi     = 0.0
phi_old = 0.0

s0 = 0.0


!Start source iteration
do while((iter .lt. maxIter).and.(error .gt. conv))

  !-------Start sweeps in 1st octant: xi > 0, eta > 0, mu > 0
  do n=1,n1

    do j=1,ny  !--loop over y starting at y=0

      do i=1,nx  !--loop over x starting at x=0 

        !--Get previous scattering source + external source at (i,j,n)
	call source(i,j,nx,ny,dx,dy,pi,c,phi_new(i,j),ss)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0*xi(n)/dx)*left(i-1,j,n) +  (2.0*eta(n)/dy)*bottom(i,j-1,n) + ss) /  &
                      ( 2.0*xi(n)/dx + 2.0*eta(n)/dy + 1.0)

        !--Update left face
        left(i,j,n) = 2*psi(i,j,n) - left(i-1,j,n)
  
      end do

      !--Update bottom face
      do i=1,nx
        bottom(i,j,n) = 2*psi(i,j,n) - bottom(i,j-1,n)
      end do

    end do

  end do
  !
  !-------Start sweeps in 2nd octant: xi < 0, eta > 0, mu > 0
  !
  do n=n1+1,n1+n2

    do j=1,ny       !--loop over y starting at y=0

      do i=nx,1,-1  !--loop over x starting at x=Lx 

        !--Get previous scattering source + external source at (i,j,n)
	call source(i,j,nx,ny,dx,dy,pi,c,phi_new(i,j),ss)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0*eta(n)/dy)*bottom(i,j-1,n) + ss) /  &
                      ( 2.0*abs(xi(n))/dx + 2.0*eta(n)/dy + 1.0)

        !--Update right face
        right(i,j,n) = 2*psi(i,j,n) - right(i+1,j,n)      

      end do

      !--Update bottom face
      do i=1,nx
        bottom(i,j,n) = 2*psi(i,j,n) - bottom(i,j-1,n)
      end do

    end do

  end do
  !
  !-------Start sweeps in 3rd octant: xi < 0, eta < 0, mu > 0
  !
  do n=n1+n2+1,n1+n2+n3

    do j=ny,1,-1   !--loop over y starting at y=Ly

      do i=nx,1,-1   !--loop over x starting at x=Lx 

        !--Get previous scattering source + external source at (i,j,n)
	call source(i,j,nx,ny,dx,dy,pi,c,phi_new(i,j),ss)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0*abs(eta(n))/dy)*top(i,j+1,n) + ss) /  &
                      ( 2.0*abs(xi(n))/dx + 2.0*abs(eta(n))/dy + 1.0)

        !--Update right face
	right(i,j,n) = 2*psi(i,j,n) - right(i+1,j,n)	

      end do

      !--Update top face
      do i=1,nx
        top(i,j,n) = 2*psi(i,j,n) - top(i,j+1,n) 
      end do

    end do

  end do
  !
  !-------Start sweeps in 4th octant: xi > 0, eta < 0, mu > 0
  !
  do n=n1+n2+n3+1,n1+n2+n3+n4

    do j=ny,1,-1  !--loop over y starting at y=Ly

      do i=1,nx    !--loop over x starting at x=0 

        !--Get previous scattering source + external source at (i,j,n)
	call source(i,j,nx,ny,dx,dy,pi,c,phi_new(i,j),ss)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0*abs(xi(n))/dx)*left(i-1,j,n) +  (2.0*abs(eta(n))/dy)*top(i,j+1,n) + ss) /  &
                      ( 2.0*abs(xi(n))/dx + 2.0*abs(eta(n))/dy + 1.0)

        !--Update right face
        left(i,j,n) = 2*psi(i,j,n) - left(i-1,j,n)      
        
      end do

      !--Update top face
      do i=1,nx
        top(i,j,n) = 2*psi(i,j,n) - top(i,j+1,n)    
      end do

    end do

  end do



  !--Calculate new estimate for phi
  do i=1,nx

    do j=1,ny
      
        phi_new(i,j) =  DOT_PRODUCT(w, psi(i,j,1:nOrds))
      
    end do

  end do

  error = MAXVAL(phi_new - phi_old)
 
  write(*,*) 'Error: ',MAXVAL(phi_new - phi_old)

  phi_old = phi_new


  iter = iter + 1

end do

!Calculate total absorption rate -- should be unity for normalized source
write(*,*) 'Total absorption: ', (1.0 - c) * SUM(phi_old)*dx*dy


open(unit=6,file='test.txt',status='unknown')
do i=1,nx
  write(6,*) (phi_new(i,j),j=1,ny)
end do
close(6)




deallocate(psi)
deallocate(sSource)
deallocate(phi_new, phi_old)

deallocate(xi,eta,mu,w)
deallocate(left,right,top,bottom)

end program twoDSn

!--Subroutine for calculating scattering source and external source
!--Needs modification to allow for anisotropic sources

subroutine source(i,j,nx,ny,dx,dy,pi,c,phi,ss)

implicit none

double precision :: phi
double precision :: ss, pi, c

!Mesh parameters
double precision :: dx, dy
integer :: nx, ny

!Integers
integer :: i, j

        ss = c * phi/(4.0*pi) !sSource(i,j,n)
        if ((i .eq. nx/2) .and. (j .eq. ny/2)) then
          ss = ss + 1.0/(dx*dy*4.0*pi)
        end if
return
end subroutine source

!--Subroutines below are not yet used...

subroutine bottom(x,xi,eta,value)
  !--Boundary data on bottom face--!

  if(eta .lt. 0.0) then
    write(*,*) "Can not specify outgoing bottom flux"
    stop
  else
    value = 0.0
  end if

end 


subroutine top(x,xi,eta,value)
  !--Boundary data on top face--!

  if(xi .gt. 0.0) then
    write(*,*) "Can not specify outgoing top flux"
    stop
  else
    value = 0.0
  end if

end 

subroutine left(y,xi,eta,value)
  !--Boundary data on left face--!

  if(xi .lt. 0.0) then
    write(*,*) "Can not specify outgoing left flux"
    stop
  else
    value = 0.0
  end if

end 


subroutine right(y,xi,eta,value)
  !--Boundary data on right face--!

  if(xi .gt. 0.0) then
    write(*,*) "Can not specify outgoing right flux"
    stop
  else
    value = 0.0
  end if

end 


