Program twoDSn

!---------------------------------------------------------!
! Program to solve one-group, isotropic, 2D S_n equations !
!
!
!
!---------------------------------------------------------!

use mpi

implicit none

double precision, allocatable, dimension(:,:,:) :: psi, sSource, ss, ext_ss
double precision, allocatable, dimension(:,:)   :: phi_new, phi_mfd
double precision, allocatable, dimension(:,:,:) :: left, right, top, bottom
double precision, allocatable, dimension(:)     :: legP_coefs
double precision, allocatable, dimension(:,:,:)   :: psi_mfd

double precision :: error, conv, t_ss, pi, wsum

!Quadrature (x,y,z,w) - (xi,eta,mu,w)
double precision, allocatable, dimension(:) :: xi, eta, mu, w 
double precision                            :: cosphi, sinphi

double precision sigmaA, sigmaS, sigmaT, c

!Spherical harmonics
double precision, allocatable, dimension(:,:) :: p, phir, phii
double precision, allocatable, dimension(:)   :: phire, phiim
double precision, allocatable, dimension(:,:,:) :: sphre

!Angular flux moments
double precision, allocatable, dimension(:,:,:,:) :: fmre, fmim, fmre_old
double precision, allocatable, dimension(:) :: fm

!Mesh parameters
double precision :: Lx, Ly
double precision :: dx, dy
integer :: nx, ny

!Integers
integer :: i, j, m, n, iter, nOrds, n_moments, n_max, n_min, n_sph ! k, l 
integer :: n1, n2, n3, n4 ! n5, n6, n7, n8

!Maximum number of iterations
integer :: maxIter

!MPI variables
integer :: ierror, rank, num_cores, master

!===========================================================!
!===========================================================!

error   = 10.0d0
conv    = 1.0d-15
iter    = 0
maxIter = 5000

pi = 3.1415926535897932385d0

!Setup mesh
Lx = 20.0
Ly = 20.0

nx = 20
ny = 20

dx = Lx/(1.0*nx)
dy = Ly/(1.0*ny)

!Cross section information
sigmaA = 0.50d0
sigmaS = 0.50d0
sigmaT = sigmaA + sigmaS
c      = sigmaS/sigmaT

n_moments = 0 !max moments for mfd solution

!Legendre scattering coefficients
n_max = 0
allocate(legP_coefs(0:n_max))
  do i=0,n_max
  legP_coefs(i) = 1.0d0 / (2.0d0**i)
  end do
  
write(*,*) 'c = ', c
write(*,*) 'sigmaT*mesh size: ', sigmaT * dx, sigmaT * dy

!Read quadrature set
open(unit=7,file='../../quadratures/level-symmetric/s4-ords.txt',status='old')
read(7,*) n1, n2, n3, n4
close(7)

!Number of discrete ordinates
nOrds = 2 * (n1+n2+n3+n4)
write(*,*) 'nOrds=', nOrds

!Allocate memory
allocate(    psi(nx,ny,nOrds))
allocate(sSource(nx,ny,nOrds))
allocate(phi_new(nx,ny))

allocate( left(0:nx+1,0:ny+1,nOrds), right(0:nx+1,0:ny+1,nOrds), &
           top(0:nx+1,0:ny+1,nOrds),bottom(0:nx+1,0:ny+1,nOrds))
allocate(    ss(nx,ny,nOrds))
allocate(ext_ss(nx,ny,nOrds))
allocate(psi_mfd(nx,ny,nOrds))
allocate(phi_mfd(nx,ny))

!Space for quadrature
allocate( xi(nOrds))
allocate(eta(nOrds))
allocate( mu(nOrds))
allocate(  w(nOrds))

n_sph = MAX(n_moments, n_max)
n_min = MIN(n_moments, n_max)

!Space for spherical harmonics
allocate(p(0:n_sph,0:n_sph))
allocate(phire(0:n_sph), phiim(0:n_sph))
allocate(sphre(nOrds,0:n_sph,0:n_sph))
allocate(phir(nOrds,0:n_sph))
allocate(phii(nOrds,0:n_sph))

!Angular flux moments
allocate(fm(0:n_max))
allocate(fmre(nx,ny,0:n_max,0:n_max))
allocate(fmim(nx,ny,0:n_max,0:n_max))

!Read quadrature
wsum = 0.0d0
open(unit=7,file='../../quadratures/level-symmetric/s4.txt',status='old')
do i=1,nOrds
  read(7,*) xi(i), eta(i), mu(i), w(i)
  wsum = wsum + w(i)
end do
close(7)
w = w*(4.0d0*pi/wsum)
wsum = 0.0d0
wsum = SUM(w)
write(*,*) 'Sum of weights divided by 4pi: ', wsum/(4.0d0*pi)

!Initialize psi and phi
psi     = 0.0

ss = 0.0

!Precompute spherical harmonics at the quadrature points

do j=1,nOrds

  cosphi =  xi(j) / (1.0d0 - mu(j)*mu(j) ) ** 0.5d0

  sinphi = eta(j) / (1.0d0 - mu(j)*mu(j) ) ** 0.5d0

  call spharm(n_sph,mu(j),sinphi,cosphi,p,phire,phiim)

  !populate a 3-index array with p, phire and phiim
  do n=0,n_sph

    do m=0,n
      sphre(j,m,n) = p(m,n) 
     
    phir(j,m) = phire(m)                
    phii(j,m) = phiim(m)
    
    end do                   
  end do

end do

!--Call mfd solution, bondary conditions, and external sources
call mfd(nx,ny,dx,dy,xi,eta,nOrds,n_moments,n_max,n_min,n_sph,legP_coefs,c,sphre,phir,phii,left,right,&
         bottom,top,ext_ss,psi_mfd,phi_mfd)

!Begin MPI
call MPI_Init(ierror)
call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierror)
call MPI_Comm_size(MPI_COMM_WORLD,num_cores,ierror)
write(*,*) rank
!Start source iteration
do while((iter .lt. maxIter).and.(error .gt. conv))

  !-------Start sweeps in 1st octant: xi > 0, eta > 0, mu > 0
  do n=1,n1

    do j=1,ny  !--loop over y starting at y=0

      do i=1,nx  !--loop over x starting at x=0 

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)

      end do

      !--Update bottom face
      do i=1,nx
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)
      end do

    end do

  end do
  !
  !-------Start sweeps in 2nd octant: xi < 0, eta > 0, mu > 0
  !
  do n=n1+1,n1+n2

    do j=1,ny       !--loop over y starting at y=0

      do i=nx,1,-1  !--loop over x starting at x=Lx 

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)      

      end do

      !--Update bottom face
      do i=1,nx
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)
      end do

    end do

  end do
  !
  !-------Start sweeps in 3rd octant: xi < 0, eta < 0, mu > 0
  !
  do n=n1+n2+1,n1+n2+n3

    do j=ny,1,-1   !--loop over y starting at y=Ly

      do i=nx,1,-1   !--loop over x starting at x=Lx 

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)

        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)
	
      end do

      !--Update top face
      do i=1,nx
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n) 
      end do

    end do

  end do
  !
  !-------Start sweeps in 4th octant: xi > 0, eta < 0, mu > 0
  !
  do n=n1+n2+n3+1,n1+n2+n3+n4

    do j=ny,1,-1  !--loop over y starting at y=Ly

      do i=1,nx    !--loop over x starting at x=0 

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)      
        
      end do

      !--Update top face
      do i=1,nx
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n)    
      end do

    end do

  end do
  !--------------------------------------!
  !------Done with transport sweeps------!
  !--------------------------------------!  

  !--Assign angular flux for lower hemisphere using parity condition for 2D
  !--psi(x,y,mu,phi) = psi(x,y,-mu,phi)
  do i=1,nx
  
    do j=1,ny
    
      do n=nOrds/2+1,nOrds
      
      psi(i,j,n) = psi(i,j,n-nOrds/2)
  
      end do
      
    end do
    
  end do

  !--Calculate moments
  fmre_old = fmre

  call f_moments(nx,ny,psi,w,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim)

  !--scattering source via spherical harmonics
  call source(nx,ny,legP_coefs,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim,ss)

  error = abs(MAXVAL(fmre - fmre_old))
  
  !write(*,*) 'Diff: ', error

  iter = iter + 1

end do

  !--Calculate scalar flux
  phi_new = 0.0d0

  do i=1,nx

    do j=1,ny
    
      do n=1,nOrds
      
      phi_new(i,j) =  phi_new(i,j) + w(n) * psi(i,j,n)
 
      end do

    end do
    
  end do


call output(nx,ny,dx,dy,nOrds,psi,psi_mfd,phi_new,phi_mfd)

call MPI_Finalize(ierror)

deallocate(psi, ss)
deallocate(sSource, ext_ss)
deallocate(phi_new)
deallocate(legP_coefs, psi_mfd, phi_mfd)


deallocate(p, phir, phii)
deallocate(phire, phiim)
deallocate(sphre)
deallocate(fmre, fmim, fmre_old)
deallocate(fm)


deallocate(xi,eta,mu,w)
deallocate(left,right,top,bottom)

end program twoDSn
