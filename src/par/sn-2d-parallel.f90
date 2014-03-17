Program twoDSn

use mpi
use func_ang
use func_mfd
use func_out

double precision, allocatable, dimension(:,:,:) :: psi, ss, ext_ss, left, right
double precision, allocatable, dimension(:,:)   :: phi_new, phi_mfd
double precision, allocatable, dimension(:,:,:) :: top, bottom, ferror
double precision, allocatable, dimension(:,:) :: b_r, t_r, l_r, r_r, f_r, phi_tot
double precision, allocatable, dimension(:)     :: legP_coefs
double precision, allocatable, dimension(:,:,:)   :: psi_mfd, dfr, strm

double precision :: error, conv, t_ss, pi, wsum, errormax

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
double precision, allocatable, dimension(:,:) :: fm,fmi

!Mesh parameters
double precision :: Lx, Ly
double precision :: dx, dy
integer :: nx, ny

!Error parameters
double precision, allocatable, dimension(:,:,:) :: diff
double precision, allocatable, dimension(:,:) :: diff2
double precision :: infin, smax
double precision :: L_infin
double precision :: sL_infin

!Integers
integer :: i, j, m, n, iter, nOrds, n_moments, n_max, n_min, n_sph ! k, l 
integer :: n1, n2, n3, n4 ! n5, n6, n7, n8

!Maximum number of iterations
integer :: maxIter

!MPI variables
integer :: ierror, rank, num_cores, master, my_status, half, n_size
double precision, allocatable, dimension(:) :: snd_bot, snd_top
double precision :: start, finish
logical :: lower, spatial

!===========================================================!
!===========================================================!

  !Begin MPI
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierror)
  call MPI_Comm_size(MPI_COMM_WORLD,num_cores,ierror)

spatial = .true.

start = MPI_Wtime()

master = 0
half = num_cores/2
if (rank .lt. half) then
  lower = .true.
else
  lower = .false.
end if

errormax   = 10.0d0
conv    = 1.0d-14
iter    = 0
maxIter = 5000
pi = 3.1415926535897932385d0

!Setup mesh
Lx = 10.0d0
Ly = 10.0d0

nx = 200
ny = 200

dx = Lx/(1.0d0*nx)
dy = Ly/(1.0d0*ny)

if ( MOD(ny,num_cores) .ne. 0 ) then
  write(*,*) 'Number of rows must be divisible by number of cores'
  call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
end if

n_size = ny/num_cores

!Cross section information
sigmaA = 0.50d0
sigmaS = 0.50d0
sigmaT = sigmaA + sigmaS
c      = sigmaS/sigmaT

n_moments = 20 !max moments for mfd solution

!Legendre scattering coefficients
n_max = 20
allocate(legP_coefs(0:n_max))
  do i=0,n_max
  legP_coefs(i) = 1.0d0 / (2.0d0**i)
  end do

if (rank .eq. master) then  
write(*,*) 'c = ', c
write(*,*) 'sigmaT*mesh size: ', sigmaT * dx, sigmaT * dy
end if

!Read quadrature set
open(unit=7,file='../../quadratures/level-symmetric/s24-ords.txt',status='old')
read(7,*) n1, n2, n3, n4
close(7)

!Number of discrete ordinates
nOrds = 2 * (n1+n2+n3+n4)
if(rank .eq. master) then
write(*,*) 'nOrds=', nOrds
end if

!Allocate memory
allocate(    psi(nx,n_size,nOrds))
allocate(phi_new(nx,n_size))

allocate( left(0:nx,n_size,nOrds), right(nx+1,n_size,nOrds), &
           top(nx,n_size+1,nOrds),bottom(nx,0:n_size,nOrds))
allocate( l_r(0:nx,n_size), r_r(nx+1,n_size), t_r(nx,n_size+1),&
          b_r(nx,0:n_size), f_r(nx,n_size), dfr(nx,n_size,nOrds))
allocate(    ss(nx,n_size,nOrds))
allocate(ext_ss(nx,n_size,nOrds))
allocate(psi_mfd(nx,n_size,nOrds))
allocate(phi_mfd(nx,n_size))
allocate(snd_bot(nx))
allocate(snd_top(nx))
allocate(strm(nx,n_size,nOrds))

!Error
allocate(diff(nx,n_size,nOrds))
allocate(diff2(nx,n_size))
allocate(ferror(nx,ny,nOrds))
allocate(phi_tot(nx,ny))

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
allocate(fm(0:n_moments,0:n_moments))
allocate(fmi(0:n_moments,0:n_moments))
allocate(fmre(nx,n_size,0:n_max,0:n_max))
allocate(fmim(nx,n_size,0:n_max,0:n_max))
allocate(fmre_old(nx,n_size,0:n_max,0:n_max))

!Read quadrature
wsum = 0.0d0
open(unit=7,file='../../quadratures/level-symmetric/s24.txt',status='old')
do i=1,nOrds
  read(7,*) xi(i), eta(i), mu(i), w(i)
  wsum = wsum + w(i)
end do
close(7)
w = w*(4.0d0*pi/wsum)
wsum = 0.0d0
wsum = SUM(w)
if (rank .eq. master) then
write(*,*) 'Sum of weights divided by 4pi: ', wsum/(4.0d0*pi)
end if

!Initialize psi and phi
psi = 0.0d0
ss = 0.0d0
ext_ss = 0.0d0

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
call mfd(spatial,rank,master,num_cores,nx,ny,n_size,dx,dy,xi,eta,nOrds,n_moments,&
         l_r,r_r,b_r,t_r,f_r,dfr,fm,fmi,phi_mfd)

!--Calculate manufactured solution
call solution(rank,master,num_cores,nx,n_size,n_moments,n_sph,nOrds,f_r,l_r,r_r,b_r,t_r,dfr,&
              fm,fmi,sphre,phir,phii,strm,psi_mfd,left,right,top,bottom)

!--external source via manufactured solution
call ext_source(nx,n_size,legP_coefs,c,n_min,n_max,n_moments,n_sph,&
                nOrds,strm,psi_mfd,fm,fmi,f_r,sphre,phir,phii,ext_ss)

snd_bot = 0.0d0
snd_top = 0.0d0

!Start source iteration
do while((iter .lt. maxIter).and.(errormax .gt. conv))

select case (lower)

case (.true.)

  !-------Start sweeps in 1st octant: xi > 0, eta > 0, mu > 0
  do n=1,n1

    do j=1,n_size

      do i=1,nx  !--loop over x starting at x=0 

        !Receive face-centered flux from previous processor
        if (rank .ne. master .and. j .eq. 1) then
        call MPI_Recv( snd_bot, nx, MPI_DOUBLE_PRECISION, rank-1, (rank-1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        bottom(:,0,n) = snd_bot
        end if  

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)

        !--Update bottom face
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)

        !Send face-centered flux to next processor
        if (j .eq. n_size) then
        snd_bot = bottom(:,j,n)
        call MPI_Send(snd_bot, nx, MPI_DOUBLE_PRECISION, rank+1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do
  !
  !-------Start sweeps in 2nd octant: xi < 0, eta > 0, mu > 0
  !
  do n=n1+1,n1+n2

    do j=1,n_size

      do i=nx,1,-1  !--loop over x starting at x=Lx 

        !Receive face-centered flux from previous processor
        if (rank .ne. master .and. j .eq. 1) then
        call MPI_Recv(snd_bot, nx, MPI_DOUBLE_PRECISION, rank-1, (rank-1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        bottom(:,0,n) = snd_bot
        end if  

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)      

        !--Update bottom face
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)

        !Send face-centered flux to next processor
        if (j .eq. n_size) then
        snd_bot = bottom(:,j,n)
        call MPI_Send(snd_bot, nx, MPI_DOUBLE_PRECISION, rank+1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do

case(.false.)

  !
  !-------Start sweeps in 3rd octant: xi < 0, eta < 0, mu > 0
  !
  do n=n1+n2+1,n1+n2+n3

    do j=n_size,1,-1

      do i=nx,1,-1   !--loop over x starting at x=Lx 

        !Receive face-centered flux from previous processor
        if (rank .ne. num_cores-1 .and. j .eq. n_size) then
        call MPI_Recv(snd_top, nx, MPI_DOUBLE_PRECISION, rank+1, (rank+1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        top(:,n_size+1,n) = snd_top
        end if  

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)

        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)

        !--Update top face
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n) 

        !Send face-centered flux to next processor
        if (j .eq. 1) then
        snd_top = top(:,j,n)
        call MPI_Send(snd_top, nx, MPI_DOUBLE_PRECISION, rank-1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do
  !
  !-------Start sweeps in 4th octant: xi > 0, eta < 0, mu > 0
  !
  do n=n1+n2+n3+1,n1+n2+n3+n4

    do j=n_size,1,-1

      do i=1,nx    !--loop over x starting at x=0 

        !Receive face-centered flux from previous processor
        if (rank .ne. num_cores-1 .and. j .eq. n_size) then
        call MPI_Recv(snd_top, nx, MPI_DOUBLE_PRECISION, rank+1, (rank+1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        top(:,n_size+1,n) = snd_top
        end if  

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)      

        !--Update top face
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n)    

        !Send face-centered flux to next processor
        if (j .eq. 1) then
        snd_top = top(:,j,n)
        call MPI_Send(snd_top, nx, MPI_DOUBLE_PRECISION, rank-1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do

end select

select case(lower)

case (.true.)

  !
  !-------Start sweeps in 3rd octant: xi < 0, eta < 0, mu > 0
  !
  do n=n1+n2+1,n1+n2+n3

    do j=n_size,1,-1

      do i=nx,1,-1   !--loop over x starting at x=Lx 

        !Receive face-centered flux from previous processor
        if (j .eq. n_size) then
        call MPI_Recv(snd_top, nx, MPI_DOUBLE_PRECISION, rank+1, (rank+1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        top(:,n_size+1,n) = snd_top
        end if

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)
if (i .eq. 1) then
write(*,*) ' '
end if
        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)

        !--Update top face
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n) 

        !Send face-centered flux to next processor
        if (rank .ne. master .and. j .eq. 1) then
        snd_top = top(:,j,n)
        call MPI_Send(snd_top, nx, MPI_DOUBLE_PRECISION, rank-1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do
  !
  !-------Start sweeps in 4th octant: xi > 0, eta < 0, mu > 0
  !
  do n=n1+n2+n3+1,n1+n2+n3+n4

    do j=n_size,1,-1

      do i=1,nx    !--loop over x starting at x=0 

        !Receive face-centered flux from previous processor
        if (j .eq. n_size) then
        call MPI_Recv(snd_top, nx, MPI_DOUBLE_PRECISION, rank+1, (rank+1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        top(:,n_size+1,n) = snd_top
        end if

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*abs(eta(n))/dy)*top(i,j+1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*abs(eta(n))/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)      

        !--Update top face
        top(i,j,n) = 2.0d0*psi(i,j,n) - top(i,j+1,n)    

        !Send face-centered flux to next processor
        if (rank .ne. master .and. j .eq. 1) then
        snd_top = top(:,j,n)
        call MPI_Send(snd_top, nx, MPI_DOUBLE_PRECISION, rank-1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do

case(.false.)

  !-------Start sweeps in 1st octant: xi > 0, eta > 0, mu > 0
  do n=1,n1

    do j=1,n_size

      do i=1,nx  !--loop over x starting at x=0 

        !Receive face-centered flux from previous processor
        if (j .eq. 1) then
        call MPI_Recv(snd_bot, nx, MPI_DOUBLE_PRECISION, rank-1, (rank-1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        bottom(:,0,n) = snd_bot
        end if

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)

        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*xi(n)/dx)*left(i-1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*xi(n)/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update left face
        left(i,j,n) = 2.0d0*psi(i,j,n) - left(i-1,j,n)

        !--Update bottom face
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)

        !Send face-centered flux to next processor
        if (rank .ne. num_cores-1 .and. j .eq. n_size) then
        snd_bot = bottom(:,j,n)
        call MPI_Send(snd_bot, nx, MPI_DOUBLE_PRECISION, rank+1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do
  !
  !-------Start sweeps in 2nd octant: xi < 0, eta > 0, mu > 0
  !
  do n=n1+1,n1+n2

    do j=1,n_size

      do i=nx,1,-1  !--loop over x starting at x=Lx 

        !Receive face-centered flux from previous processor
        if (j .eq. 1) then
        call MPI_Recv(snd_bot, nx, MPI_DOUBLE_PRECISION, rank-1, (rank-1)*i*n, MPI_COMM_WORLD, my_status, ierror)
        bottom(:,0,n) = snd_bot
        end if

        !--total scattering source
        t_ss = ext_ss(i,j,n) + c*ss(i,j,n)
		
        !--Update cell centered quantity
        psi(i,j,n) = ( (2.0d0*abs(xi(n))/dx)*right(i+1,j,n) +  (2.0d0*eta(n)/dy)*bottom(i,j-1,n) + t_ss) /  &
                      ( 2.0d0*abs(xi(n))/dx + 2.0d0*eta(n)/dy + 1.0d0)

        !--Update right face
        right(i,j,n) = 2.0d0*psi(i,j,n) - right(i+1,j,n)      

        !--Update bottom face
        bottom(i,j,n) = 2.0d0*psi(i,j,n) - bottom(i,j-1,n)

        !Send face-centered flux to next processor
        if (rank .ne. num_cores-1 .and. j .eq. n_size) then
        snd_bot = bottom(:,j,n)
        call MPI_Send(snd_bot, nx, MPI_DOUBLE_PRECISION, rank+1, rank*i*n, MPI_COMM_WORLD, ierror)
        end if

      end do

    end do

  end do

end select
  !--------------------------------------!
  !------Done with transport sweeps------!
  !--------------------------------------!  

  !--Assign angular flux for lower hemisphere using parity condition for 2D
  !--psi(x,y,mu,phi) = psi(x,y,-mu,phi)
  do j=1,n_size

    do i=1,nx
    
      do n=nOrds/2+1,nOrds
      
      psi(i,j,n) = psi(i,j,n-nOrds/2)
  
      end do

    end do
    
  end do

  !--Calculate moments
  call f_moments(nx,n_size,psi,w,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim)

  !--scattering source via spherical harmonics
  call source(nx,n_size,legP_coefs,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim,ss)

  error = abs(MAXVAL(fmre - fmre_old))

  call MPI_Allreduce( error, errormax, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC,&
                   MPI_COMM_WORLD, ierror)

  iter = iter + 1
  fmre_old = fmre

call MPI_Barrier( MPI_COMM_WORLD, ierror)

end do

  !--Calculate scalar flux
  phi_new = 0.0d0

  do j=1,n_size

    do i=1,nx
    
      do n=1,nOrds
      
      phi_new(i,j) =  phi_new(i,j) + w(n) * psi(i,j,n)
 
      end do

    end do
    
  end do

call output(nx,n_size,nOrds,psi,psi_mfd,phi_new,phi_mfd,diff,diff2)

    !Find L_inf
    call error_calc3(nx,n_size,dx,dy,nOrds,diff,10,infin)
    call MPI_Reduce( infin, L_infin, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, master, MPI_COMM_WORLD, ierror)

    !Find Linf for scalar flux
    call error_calc2(nx,n_size,dx,dy,diff2,10,smax)
    call MPI_Reduce( smax, sL_infin, 1, MPI_DOUBLE_PRECISION,&
                     MPI_SUM, master, MPI_COMM_WORLD, ierror)

    !Print results
    call MPI_Gather( diff, nx*n_size*nOrds, MPI_DOUBLE_PRECISION, ferror, nx*n_size*nOrds,&
                  MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)
  
    call MPI_Gather( phi_new, nx*n_size, MPI_DOUBLE_PRECISION, phi_tot, nx*n_size,&
                  MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierror)

  if (rank .eq. master) then

    write(*,*) 'ainf: ', infin
    write(*,*) 'sinf: ', smax

    !Print out results
    open(unit=10, file='errora.txt', access='append', status='unknown')
    write(10,*) infin
    close(10)
  
    open(unit=11, file='errors.txt', access='append', status='unknown')
    write(11,*) smax
    close(11)

    open(unit=12,file='ferror.txt', status='unknown')
    do j=1,nx
      write(12,'(ES16.8)') (Abs(Maxval(ferror(i,j,:))), i=1,nx)
    end do
    close(12)

    open(unit=13,file='sflux.txt', status='unknown')
    do j=1,nx
      write(13,'(ES16.8)') (phi_tot(i,j), i=1,nx)
    end do
    close(13)

  end if

deallocate(psi, ss)
deallocate(ext_ss)
deallocate(phi_new)
deallocate(legP_coefs, psi_mfd, phi_mfd)

deallocate(p, phir, phii)
deallocate(phire, phiim)
deallocate(sphre,ferror,phi_tot)
deallocate(fmre, fmim, fmre_old)
deallocate(fm,fmi,snd_bot,snd_top)

deallocate(diff,diff2)
deallocate(xi,eta,mu,w)
deallocate(left,right,top,bottom)
deallocate(l_r,r_r,b_r,t_r,f_r,dfr)

finish = MPI_Wtime()
  if (rank .eq. master) then
  open(unit=14,file='time.txt', access='append', status='unknown')
  write(14,'(A5,I5,A8,ES16.8)') 'P:', num_cores, 'Time:', finish-start
  close(14)
  end if

call MPI_Barrier( MPI_COMM_WORLD, ierror)

call MPI_Finalize( ierror)

end program twoDSn
