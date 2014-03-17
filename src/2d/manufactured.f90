subroutine mfd(nx,ny,dx,dy,xi,eta,nOrds,n_moments,n_max,n_min,n_sph,legP_coefs,c,sphre,phir,phii,left,right,bottom,&
               top,ext_ss,psi_mfd,phi_mfd) 

!---------------------------------------------------------!
! Routine to calculate manufactured solution using        !
! spherical harmonics.  Also evaluates boundary values and!
! external sources.                                       !
!                                                         !
!---------------------------------------------------------!
implicit none

integer, intent(in) :: nx, ny, nOrds, n_moments, n_max, n_sph, n_min
double precision, intent(in) :: c, dx, dy
double precision, intent(in), dimension(nOrds) :: xi, eta
double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
double precision, intent(in), dimension(0:n_max) :: legP_coefs

double precision, dimension(nx,ny,nOrds) :: psi_mfd, strm
double precision, dimension(nx,ny,nOrds) :: ext_ss
double precision, dimension(nx,ny) :: phi_mfd
double precision, dimension(0:nx+1,0:ny+1,nOrds) :: left, right, top, bottom
double precision, allocatable, dimension(:,:,:) :: dfr
double precision, allocatable, dimension(:,:)   :: f_r, l_r, r_r, b_r, t_r

double precision :: pi

!Angular flux moments
double precision, allocatable, dimension(:,:) :: fm, fmi

!Integers
integer :: i, j, n 

!===========================================================!
!===========================================================!

pi = 3.1415926535897932385d0

!Allocate memory
allocate( f_r(nx,ny))
allocate( l_r(0:nx+1,0:ny+1))
allocate( r_r(0:nx+1,0:ny+1))
allocate( b_r(0:nx+1,0:ny+1))
allocate( t_r(0:nx+1,0:ny+1))
allocate( dfr(nx,ny,nOrds))

!Angular flux moments
allocate(fm(0:n_moments,0:n_moments))
allocate(fmi(0:n_moments,0:n_moments))

!Angular flux moments (prescribed)
fm = 0.0d0
fm(0,0) = 1.0d0
fm(1,1) = 1.0d0/3.0d0
fm(0,2) = ( 3.0d0 ) ** (-2.0d0)
fm(2,2) = ( 3.0d0 ) ** (-2.0d0)
fm(1,3) = ( 3.0d0 ) ** (-3.0d0)
fm(3,3) = ( 3.0d0 ) ** (-3.0d0)
fm(0,4) = ( 3.0d0 ) ** (-4.0d0)
fm(2,4) = ( 3.0d0 ) ** (-4.0d0)
fm(4,4) = ( 3.0d0 ) ** (-4.0d0)
fm(1,5) = ( 3.0d0 ) ** (-5.0d0)
fm(3,5) = ( 3.0d0 ) ** (-5.0d0)
fm(5,5) = ( 3.0d0 ) ** (-5.0d0)
fm(0,6) = ( 3.0d0 ) ** (-6.0d0)
fm(2,6) = ( 3.0d0 ) ** (-6.0d0)
fm(4,6) = ( 3.0d0 ) ** (-6.0d0)
fm(6,6) = ( 3.0d0 ) ** (-6.0d0)
fm(1,7) = ( 3.0d0 ) ** (-7.0d0)
fm(3,7) = ( 3.0d0 ) ** (-7.0d0)
fm(5,7) = ( 3.0d0 ) ** (-7.0d0)
fm(7,7) = ( 3.0d0 ) ** (-7.0d0)
fm(0,8) = ( 3.0d0 ) ** (-8.0d0)
fm(2,8) = ( 3.0d0 ) ** (-8.0d0)
fm(4,8) = ( 3.0d0 ) ** (-8.0d0)
fm(6,8) = ( 3.0d0 ) ** (-8.0d0)
fm(8,8) = ( 3.0d0 ) ** (-8.0d0)
fm(1,9) = ( 3.0d0 ) ** (-9.0d0)
fm(3,9) = ( 3.0d0 ) ** (-9.0d0)
fm(5,9) = ( 3.0d0 ) ** (-9.0d0)
fm(7,9) = ( 3.0d0 ) ** (-9.0d0)
fm(9,9) = ( 3.0d0 ) ** (-9.0d0)
fm(0,10)= ( 3.0d0 ) ** (-10.0d0)
fm(2,10)= ( 3.0d0 ) ** (-10.0d0)
fm(4,10)= ( 3.0d0 ) ** (-10.0d0)
fm(6,10)= ( 3.0d0 ) ** (-10.0d0)
fm(8,10)= ( 3.0d0 ) ** (-10.0d0)
fm(10,10)= ( 3.0d0 ) ** (-10.0d0)
fm(1,11)= ( 3.0d0 ) ** (-11.0d0)
fm(3,11)= ( 3.0d0 ) ** (-11.0d0)
fm(5,11)= ( 3.0d0 ) ** (-11.0d0)
fm(7,11)= ( 3.0d0 ) ** (-11.0d0)
fm(9,11)= ( 3.0d0 ) ** (-11.0d0)
fm(11,11)= ( 3.0d0 ) ** (-11.0d0)
fm(0,12)= ( 3.0d0 ) ** (-12.0d0)
fm(2,12)= ( 3.0d0 ) ** (-12.0d0)
fm(4,12)= ( 3.0d0 ) ** (-12.0d0)
fm(6,12)= ( 3.0d0 ) ** (-12.0d0)
fm(8,12)= ( 3.0d0 ) ** (-12.0d0)
fm(10,12)= ( 3.0d0 ) ** (-12.0d0)
fm(12,12)= ( 3.0d0 ) ** (-12.0d0)
fm(1,13)= ( 3.0d0 ) ** (-13.0d0)
fm(3,13)= ( 3.0d0 ) ** (-13.0d0)
fm(5,13)= ( 3.0d0 ) ** (-13.0d0)
fm(7,13)= ( 3.0d0 ) ** (-13.0d0)
fm(9,13)= ( 3.0d0 ) ** (-13.0d0)
fm(11,13)= ( 3.0d0 ) ** (-13.0d0)
fm(13,13)= ( 3.0d0 ) ** (-13.0d0)
fm(0,14)= ( 3.0d0 ) ** (-14.0d0)
fm(2,14)= ( 3.0d0 ) ** (-14.0d0)
fm(4,14)= ( 3.0d0 ) ** (-14.0d0)
fm(6,14)= ( 3.0d0 ) ** (-14.0d0)
fm(8,14)= ( 3.0d0 ) ** (-14.0d0)
fm(10,14)= ( 3.0d0 ) ** (-14.0d0)
fm(12,14)= ( 3.0d0 ) ** (-14.0d0)
fm(14,14)= ( 3.0d0 ) ** (-14.0d0)
fm(1,15)= ( 3.0d0 ) ** (-15.0d0)
fm(3,15)= ( 3.0d0 ) ** (-15.0d0)
fm(5,15)= ( 3.0d0 ) ** (-15.0d0)
fm(7,15)= ( 3.0d0 ) ** (-15.0d0)
fm(9,15)= ( 3.0d0 ) ** (-15.0d0)
fm(11,15)= ( 3.0d0 ) ** (-15.0d0)
fm(13,15)= ( 3.0d0 ) ** (-15.0d0)
fm(15,15)= ( 3.0d0 ) ** (-15.0d0)
fm(0,16)= ( 3.0d0 ) ** (-16.0d0)
fm(2,16)= ( 3.0d0 ) ** (-16.0d0)
fm(4,16)= ( 3.0d0 ) ** (-16.0d0)
fm(6,16)= ( 3.0d0 ) ** (-16.0d0)
fm(8,16)= ( 3.0d0 ) ** (-16.0d0)
fm(10,16)= ( 3.0d0 ) ** (-16.0d0)
fm(12,16)= ( 3.0d0 ) ** (-16.0d0)
fm(14,16)= ( 3.0d0 ) ** (-16.0d0)
fm(16,16)= ( 3.0d0 ) ** (-16.0d0)
fm(1,17)= ( 3.0d0 ) ** (-17.0d0)
fm(3,17)= ( 3.0d0 ) ** (-17.0d0)
fm(5,17)= ( 3.0d0 ) ** (-17.0d0)
fm(7,17)= ( 3.0d0 ) ** (-17.0d0)
fm(9,17)= ( 3.0d0 ) ** (-17.0d0)
fm(11,17)= ( 3.0d0 ) ** (-17.0d0)
fm(13,17)= ( 3.0d0 ) ** (-17.0d0)
fm(15,17)= ( 3.0d0 ) ** (-17.0d0)
fm(17,17)= ( 3.0d0 ) ** (-17.0d0)
fm(0,18)= ( 3.0d0 ) ** (-18.0d0)
fm(2,18)= ( 3.0d0 ) ** (-18.0d0)
fm(4,18)= ( 3.0d0 ) ** (-18.0d0)
fm(6,18)= ( 3.0d0 ) ** (-18.0d0)
fm(8,18)= ( 3.0d0 ) ** (-18.0d0)
fm(10,18)= ( 3.0d0 ) ** (-18.0d0)
fm(12,18)= ( 3.0d0 ) ** (-18.0d0)
fm(14,18)= ( 3.0d0 ) ** (-18.0d0)
fm(16,18)= ( 3.0d0 ) ** (-18.0d0)
fm(18,18)= ( 3.0d0 ) ** (-18.0d0)
fm(1,19)= ( 3.0d0 ) ** (-19.0d0)
fm(3,19)= ( 3.0d0 ) ** (-19.0d0)
fm(5,19)= ( 3.0d0 ) ** (-19.0d0)
fm(7,19)= ( 3.0d0 ) ** (-19.0d0)
fm(9,19)= ( 3.0d0 ) ** (-19.0d0)
fm(11,19)= ( 3.0d0 ) ** (-19.0d0)
fm(13,19)= ( 3.0d0 ) ** (-19.0d0)
fm(15,19)= ( 3.0d0 ) ** (-19.0d0)
fm(17,19)= ( 3.0d0 ) ** (-19.0d0)
fm(19,19)= ( 3.0d0 ) ** (-19.0d0)
fm(0,20)= ( 3.0d0 ) ** (-20.0d0)
fm(2,20)= ( 3.0d0 ) ** (-20.0d0)
fm(4,20)= ( 3.0d0 ) ** (-20.0d0)
fm(6,20)= ( 3.0d0 ) ** (-20.0d0)
fm(8,20)= ( 3.0d0 ) ** (-20.0d0)
fm(10,20)= ( 3.0d0 ) ** (-20.0d0)
fm(12,20)= ( 3.0d0 ) ** (-20.0d0)
fm(14,20)= ( 3.0d0 ) ** (-20.0d0)
fm(16,20)= ( 3.0d0 ) ** (-20.0d0)
fm(18,20)= ( 3.0d0 ) ** (-20.0d0)
fm(20,20)= ( 3.0d0 ) ** (-20.0d0)
fmi = fm

!Designate f(r) and Omega*del(f(r)) from manufactured solution
!f_r = ( 4.0d0 * pi ) ** (-0.50d0)
!l_r = ( 4.0d0 * pi ) ** (-0.50d0)
!r_r = ( 4.0d0 * pi ) ** (-0.50d0)
!b_r = ( 4.0d0 * pi ) ** (-0.50d0)
!t_r = ( 4.0d0 * pi ) ** (-0.50d0)
!dfr = 0.0d0

do i=1,nx
  do j=1,ny
  f_r(i,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0+(((i-0.5d0)*dx)/10.d0)**6.0d0 + &
              (((j-0.5d0)*dy)/10.0d0)**6.0d0)
  end do
end do

do j=1,ny
  l_r(0,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + (((j-0.5d0)*dy)/10.0d0)**6.0d0)
  r_r(nx+1,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0+((nx*dx)/10.0d0)**6.0d0 + &
                (((j-0.5d0)*dy)/10.0d0)**6.0d0) 
end do
do i=1,nx
  b_r(i,0) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + (((i-0.5d0)*dx)/10.0d0)**6.0d0)
  t_r(i,ny+1) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + (((i-0.5d0)*dx)/10.0d0)**6.0d0 +&
                ((ny*dy)/10.0d0)**6.0d0)
end do

do i=1,nx
  do j=1,ny
    do n=1,nOrds
    dfr(i,j,n) = ( 4.0d0 * pi ) ** (-0.50d0) * 6.0d0/10.0d0 * ( xi(n) * (((i-0.5d0)*dx)/10.0d0)**5.0d0 +&
                 eta(n) * (((j-0.5d0)*dy)/10.0d0)**5.0d0 )
    end do
  end do
end do

  !--Calculate manufactured solution
  call solution(nx,ny,n_moments,n_sph,nOrds,f_r,l_r,r_r,b_r,t_r,dfr,fm,fmi,sphre,phir,phii,strm,psi_mfd,left,right,top,bottom)

  !--external source via manufactured solution
  call ext_source(nx,ny,legP_coefs,c,n_min,n_max,n_moments,n_sph,nOrds,strm,psi_mfd,fm,fmi,f_r,sphre,phir,phii,ext_ss)

  !--Calculate manufactured scalar flux

  do i=1,nx
  
    do j=1,ny

      phi_mfd(i,j) = ( 4.0d0 * pi ) ** (0.5d0) * fm(0,0) * f_r(i,j)

    end do
    
  end do 

!--Deallocate memory
deallocate(f_r,dfr,fm,fmi)
deallocate(l_r,r_r,b_r,t_r)
 
end subroutine mfd
