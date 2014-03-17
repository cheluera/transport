module func_mfd

contains

  subroutine mfd(spatial,rank,master,num_cores,nx,ny,n_size,dx,dy,xi,eta,&
                 nOrds,n_moments,l_r,r_r,b_r,t_r,f_r,dfr,fm,fmi,phi_mfd) 

  !---------------------------------------------------------!
  ! Routine to calculate manufactured solution using        !
  ! spherical harmonics.  Also evaluates boundary values and!
  ! external sources.                                       !
  !                                                         !
  !---------------------------------------------------------!

  implicit none

  integer, intent(in) :: nx, ny, n_size, nOrds, n_moments
  integer, intent(in) :: rank, master, num_cores
  double precision, intent(in) :: dx, dy
  double precision, intent(in), dimension(nOrds) :: xi, eta
  logical, intent(in) :: spatial

  double precision, dimension(nx,ny) :: phi_mfd
  double precision, dimension(0:nx,n_size) :: l_r
  double precision, dimension(nx+1,n_size) :: r_r
  double precision, dimension(nx,0:n_size) :: b_r
  double precision, dimension(nx,n_size+1) :: t_r
  double precision, dimension(nx,n_size) :: f_r
  double precision, dimension(nx,n_size,nOrds) :: dfr
  double precision :: pi

  !Angular flux moments
  double precision, dimension(0:n_moments,0:n_moments) :: fm, fmi

  !Integers
  integer :: i, j, n 

  !===========================================================!
  !===========================================================!

  pi = 3.1415926535897932385d0

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
  fm(0,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(2,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(4,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(6,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(8,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(10,10) = ( 3.0d0 ) ** (-10.0d0)
  fm(1,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(3,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(5,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(7,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(9,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(11,11) = ( 3.0d0 ) ** (-11.0d0)
  fm(0,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(2,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(4,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(6,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(8,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(10,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(12,12) = ( 3.0d0 ) ** (-12.0d0)
  fm(1,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(3,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(5,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(7,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(9,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(11,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(13,13) = ( 3.0d0 ) ** (-13.0d0)
  fm(0,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(2,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(4,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(6,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(8,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(10,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(12,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(14,14) = ( 3.0d0 ) ** (-14.0d0)
  fm(1,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(3,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(5,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(7,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(9,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(11,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(13,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(15,15) = ( 3.0d0 ) ** (-15.0d0)
  fm(0,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(2,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(4,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(6,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(8,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(10,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(12,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(14,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(16,16) = ( 3.0d0 ) ** (-16.0d0)
  fm(1,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(3,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(5,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(7,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(9,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(11,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(13,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(15,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(17,17) = ( 3.0d0 ) ** (-17.0d0)
  fm(0,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(2,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(4,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(6,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(8,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(10,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(12,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(14,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(16,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(18,18) = ( 3.0d0 ) ** (-18.0d0)
  fm(1,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(3,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(5,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(7,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(9,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(11,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(13,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(15,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(17,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(19,19) = ( 3.0d0 ) ** (-19.0d0)
  fm(0,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(2,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(4,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(6,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(8,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(10,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(12,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(14,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(16,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(18,20) = ( 3.0d0 ) ** (-20.0d0)
  fm(20,20) = ( 3.0d0 ) ** (-20.0d0)

  fmi = fm

  !Designate f(r) and Omega*del(f(r)) from manufactured solution
  l_r = 0.0d0
  r_r = 0.0d0
  b_r = 0.0d0
  t_r = 0.0d0

  select case (spatial)

  case (.false.)

  f_r = ( 4.0d0 * pi ) ** (-0.50d0)
  l_r(0,:) = ( 4.0d0 * pi ) ** (-0.50d0)
  r_r(nx+1,:) = ( 4.0d0 * pi ) ** (-0.50d0)
  dfr = 0.0d0

  if (rank .eq. master) then  
  
    b_r(:,0) = ( 4.0d0 * pi ) ** (-0.50d0)

  else if (rank .eq. num_cores-1) then

  t_r(:,n_size+1) = ( 4.0d0 * pi ) ** (-0.50d0)

  end if

  case (.true.)

  do j=1,n_size

    do i=1,nx

      f_r(i,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0+(((i-0.5d0)*dx)/10.d0)**6.0d0 + &
                  ((((j+rank*n_size)-0.5d0)*dy)/10.0d0)**6.0d0)

    end do

  end do

  do j=1,n_size

    l_r(0,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + ((((j+rank*n_size)-0.5d0)*dy)/10.0d0)**6.0d0)
    r_r(nx+1,j) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0+((nx*dx)/10.0d0)**6.0d0 + &
                  ((((j+rank*n_size)-0.5d0)*dy)/10.0d0)**6.0d0) 

  end do

  if (rank .eq. master) then

  do i=1,nx
    b_r(i,0) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + (((i-0.5d0)*dx)/10.0d0)**6.0d0)
  end do

  else if (rank .eq. num_cores-1) then

  do i=1,nx

    t_r(i,n_size+1) = ( 4.0d0 * pi) ** (-0.50d0) * (1.0d0 + (((i-0.5d0)*dx)/10.0d0)**6.0d0 +&
                      ((ny*dy)/10.0d0)**6.0d0)
  end do

  end if

  do j=1,n_size

    do i=1,nx

      do n=1,nOrds

        dfr(i,j,n) = ( 4.0d0 * pi ) ** (-0.50d0) * 6.0d0/10.0d0 * ( xi(n) * (((i-0.5d0)*dx)/10.0d0)**5.0d0 +&
                    eta(n) * ((((j+rank*n_size)-0.5d0)*dy)/10.0d0)**5.0d0 )

      end do

    end do

  end do

  end select

    !--Calculate manufactured scalar flux

  do j=1,n_size

    do i=1,nx

        phi_mfd(i,j) = ( 4.0d0 * pi ) ** (0.5d0) * fm(0,0) * f_r(i,j)
    
    end do

  end do 

  end subroutine mfd
!--------------------------------------------------------------------------------
  subroutine  solution(rank,master,num_cores,nx,n_size,n_moments,n_sph,nOrds,f_r,l_r,r_r,b_r,t_r,dfr,&
                       fm,fmi,sphre,phir,phii,strm,psi_mfd,left,right,top,bottom)

    !------------------------------------------!
    ! Calculates the manufactured solution     !
    !------------------------------------------!

    implicit none

    double precision, intent(in), dimension(0:nx,n_size) :: l_r
    double precision, intent(in), dimension(nx+1,n_size) :: r_r
    double precision, intent(in), dimension(nx,0:n_size) :: b_r
    double precision, intent(in), dimension(nx,n_size+1) :: t_r
    double precision, intent(in), dimension(nx,n_size) :: f_r
    double precision, intent(in), dimension(nx,n_size,nOrds) :: dfr
    double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
    double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
    double precision, intent(in), dimension(0:n_moments,0:n_moments) :: fm, fmi
    integer, intent(in) :: rank, master, num_cores, n_size
  
    double precision, dimension(0:nx,n_size,nOrds) :: left
    double precision, dimension(nx+1,n_size,nOrds) :: right
    double precision, dimension(nx,0:n_size,nOrds) :: bottom
    double precision, dimension(nx,n_size+1,nOrds) :: top
    double precision, dimension(nx,n_size,nOrds) :: psi_mfd, strm
    double precision, allocatable, dimension(:) :: ss
    double precision :: calc, calc2
  
    !Mesh parameters
    integer, intent(in) :: nx

    !Integers
    integer :: i, j, k, n, m, n_moments, n_sph, nOrds
  
    !Allocate memory
    allocate(ss(nOrds))

  ss = 0.0d0

  !Calculate solution from spherical harmonics and flux moments

      do k=1,nOrds

        do n=0,n_moments

        calc = fm(0,n) * sphre(k,0,n)

        calc2 = 0.0d0

          do m=1,n
            
          calc2 = calc2 + sphre(k,m,n) * ( fm(m,n) * phir(k,m) + fmi(m,n) * phii(k,m) )

          end do
	
        ss(k) = ss(k) + ( calc + 2.0d0 * calc2 )

        end do

      end do
  
  left = 0.0d0
  right = 0.0d0
  top = 0.0d0
  bottom = 0.0d0

  do j=1,n_size

    do i=1,nx

      do k=1,nOrds

      psi_mfd(i,j,k) = f_r(i,j) * ss(k)

      strm(i,j,k) = dfr(i,j,k) * ss(k)  !Streaming operator

      left(0,j,k) = l_r(0,j) * ss(k)
    
      right(nx+1,j,k) = r_r(nx+1,j) * ss(k)
    
      end do

    end do

  end do

  if (rank .eq. master) then

    do i=1,nx

      do k=1,nOrds

      bottom(i,0,k) = b_r(i,0) * ss(k)

      end do

    end do

  else if (rank .eq. num_cores-1) then

    do i=1,nx

      do k=1,nOrds

      top(i,n_size+1,k) = t_r(i,n_size+1) * ss(k)

      end do

    end do

  end if

  deallocate(ss)

  end subroutine solution
!-----------------------------------------------------------------------
  subroutine  ext_source(nx,n_size,legP_coefs,c,n_min,n_max,n_moments,n_sph,&
                         nOrds,strm,psi_mfd,fm,fmi,f_r,sphre,phir,phii,ext_ss)

    !------------------------------------------!
    ! Calculates the external source term      !
    !------------------------------------------!

    implicit none

    double precision, intent(in), dimension(0:n_max)  :: legP_coefs
    double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
    double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
    double precision, intent(in), dimension(nx,n_size,nOrds) :: strm, psi_mfd
    double precision, intent(in), dimension(nx,n_size) :: f_r
  
    double precision :: calc, c, calc2
    double precision, dimension(nx,n_size,nOrds) :: ext_ss
    double precision, intent(in), dimension(0:n_moments,0:n_moments) :: fm, fmi
    double precision, allocatable, dimension(:) :: sk
  
    !Mesh parameters
    integer, intent(in) :: nx, n_size

    !Integers
    integer :: i, j, k, n, m, n_max, n_min, n_moments, n_sph, nOrds

  !--allocate scattering source
  allocate(sk(nOrds))

  sk = 0.0d0

  ! Scattering source component of L operator

      do k=1,nOrds
  
        do n=0,n_min
    
        calc = fm(0,n) * sphre(k,0,n)

        calc2 = 0.0d0
	   
          do m=1,n

          calc2 = calc2 + sphre(k,m,n) * ( fm(m,n) * phir(k,m) + fmi(m,n) * phii(k,m) )

          end do

        sk(k) = sk(k) + legP_coefs(n) * ( calc + 2.0d0 * calc2 )

        end do
    
      end do

  do j=1,n_size

    do i=1,nx

      do k=1,nOrds

      ext_ss(i,j,k) = strm(i,j,k) + psi_mfd(i,j,k) - c * f_r(i,j) * sk(k)

      end do

    end do

  end do

  !--Deallocate Memory
  deallocate(sk)

  end subroutine ext_source

end module func_mfd
