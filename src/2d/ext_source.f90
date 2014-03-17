subroutine  ext_source(nx,ny,legP_coefs,c,n_min,n_max,n_moments,n_sph,nOrds,strm,psi_mfd,fm,fmi,f_r,sphre,phir,phii,ext_ss)

  !------------------------------------------!
  ! Calculates the external source term      !
  !------------------------------------------!

  implicit none

  double precision, intent(in), dimension(0:n_max)  :: legP_coefs
  double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
  double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
  double precision, intent(in), dimension(nx,ny,nOrds) :: strm, psi_mfd
  double precision, intent(in), dimension(nx,ny) :: f_r
  
  double precision :: calc, c, calc2
  double precision, dimension(nx,ny,nOrds) :: ext_ss
  double precision, intent(in), dimension(0:n_moments,0:n_moments) :: fm, fmi
  double precision, allocatable, dimension(:) :: sk
  
  !Mesh parameters
  integer :: nx, ny

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

do i =1,nx

  do j=1,ny

    do k=1,nOrds

    ext_ss(i,j,k) = strm(i,j,k) + psi_mfd(i,j,k) - c * f_r(i,j) * sk(k)

    end do
    
  end do

end do

!--Deallocate Memory
deallocate(sk)

end subroutine ext_source
