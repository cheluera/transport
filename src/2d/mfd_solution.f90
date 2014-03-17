subroutine  solution(nx,ny,n_moments,n_sph,nOrds,f_r,l_r,r_r,b_r,t_r,dfr,fm,fmi,sphre,phir,phii,strm,psi_mfd,left,right,top,bottom)

  !------------------------------------------!
  ! Calculates the manufactured solution     !
  !------------------------------------------!

  implicit none

  double precision, intent(in), dimension(0:nx+1,0:ny+1) :: l_r, r_r, b_r, t_r
  double precision, intent(in), dimension(nx,ny) :: f_r
  double precision, intent(in), dimension(nx,ny,nOrds) :: dfr
  double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
  double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
  double precision, intent(in), dimension(0:n_moments,0:n_moments) :: fm, fmi
  
  double precision, dimension(0:nx+1,0:ny+1,nOrds) :: left, right, top, bottom
  double precision, dimension(nx,ny,nOrds) :: psi_mfd, strm
  double precision, allocatable, dimension(:) :: ss
  double precision :: calc, calc2
  
  !Mesh parameters
  integer :: nx, ny

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

do i=1,nx

  do j=1,ny

    do k=1,nOrds

    psi_mfd(i,j,k) = f_r(i,j) * ss(k)

    strm(i,j,k) = dfr(i,j,k) * ss(k)  !Streaming operator

    left(0,j,k) = l_r(0,j) * ss(k)
    
    right(nx+1,j,k) = r_r(nx+1,j) * ss(k)
    
    top(i,ny+1,k) = t_r(i,ny+1) * ss(k)
    
    bottom(i,0,k) = b_r(i,0) * ss(k)
    
    end do
    
  end do

end do

deallocate(ss)

end subroutine solution
