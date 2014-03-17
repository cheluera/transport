subroutine source(nx,ny,legP_coefs,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim,ss)
  !--------------------------------------------------!
  ! Calculates scattering and external sources. The  !
  ! scattering source is calculated based on the 
  ! spherical harmonic expansion of the flux and 
  ! scattering x-section. 
  !--------------------------------------------------!
  implicit none

  double precision, intent(in), dimension(0:n_max)  :: legP_coefs
  double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
  double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre
  double precision, intent(in), dimension(nx,ny,0:n_max,0:n_max) :: fmre, fmim

  double precision :: calc, calc2
  double precision, dimension(nx,ny,nOrds) :: ss
  
  !Mesh parameters
  !double precision :: dx, dy
  integer :: nx, ny

  !Integers
  integer :: i, j, k, n, m, n_max, n_min, n_sph, nOrds
  
ss = 0.0d0

! Scattering source
do i=1,nx
  
  do j=1,ny
    
    do k=1,nOrds
      
      do n=0,n_min
  
      calc2 = fmre(i,j,0,n) * sphre(k,0,n)
      
      calc = 0

        do m=1,n

        calc = calc + sphre(k,m,n) * ( fmre(i,j,m,n) * phir(k,m) + fmim(i,j,m,n) * phii(k,m) )
 
        end do

      ss(i,j,k) = ss(i,j,k) + legP_coefs(n) * ( calc2 + 2.0d0 * calc )

      end do

    end do
     
  end do
   
end do  

end subroutine source
