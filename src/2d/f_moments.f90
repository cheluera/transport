subroutine f_moments(nx,ny,psi,w,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim)
  !--------------------------------------------------!
  ! Calculates the angular flux moments. The real    !
  ! and imaginary components are calculated          ! 
  ! seperately using spherical harmonics expanded    ! 
  ! to cosine and sine terms.                        ! 
  !--------------------------------------------------!
  implicit none

  double precision, intent(in), dimension(nOrds)  :: w
  double precision, intent(in), dimension(nx,ny,nOrds)  :: psi
  double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
  double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre

  double precision, dimension(nx,ny,0:n_max,0:n_max) :: fmre, fmim
  
  !Mesh parameters
  integer :: nx,ny
  
  !Integers
  integer :: i,j,k,n,m,n_max,n_min,n_sph,nOrds

  fmre = 0.0d0
  fmim = 0.0d0

!Calculate angular flux moments

do i=1,nx
  
  do j=1,ny
    
    do n=0,n_min
      
      do m=0,n
	
        ! Real component of angular flux moment
        do k=1,nOrds
	
        fmre(i,j,m,n) = fmre(i,j,m,n) + w(k) * sphre(k,m,n) * phir(k,m) * psi(i,j,k)        

        ! Imaginary component of angular flux moment
	
        fmim(i,j,m,n) = fmim(i,j,m,n) + w(k) * sphre(k,m,n) * phii(k,m) * psi(i,j,k)
	
        end do
        
      end do
	
    end do
      
  end do
    
end do

end subroutine f_moments
