subroutine source(psi,legP_coefs,n_moments,ss)
  !--------------------------------------------------!
  ! Calculates scattering and external sources. The  !
  ! scattering source is calculated based on the 
  ! spherical harmonic expansion of the flux and 
  ! scattering x-section. 
  !--------------------------------------------------!
  implicit none

  double precision, intent(in)  :: psi

  double precision, intent(in)  :: legP_coefs


  double precision :: phi
  double precision :: ss, pi, c
  
  !Mesh parameters
  double precision :: dx, dy
  integer :: nx, ny

  !Integers
  integer :: i, j, n_moments

  pi = dacos(-1.0d0)

! Scattering source
! do n=0,n_moments
!
!    do m=-n,n
!            
!      calculate (n,m)th spherical harmonic moment of psi    --- should do this at the 
!      
!    end do
!
!  end do

  ss = c * phi/(4.0*pi) 


  return

end subroutine source
