subroutine source(i,j,nx,ny,dx,dy,c,phi,legP_coefs,n_moments,ss)
  !--------------------------------------------------!
  ! Calculates scattering and external sources. The  !
  ! scattering source is calculated based on the 
  ! spherical harmonic expansion of the flux and 
  ! scattering x-section. 
  !--------------------------------------------------!

  implicit none

  double precision :: phi, legP_coefs
  double precision :: ss, pi, c
  
  !Mesh parameters
  double precision :: dx, dy
  integer :: nx, ny

  !Integers
  integer :: i, j,n_moments

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

! External source
  if ((i .eq. nx/2) .and. (j .eq. ny/2)) then
    ss = ss + 1.0/(dx*dy*4.0*pi)
  end if

  return

end subroutine source
