subroutine  ext_source(i,j,nx,ny,dx,dy,xi,eta,mu,ext_ss)
  !------------------------------------------!
  ! Calculates the external source term      !
  !                                          !
  ! Written by: CDA                          !
  ! Last modified: 4/30/12                   !
  ! Comments: initial writing                !        
  !------------------------------------------!
  implicit none

  double precision, intent(out) :: ext_ss
  double precision              :: pi
  
  !Mesh parameters
  double precision, intent(in)  :: dx, dy
  integer,          intent(in)  :: nx, ny
  double precision, intent(in)  :: xi, eta, mu

  !Integers
  integer :: i, j

  pi = dacos(-1.0d0)


! External source
  if ((i .eq. nx/2) .and. (j .eq. ny/2)) then
    ext_ss = 1.0d0/(dx*dy*4.0d0*pi)
  else
    ext_ss = 0.0d0
  end if

end subroutine ext_source
