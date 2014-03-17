subroutine error_calc2(nx,ny,dx,dy,f,p,norm)
  !--------------------------------------------------!
  ! Calculates the L1-norm, L2-norm, and L_inf       !
  ! relative errors.                                 ! 
  !--------------------------------------------------!

  implicit none

  double precision, intent(in) :: dx, dy
  double precision, intent(in), dimension(nx,ny) :: f
  integer, intent(in) :: nx, ny, p
  
  double precision :: rsum, norm
  
  !Integers
  integer :: i,j

  rsum = 0.0d0

  !Do maximum error for L-infinity
  
  if (p .gt. 2) then
      
    norm = ABS(MAXVAL(f))
  
  else
  
    !Set up do loops to calculate general norms with a given order (p)

    do i=1,nx
  
      do j=1,ny
    
      rsum = rsum + ABS(f(i,j))**(1.0d0*p)
    
      end do
    
    end do

    rsum = rsum*dx*dy
  
    norm = rsum**(1.0d0/(1.0d0*p))
  
  end if
  
end subroutine error_calc2
