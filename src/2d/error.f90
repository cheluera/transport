subroutine error_calc(nx,ny,dx,dy,w,nOrds,f,g,p,dim,norm)
  !--------------------------------------------------!
  ! Calculates the L1-norm, L2-norm, and L_inf       !
  ! relative errors.                                 ! 
  !--------------------------------------------------!

  implicit none

  double precision, intent(in) :: dx, dy
  double precision, intent(in), dimension(nOrds) :: w
  double precision, intent(in), dimension(nx,ny) :: g
  double precision, intent(in), dimension(nx,ny,nOrds) :: f
  integer, intent(in) :: nx, ny, nOrds, p, dim
  
  double precision :: rsum, norm
  
  !Integers
  integer :: i,j,n

  rsum = 0.0d0

  !Do maximum error for L-infinity
  if (p .gt. 2) then

    if (dim .eq. 3) then
      
    norm = MAXVAL(ABS(f))
    
    else
    
    norm = Maxval(abs(g))
    
    end if
  
  else
    !Set up do loops to calculate general norms with a given order (p)

    if (dim .eq. 3) then

    do i=1,nx

      do j=1,ny

        do n=1,nOrds

          rsum = rsum + w(n)*ABS(f(i,j,n))**(1.0d0*p)

        end do

      end do

    end do
  
    else

    do i=1,nx
  
      do j=1,ny
    
      rsum = rsum + ABS(g(i,j))**(1.0d0*p)
    
      end do
    
    end do
  
    end if

    rsum = rsum*dx*dy
  
    norm = rsum**(1.0d0/(1.0d0*p))
  
  end if
  
end subroutine error_calc
