module func_out

contains

  subroutine output(nx,n_size,nOrds,psi,psi_mfd,phi_new,phi_mfd,diff,diff2)

    !-----------------------------------------------------!
    !  Creates error plot data, and scalar flux data.     !
    !-----------------------------------------------------!
  
    implicit none
  
    double precision, intent(in), dimension(nx,n_size,nOrds) :: psi, psi_mfd
    double precision, intent(in), dimension(nx,n_size) :: phi_new, phi_mfd
    integer, intent(in) :: nx, n_size, nOrds
  
    double precision, dimension(nx,n_size) :: diff2
    double precision, dimension(nx,n_size,nOrds) :: diff
    integer :: i,j,n

    !Calculate psi_mfd - psi_num

    do j=1,n_size

      do i=1,nx
    
        do n=1,nOrds
      
        diff(i,j,n) = psi_mfd(i,j,n) - psi(i,j,n)
      
        end do

      end do
    
    end do

    !Calculate phi_mfd - phi_num

    do j=1,n_size

      do i=1,nx
    
        diff2(i,j) = phi_mfd(i,j) - phi_new(i,j)

      end do
    
    end do
  
  end subroutine output
!----------------------------------------------------------------------------
  subroutine error_calc2(nx,n_size,dx,dy,f,p,norm)
    !--------------------------------------------------!
    ! Calculates the L1-norm, L2-norm, and L_inf       !
    ! relative errors.                                 ! 
    !--------------------------------------------------!

    implicit none

    double precision, intent(in) :: dx, dy
    double precision, intent(in), dimension(nx,n_size) :: f
    integer, intent(in) :: nx, n_size, p
  
    double precision :: rsum, norm
  
    !Integers
    integer :: i, j

    rsum = 0.0d0

    !Do maximum error for L-infinity
  
    if (p .gt. 2) then
      
      norm = ABS(MAXVAL(f))
  
    else
  
      !Set up do loops to calculate general norms with a given order (p)

      do j=1,n_size

        do i=1,nx
    
        rsum = rsum + ABS(f(i,j))**(1.0d0*p)

        end do
    
      end do

      rsum = rsum*dx*dy
  
      norm = rsum**(1.0d0/(1.0d0*p))
  
    end if
  
  end subroutine error_calc2
!---------------------------------------------------------------------------
  subroutine error_calc3(nx,n_size,dx,dy,nOrds,f,p,norm)
    !--------------------------------------------------!
    ! Calculates the L1-norm, L2-norm, and L_inf       !
    ! relative errors.                                 ! 
    !--------------------------------------------------!

    implicit none

    double precision, intent(in) :: dx, dy
    double precision, intent(in), dimension(nx,n_size,nOrds) :: f
    integer, intent(in) :: nx, n_size, nOrds, p
  
    double precision :: rsum, norm
  
    !Integers
    integer :: i,j,n

    rsum = 0.0d0

    !Do maximum error for L-infinity
    if (p .gt. 2) then

      norm = ABS(MAXVAL(f))
  
    else

      !Set up do loops to calculate general norms with a given order (p)

      do j=1,n_size

        do i=1,nx

          do n=1,nOrds

            rsum = rsum + ABS(f(i,j,n))**(1.0d0*p)

          end do

        end do

      end do

      rsum = rsum*dx*dy
  
      norm = rsum**(1.0d0/(1.0d0*p))
  
    end if
  
  end subroutine error_calc3

end module func_out
