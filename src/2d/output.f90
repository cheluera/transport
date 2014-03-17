subroutine output(nx,ny,dx,dy,nOrds,psi,psi_mfd,phi_new,phi_mfd)

  !-----------------------------------------------------!
  !  Creates error plot data, and scalar flux data.     !
  !-----------------------------------------------------!
  
  implicit none
  
  double precision, intent(in) :: dx, dy
  double precision, intent(in), dimension(nx,ny,nOrds) :: psi, psi_mfd
  double precision, intent(in), dimension(nx,ny) :: phi_new, phi_mfd
  integer, intent(in) :: nx, ny, nOrds
  
  double precision, allocatable, dimension(:,:) :: diff2
  double precision, allocatable, dimension(:,:,:) :: diff
  double precision :: mnorm, infin
  double precision :: snorm3, smax
  integer :: i,j,n

  !Allocate memory
  allocate(diff(nx,ny,nOrds))
  allocate(diff2(nx,ny))

  !Calculate psi_mfd - psi_num

  do i=1,nx
  
    do j=1,ny
    
      do n=1,nOrds
      
      diff(i,j,n) = psi_mfd(i,j,n) - psi(i,j,n)
      
      end do
      
    end do
    
  end do

  !Calculate phi_mfd - phi_num

  do i=1,nx
  
    do j=1,ny
    
    diff2(i,j) = phi_mfd(i,j) - phi_new(i,j)
    
    end do
    
  end do


  !Find L_inf
  
  !Max(|psi_mfd(i,j,n) - psi(i,j,n)|)
  call error_calc3(nx,ny,dx,dy,nOrds,diff,10,mnorm)
  
  infin = mnorm

  
  !Find Linf for scalar flux
  
  call error_calc2(nx,ny,dx,dy,diff2,10,snorm3)
  
  smax = snorm3

  write(*,*) 'ainf: ', infin
  write(*,*) 'sinf: ', smax

  !Print out results
  open(unit=7,file='ferror.txt',status='unknown')
  do i=1,nx
    write(7,'(ES16.8)') (Abs(Maxval(diff(i,j,:))),j=1,ny)
  end do
  close(7)
  
  open(unit=7,file='sflux.txt',status='unknown')
  do i=1,nx
    write(7,'(ES16.8)') (Abs(phi_new(i,j)),j=1,ny)
  end do
  close(7)
  
  open(unit=7, file='errora.txt', access='append', status='unknown')
  write(7,*) infin
  close(7)
  
  open(unit=7, file='errors.txt', access='append', status='unknown')
  write(7,*) smax
  close(7)

  !Deallocate memory
  deallocate(diff,diff2)
  
end subroutine output
