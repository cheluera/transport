module func_ang

contains

  subroutine spharm(nmax,ctheta,sphi,cphi,p,phire,phiim)
  !----------------------------------------------------------------
  !	Computes
  !  
  !	Sqrt[(2 n +1)/(4 pi)]*Sqrt[(n-m)!/(n+m)!] P_n^m(ctheta)
  !
  !
  !       as well as real and imaginary parts of E^(i m phi) 
  !       without multiplying them (to form spherical harmonics)
  !
  !	========Input=======
  !
  !	ctheta	---	cos(theta)
  !	sphi	---	sin(phi)
  !	cphi	---	cos(phi)
  !	nmax	---	maximum degree
  !
  !	========output=======
  !
  !	values p[m,n], and Re(Exp(m*phi)), Im(Exp(m*phi))
  !
  !
  !       constant below is 1/Sqrt[4 pi]     
  !-------------------------------------------------------------
          implicit real*8 (a-h,o-z)
          real *8 p(0:nmax,0:nmax)
          real *8 twonplus(0:nmax),re(2),re0(2),srint(0:2*nmax+1)
          real *8 phire(0:nmax),phiim(0:nmax)
          complex *16 zz0,zz
          equivalence(zz,re(1))
          equivalence(zz0,re0(1))
  !
  !       -------------------------
  !
          zero = 0.0d0
          one  = 1.0d0
          two  = 2.0d0
          oversqrtpi  = 0.2820947917738781434740397257803863d0
          pi  = 4.0d0*atan(1.0d0) 
  !
  !-----  compute factor from E^(m phi)  
  !
          re0(1)   = cphi
          re0(2)   = sphi
          phire(0) = one
          phiim(0) = zero
          re(1)    = one
          re(2)    = zero
	
  !
  !------ compute the azimuthal part
  !
          do m=1,nmax

            zz = zz*zz0
            phire(m) = re(1)
            phiim(m) = re(2)

          enddo     
  !
  !------ compute common factors for what follows
  !
          do m=0,2*nmax + 1       

            srint(m) = sqrt(1.0d0*m)

          enddo
  !
          do m=0,nmax

            twonplus(m) = two*m + one

          enddo
  !
  !------ initialize
  !
          ss = sqrt(one - ctheta*ctheta)
          do n=0,nmax

            do m=0,nmax

               p(m,n) = zero 

            enddo

          enddo
	
          p(0,0) = one
  !
  !------ start recursion along diagonal
  !
          do m = 0,nmax - 1
	
            p(m+1,m+1) = -ss*p(m,m)*srint(2*m+1)/srint(2*m+2)
	  
          enddo
  !
  !------ work off diagonal, to the right
  !
          do m = 0,nmax - 1
  ! 
            p(m,m+1) = ctheta*p(m,m)*srint(2*m+1)
  !
            do n = m + 1,nmax - 1
	  
              p(m,n+1) = ctheta*p(m,n)*twonplus(n)/srint(n+1+m)/srint(n+1-m) - &
        &                p(m,n-1)*srint(n+m)*srint(n-m)/srint(n+1+m)/srint(n+1-m)
     
            enddo
  !
          enddo
  !	
  !
  !------ normalize
  !
          do n = 0,nmax
  !
            do m = 0,n
	  
              p(m,n)  =  p(m,n)*srint(2*n+1)*oversqrtpi
	    
            enddo
  !
          enddo 
  !	
          return
	
        end subroutine spharm
!------------------------------------------------------------------------------------
  subroutine f_moments(nx,n_size,psi,w,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim)
    !--------------------------------------------------!
    ! Calculates the angular flux moments. The real    !
    ! and imaginary components are calculated          ! 
    ! seperately using spherical harmonics expanded    ! 
    ! to cosine and sine terms.                        ! 
    !--------------------------------------------------!
    
    implicit none

    double precision, intent(in), dimension(nOrds)  :: w
    double precision, intent(in), dimension(nx,n_size,nOrds)  :: psi
    double precision, intent(in), dimension(nOrds,0:n_sph) :: phir, phii
    double precision, intent(in), dimension(nOrds,0:n_sph,0:n_sph) :: sphre

    double precision, dimension(nx,n_size,0:n_max,0:n_max) :: fmre, fmim
  
    !Mesh parameters
    integer, intent(in) :: nx, n_size
  
    !Integers
    integer :: i,j,k,n,m,n_max,n_min,n_sph,nOrds

    fmre = 0.0d0
    fmim = 0.0d0

  !Calculate angular flux moments

  do i=1,nx

    do j=1,n_size
    
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
!-----------------------------------------------------------------------------------------
  subroutine source(nx,n_size,legP_coefs,n_max,n_min,n_sph,nOrds,sphre,phir,phii,fmre,fmim,ss)
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
    double precision, intent(in), dimension(nx,n_size,0:n_max,0:n_max) :: fmre, fmim

    double precision :: calc, calc2
    double precision, dimension(nx,n_size,nOrds) :: ss
  
    !Mesh parameters
    integer, intent(in) :: nx, n_size

    !Integers
    integer :: i, j, k, n, m, n_max, n_min, n_sph, nOrds
  
  ss = 0.0d0

  ! Scattering source
  do i=1,nx

    do j=1,n_size
    
      do k=1,nOrds
      
        do n=0,n_min
  
        calc2 = fmre(i,j,0,n) * sphre(k,0,n)
      
        calc = 0.0d0

          do m=1,n

          calc = calc + sphre(k,m,n) * ( fmre(i,j,m,n) * phir(k,m) + fmim(i,j,m,n) * phii(k,m) )
 
          end do

        ss(i,j,k) = ss(i,j,k) + legP_coefs(n) * ( calc2 + 2.0d0 * calc )

        end do

      end do

    end do
   
  end do  

  end subroutine source
!------------------------------------------------------------------------------------
end module func_ang
