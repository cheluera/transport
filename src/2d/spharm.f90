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
	
      end
