       subroutine spharm(nnmax,nn,mm,cc,theta,sphi,cphi,p,
     1                   phire,phiim)
c----------------------------------------------------------------
c	Computes
c
c	Sqrt[(2 n +1)/(4 pi)]*Sqrt[(n-m)!/(n+m)!] P_n^m(cc)
c
c       and its derivative with respect to cc
c
c       as well as real and imaginary parts of E^(i m phi) 
c       without multiplying them (to form spherical harmonics)
c
c	========Input=======
c
c	cc	---	cos(theta)
c	sphi	---	sin(phi)
c	cphi	---	cos(phi)
c	nnmax	---	actual dimension
c	nn	---	n=0,...,nn
c	mm	---	m=0,...,mm
c
c	========output=======
c
c	values p[m,n] and p'[m,n]
c
c       factors phire and phiim
c
c       constant below is 1/Sqrt[4 pi]     
c-------------------------------------------------------------
        implicit real*8 (a-h,o-z)
        real *8 p(0:nnmax,0:nnmax)
        real *8 twonplus(0:nnmax),re(2),re0(2),srint(0:2*nnmax+1)
        real *8 phire(0:mm),phiim(0:mm)
	
	real *8 toll
	
        complex *16 zz0,zz
        equivalence(zz,re(1))
        equivalence(zz0,re0(1))
c
c       -------------------------
c
        zero = 0.0d0
        one  = 1.0d0
	two  = 2.0d0
        oversqrtpi  = 0.2820947917738781434740397257803863d0
	pi  = 4.0d0*atan(1.0d0) 
c
c-----  compute factor from E^(m phi)
c
        re0(1)   = cphi
        re0(2)   = sphi
        phire(0) = one
        phiim(0) = zero
        re(1)    = one
        re(2)    = zero

	toll = dlamch( 's' )
	
c
c------ compute the azimuthal part
c
        do m=1,mm
          zz = zz*zz0
          phire(m) = re(1)
          phiim(m) = re(2)
        enddo     
c
c------ compute common factors for what follows
c
        do m=0,2*max(nn,mm) + 1       
          srint(m) = sqrt(1.0d0*m)
        enddo
c
        do m=0,max(nn,mm)
          twonplus(m) = two*m + one
        enddo
c
c------ initialize
c
        ss = sqrt(one - cc*cc)
        do n=0,nn
          do m=0,mm
             p(m,n) = zero 
          enddo
        enddo
	
        p(0,0) = one
c
c------ check for underflow "near" the poles
c
c        if((theta .lt. 1.0d0) .or. ((pi - theta) .lt. 1.0d0)) then
c	
c          if(pi - theta .lt. 1.0d0) then
c	  
c	    theta = pi - theta
c	    
c	  endif
c
c          if(int(log10(one/toll)/log10(one/theta)).lt.mm) then
c	  
c	    mmax = int(log10(one/toll)/log10(one/theta)) + 1
c	    
c	  else
c	  
c	    mmax = mm
c	    
c	  endif
c
c	else
c
c	  mmax = mm
c
c	endif
        mmax = mm
c
c------ start recursion along diagonal
c
        do m = 0,mmax-1
	
          p(m+1,m+1) = -ss*p(m,m)*srint(2*m+1)/srint(2*m+2)
	  
        enddo
c
c------ work off diagonal, to the right
c
        do m = 0,mmax-1
c 
          p(m,m+1) = cc*p(m,m)*srint(2*m+1)
c
          do n = m+1,nn-1
	  
            p(m,n+1) = cc*p(m,n)*twonplus(n)/srint(n+1+m)/srint(n+1-m)-
     1      p(m,n-1)*srint(n+m)*srint(n-m)/srint(n+1+m)/srint(n+1-m)
     
          enddo
c
        enddo     
c
c------ normalize
c
        do n = 0,nn
c
          do m = 0,n
	  
            p(m,n)  =  p(m,n)*srint(2*n+1)*oversqrtpi
	    
          enddo
c
        enddo 
c	
        return
	
      end
