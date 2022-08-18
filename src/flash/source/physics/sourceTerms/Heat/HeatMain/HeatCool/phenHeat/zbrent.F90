!if the values are not bracketed i.e. [>0,<0], then the routine is stopped
	FUNCTION zbrent(func,x1,x2,tol,tmp,dt,sdot,den_H,rho,kconst)

!       USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: x1,x2,tol,tmp,dt,sdot,den_H,rho,kconst
        REAL*8 :: zbrent
        INTERFACE
! get da funk 
                FUNCTION func(x,tmpf,dtf,rhof,sdotf,den_Hf,kconst)
                !USE nrtype
                IMPLICIT NONE
                REAL*8, INTENT(IN) :: x,tmpf,dtf,sdotf,den_Hf,rhof,kconst
                REAL*8 :: func
                END FUNCTION func

        END INTERFACE
        INTEGER, PARAMETER :: ITMAX=1000
        REAL*8, PARAMETER :: EPS=epsilon(x1)
        INTEGER :: iter
        REAL*8 :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm

        a=x1
        b=x2
        fa=func(a,tmp,dt,rho,sdot,den_H,kconst)
        fb=func(b,tmp,dt,rho,sdot,den_H,kconst)
	!print*,'rho',rho
        !print*,'fa,fb,x1,x2'
	!print*,fa,fb,x1,x2
        if ((fa > 0.0d0 .and. fb > 0.0d0) .or. (fa < 0.0d0 .and. fb < 0.0d0)) then
               print*,'root must be bracketed for zbrent'
               print*,'fa,fb,x1,x2,tmp,dt,rho,sdot,den_h'
		print*,fa,fb,x1,x2,tmp,dt,rho,sdot,den_h
                zbrent = -1d0
                return
        endif

        c=b
        fc=fb
        do iter=1,ITMAX
                if ((fb > 0.0d0 .and. fc > 0.0d0) .or. (fb < 0.0d0 .and. fc < 0.0d0)) then
                        c=a
                        fc=fa
                        d=b-a
                        e=d                          
		      end if
		           if (abs(fc) < abs(fb)) then
		                   a=b
		                   b=c
		                   c=a
		                   fa=fb
		                   fb=fc
		                   fc=fa
		           end if
		           tol1=2.0d0*EPS*abs(b)+0.5d0*tol
		           xm=0.5d0*(c-b)
		           if (abs(xm) <= tol1 .or. fb == 0.0) then
		                   zbrent=b
		                   RETURN
		           end if
		           if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
		                   s=fb/fa
		                   if (a == c) then
		                           p=2.0d0*xm*s
		                           q=1.0d0-s
		                   else
		                           q=fa/fc
		                           r=fb/fc
		                           p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
		                           q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
		                   end if
		                   if (p > 0.0d0) q=-q
		                   p=abs(p)
		                   if (2.0d0*p  <  min(3.0d0*xm*q-abs(tol1*q),abs(e*q))) then
		                           e=d
		                           d=p/q
		                   else
		                           d=xm
		                           e=d
		                   end if
		           else
		                   d=xm
		                   e=d
		           end if
		           a=b
		           fa=fb
		           b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
        		   fb=func(b,tmp,dt,rho,sdot,den_H,kconst)
		   end do

		   print*,'zbrent: exceeded maximum iterations'
		   zbrent=b
		   END FUNCTION zbrent
