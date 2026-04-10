  module pt_sourceUtil
    implicit none

    contains

! piecewise constant cross section + v^-3 drop off at upper end
!	subroutine getSourceEmissionH2(Teff, lowlimit, uplimit, phEnergy, phNumber)
   

!	end subroutine getSourceEmission
! calculates number photons through stellar surface
! and energy per emitted photon
!	subroutine getPhotonFluxes(starRad, starTemp, phEnergy, phFlux, ionlimit, phWEnergy)
! assumes v^-3 drop in cross section, so good for ionizing radiation
! photon number is good though for anything
  subroutine getSourceEmission(Teff, lowlimit, uplimit, sigma, phEnergy, phNumber, phActive)

    implicit none

    real, intent(in)  :: Teff, lowlimit, uplimit,sigma
    real, intent(out) :: phEnergy, phNumber, phActive

    integer*8 :: i  
    integer   :: iter
    real*8    :: ll, ul, tmp
    real :: stboltz,  energy, kT, prefac, PI2
    real :: prefac2, solarRadius, weightedE, prefac3
    real :: nPhotll, nPhotul, energyll, energyul, ratell, rateul

! hardcoded for portability:
    real, parameter :: rt_boltz   = 1.3806488e-16 
    real, parameter :: rt_speedlt = 29979245800.0 ![cm/s]
    real, parameter :: rt_planck  = 6.62606957e-27
    real, parameter :: PI         = 3.14159265358979323846264338327950288419716939937510582097494

!    real, parameter :: PI         = 
! unitless, lower limit
      ll = lowlimit/(Teff*rt_boltz)

! unitless, upper limit 
      ul = uplimit/(Teff*rt_boltz)

! k^4 T^4/c^2 h^3
!		prefac  = Teff**4d0*2d0*PI*(rt_boltz)**4d0/((rt_speedlt)**2d0*(rt_planck)**3d0)
! k^3 T^3/c^2 h^3
      prefac2 = Teff**3d0*2d0*(rt_boltz)**3d0/((rt_speedlt)**2d0*(rt_planck)**3d0)*PI
! k^3 T^3/c^2 h^3 * sigma0 *x0^3
!    prefac3 = prefac*ll**3d0

! 1000 should be enough, because of the exponential it converges quickly, 10 would probably be also fine
      iter = 1000

! this calculates int (I_v) dv, ll...inf for the photon energy 
! and int ( I_v/hv) dv, ll...inf for the photon emission rate
!  	energy = 0d0
!		rate = 0d0
!		weightedE = 0d0

      tmp     = 0d0

      nPhotll = 0d0
      nPhotul = 0d0

      energyll= 0d0
      energyul= 0d0

      ratell  = 0d0
      rateul  = 0d0

! integrates from lowerlimit -> infinity
! and upperlimit -> infinity
! then calculates value from lowerlimit - upperlimit
    do i = 1, iter
! total unweighted energy in the bin
!		  energy = energy    		+	(exp(-i*ll)*( 6d0/i**4d0  +(6d0*ll)/i**3d0 +(3d0*ll**2d0)/i**2d0 +ll**3d0/i ))
! this is 0 at infinity due to exp(-inf)
  ratell = ratell + (exp(-i*ll)*( 2d0/i**3d0  +(2d0*ll)/i**2d0 +    (ll**2d0)/i  ))

!this is from int 1/x 1/(exp(-x)-1)
    call E1XB(i*ll,tmp)

! cross section weighted number of photons 
    nPhotll  = nPhotll + tmp
! weighted energy bin
    energyll = energyll + (exp(-i*ll)/i)
  enddo

      tmp = 0d0
      if(ul .gt. 0.0) then 
        do i = 1, iter
! total unweighted energy in the bin
!        energy = energy    		+	(exp(-i*ll)*( 6d0/i**4d0  +(6d0*ll)/i**3d0 +(3d0*ll**2d0)/i**2d0 +ll**3d0/i ))
          rateul = rateul + (exp(-i*ul)*( 2d0/i**3d0  +(2d0*ul)/i**2d0 +    (ul**2d0)/i  ))

!this is from int 1/x 1/(exp(-x)-1)
          call E1XB(i*ul,tmp)
! number of photons 
          nPhotul  = nPhotul + tmp
! weighted energy bin
          energyul = energyul + (exp(-i*ul)/i)
        enddo
      endif

! actually computed int L_nu d_nu, L_nu = 4 Pi R_s^2 *Pi* Flux_s
!		phFlux    =  4d0*PI*starRad**2d0*prefac2*rate
!		print*,4d0*PI*starRad**2d0*1.204556410e24
! energy per emitted photon
!		phEnergy  =  4d0*PI*starRad**2d0*prefac*energy/phFlux

! energy per emitted photon
      if(ul .gt. 0.0) then
        phEnergy = Teff*rt_boltz*( (energyll - energyul)/ (nPhotll - nPhotul))
        phNumber = prefac2*(ratell - rateul)
        phActive = (lowlimit/(Teff*rt_boltz))**3d0*prefac2*(nPhotll - nPhotul)*sigma
      else
        phEnergy = energyll/nPhotll*Teff*rt_boltz
        phNumber = prefac2*ratell
        phActive = (lowlimit/(Teff*rt_boltz))**3d0*prefac2*nPhotll*sigma
      endif
   end subroutine getSourceEmission

! assume constant cross section sigma
! should always have upper limit
   subroutine getEmissionSigma(Teff, lowlimit, uplimit, sigma, phEnergy, phNumber)

     implicit none

     real, intent(in)  :: Teff, lowlimit, uplimit, sigma
     real, intent(out) :: phEnergy, phNumber

     integer*8:: i  
     integer ::  iter
     real*8  :: ll, ul, tmp
     real :: stboltz,  energy, kT, prefac, PI2
     real :: prefac2, solarRadius, weightedE, prefac3
     real :: nPhotll, nPhotul, energyll, energyul, ratell, rateul

! hardcoded for portability:
     real, parameter :: rt_boltz   = 1.3806488e-16 
     real, parameter :: rt_speedlt = 29979245800.0 ![cm/s]
     real, parameter :: rt_planck  = 6.62606957e-27
     real, parameter :: PI         = 3.14159265358979323846264338327950288419716939937510582097494

!    real, parameter :: PI         = 
! unitless, lower limit
      ll = lowlimit/(Teff*rt_boltz)

! unitless, upper limit 
      ul = uplimit/(Teff*rt_boltz)

! k^4 T^4/c^2 h^3
!		prefac  = Teff**4d0*2d0*PI*(rt_boltz)**4d0/((rt_speedlt)**2d0*(rt_planck)**3d0)
! k^3 T^3/c^2 h^3
! pi from solid angle integration needed for photon number 
      prefac2 = Teff**3d0*2d0*(rt_boltz)**3d0/((rt_speedlt)**2d0*(rt_planck)**3d0)*PI
! k^3 T^3/c^2 h^3 * sigma0 *x0^3
!    prefac3 = prefac*ll**3d0

! 1000 should be enough, because of the exponential it converges quickly
! less than 100 also oks
      iter = 1000

! this calculates int (I_v) dv, ll...inf for the photon energy 
! and int ( I_v/hv) dv, ll...inf for the photon emission rate
      tmp     = 0d0

      nPhotll = 0d0
      nPhotul = 0d0

      energyll= 0d0
      energyul= 0d0

! integrates from lowerlimit -> infinity
! and upperlimit -> infinity
! then calculates value from lowerlimit - upperlimit
    do i = 1, iter

      energyll = energyll +	(exp(-i*ll)*( 6d0/i**4d0  +(6d0*ll)/i**3d0 +(3d0*ll**2d0)/i**2d0 +ll**3d0/i ))

! cross section weighted number of photons 
    nPhotll  = nPhotll + (exp(-i*ll)*( 2d0/i**3d0  +(2d0*ll)/i**2d0 + (ll**2d0)/i  ))
  enddo

      tmp = 0d0
      if(ul .gt. 0.0) then 
        do i = 1, iter
!this is from int x^3 1/(exp(-x)-1)
          energyul = energyul +	(exp(-i*ul)*( 6d0/i**4d0  +(6d0*ul)/i**3d0 +(3d0*ul**2d0)/i**2d0 +ul**3d0/i ))
!this is from int x^2 1/(exp(-x)-1)

! number of photons 
          nPhotul  = nPhotul + (exp(-i*ul)*( 2d0/i**3d0  +(2d0*ul)/i**2d0 + (ul**2d0)/i  ))
        enddo
      endif

! energy per emitted photon
      if(ul .gt. 0.0) then
        phEnergy = Teff*rt_boltz*( (energyll - energyul)/ (nPhotll - nPhotul))
        phNumber = prefac2*(nPhotll - nPhotul)*sigma
      else
        phEnergy = energyll/nPhotll*Teff*rt_boltz
!        phNumber = prefac2*ratell
      endif
	  end subroutine getEmissionSigma

! piecewise constant sigma0 and v^-3 dropoff cross section
! wrappa function
subroutine directH2(Teff, lowlimit, uplimit, phEnergy, phNumber)
    implicit none

    real, intent(in)  :: Teff, lowlimit, uplimit
    real, intent(out) :: phEnergy, phNumber
    real :: totalE
    real :: totalN
    real,parameter :: eVtoErg = 1.60217657e-12
    real,parameter :: Mb      = 1e-18
    real :: low, up, sigma, tmp
    
!s[nu_] := 
! Piecewise[{{0, 0 <= nu <= 15.2}, {0.09*Mb, 
!    15.2 < nu <= 15.45}, {1.15*Mb, 15.45 < nu <= 15.7}, {3.00*Mb, 
!    15.7 < nu <= 15.95}, {5.00*Mb, 15.95 < nu <= 16.2}, {6.75*Mb, 
!    16.2 < nu <= 16.4}, {8.00*Mb, 16.4 < nu <= 16.65}, {9.0*Mb, 
!    16.65 < nu <= 16.85}, {9.5*Mb, 16.85 < nu <= 17.00}, {9.8*Mb, 
!    17. < nu <= 17.2}, {10.10*Mb, 17.2 < nu <= 17.65}, {9.85*Mb, 
!    17.65 < nu < 18.1}, {9.75*Mb*(18.1/nu)^3, 18.10 < nu}}]

    totalE = 0.0
    totalN = 0.0

! pieces 15.2
    low = lowlimit
    up  = 15.45*eVtoErg
    sigma  = 0.09*Mb 
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 15.45*eVtoErg
    up  = 15.70*eVtoErg
    sigma  = 1.15*Mb
    phEnergy = 0.0
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 15.70*eVtoErg
    up  = 15.95*eVtoErg
    sigma  = 3.00*Mb
    phEnergy = 0.0
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 15.95*eVtoErg
    up  = 16.20*eVtoErg
    sigma  = 5.00*Mb
    phEnergy = 0.0
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 16.20*eVtoErg
    up  = 16.40*eVtoErg
    sigma  = 6.75*Mb
    phEnergy = 0.0
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber
    
    low = 16.40*eVtoErg
    up  = 16.65*eVtoErg
    sigma  = 8.00*Mb
    phEnergy = 0.0
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 16.65*eVtoErg
    up  = 16.85*eVtoErg
    sigma  = 9.0*Mb
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 16.85*eVtoErg
    up  = 17.00*eVtoErg
    sigma  = 9.5*Mb
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 17.00*eVtoErg
    up  = 17.20*eVtoErg
    sigma  = 9.8*Mb
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 17.20*eVtoErg
    up  = 17.65*eVtoErg
    sigma  = 10.10*Mb
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 17.65*eVtoErg
    up  = 18.10*eVtoErg
    sigma  =  9.85*Mb
    call getEmissionSigma(Teff, low, up, sigma, phEnergy, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    low = 18.10*eVtoErg
    up  = uplimit
    sigma = 9.75*Mb
    call getSourceEmission(Teff, low, up, sigma, phEnergy, tmp, phNumber)
    totalE = totalE + phEnergy-low
    totalN = totalN + phNumber

    phNumber = totalN
    phEnergy = totalE + lowlimit
	end subroutine directH2

! some digits more accurate, but more expensive to calculate
    SUBROUTINE E1XB(X,E1)
      implicit none
!C
!C       ============================================
!C       Purpose: Compute exponential integral E1(x)
!C       Input :  x  --- Argument of E1(x)
!C       Output:  E1 --- E1(x)  ( x > 0 )
!C       ============================================
!C
      REAL*8 :: A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,Y,Z,GA,T0
      INTEGER :: K,M
      REAL*8, INTENT(IN)  :: X
      REAL*8, INTENT(OUT) :: E1 
      IF (X.EQ.0.0) THEN
         E1=1.0D+300
      ELSE IF (X.LE.1.0) THEN
         E1=1.0D0
         R=1.0D0
         DO 10 K=1,25
            R=-R*K*X/(K+1.0D0)**2
            E1=E1+R
            IF (DABS(R).LE.DABS(E1)*1.0D-15) GO TO 15
10       CONTINUE
15       GA=0.5772156649015328D0
         E1=-GA-DLOG(X)+X*E1
      ELSE
         M=20+INT(80.0/X)
         T0=0.0D0
         DO 20 K=M,1,-1
            T0=K/(1.0D0+K/(X+T0))
20       CONTINUE
         T=1.0D0/(X+T0)
         E1=DEXP(-X)*T
      ENDIF
      RETURN
    END SUBROUTINE
!-------------------------------------------------------------------------------  
end module pt_sourceUtil
