!!  ++ / 
!!  + / -   Heating/Cooling
!!   / --
!!written by C. Baczynski, 2012-2013

!! Description:
!!   uses heating rates from SN and ionising radiation to 
!!   calculate new temperature
!!   cooling according to Dalgarno&McCray 1972 is applied
!!
!! Input: 
!!   dt: current simulation timestep
!!   blockcount: number of local blocks
!!   blockList: list of local block IDs
!!   time: global simulation time
!!
!! Modifications:
!!  12/10/14 JCIM
!!
!!          modified heating and cooling for temperatures below 1000 K 
!!
!!			Temps between 1000 - 100 K: interpolation between the DalgarnoMacCray curve and the Parameterized curve.
!!										weak dependence with the local volume density
!!
!!			Below 100 K: Simon's parrameterized cooling curve
!!						 Strongly dependent on the local volume density
!!
!!			Heating: accounts for dust shielding, cosmic rays heating and variations in the metallicity.
!!					 now is density dependent.   
!!			Cooling: accounts CII excitation assuming collisions with HI dominate
!!					 conversion from C+ to CO assuming all CO in LTE and neglects opacity, freeze-out.
!!


subroutine Heat (blockCount,blockList,dt,time)

  use Heat_data	
  use Simulation_data, ONLY : sim_Z, sim_G0, sim_pe_h, sim_cr_h, sim_constant_heating, &
							  sim_stratify_heating
							
  use Grid_interface,		ONLY : Grid_fillGuardCells, Grid_getCellCoords, &
											         Grid_releaseBlkPtr,Grid_getBlkPtr, Grid_getDeltas, &
															 Grid_getBlkIndexLimits

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Eos_interface, 		ONLY : Eos_wrapped
  use Driver_data, 			ONLY : dr_simTime, dr_nStep

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Flash_mpi.h"

	! arguments  
	integer,intent(IN)	:: blockCount
	integer,dimension(blockCount),intent(IN)	:: blockList
	real,intent(IN)			:: dt,time

	! block data
  	integer	:: blockID, thisBlock
  	real, pointer, dimension(:,:,:,:)	:: solnData
  	real, allocatable, dimension(:)		:: xCoord, yCoord, zCoord
  	real, allocatable,dimension(:)		:: dx, dy, dz
  	integer														:: xSizeCoord, ySizeCoord, zSizeCoord
  	logical														:: getGuardCells = .true.
  	integer, dimension(2,MDIM)				:: blkLimits, blkLimitsGC
  	real, dimension(MDIM)							:: del

	! iterators
	! TODO see if they can be reused
  	integer :: i, j, k, l, m
  	integer :: Ncycles
  	real    :: dtsub

	! variables for the SN explosion
  	real    :: xx, yy, zz

	! heating and cooling variables
  	real :: tranheat, sdot, sheat,tdepheat
  	real :: tmp, rho, ei, ek, timestep, conf
  	real :: dt0, dt1, dt_dei, ndens, TtoEI

  	real :: radia, tcool, theat, ttherm, scdot, dei
  	real ::  convf, tmpnew, eplus, eiold, tmpold, emin
  	integer :: nstep

	real :: myheat, new_G0
	
	real :: coolhere, coolold, cooly0, cooly1, cooly2
	
	real :: emin_int, emin_new, emin_100, emin_1000
	
	real :: x0, x1, x2, y0, y1, y2
	real :: t0, t1, t2
	

  	if (.not. he_useHeat) return

  	if ( dt .lt.  he_dtThres) then
   		
  		call Timers_start("heat")
		
		if (he_meshMe == MASTER_PE ) then
			print*,'turning off cooling for this timestep, dt < dtThres'
		endif

		do thisBlock = 1, blockCount
			blockID = blockList(thisBlock)
			! Get a pointer to solution data 
			call Grid_getBlkPtr(blockID,solnData)
			do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
				do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
					do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
						solnData(PHHE_VAR,i,j,k) = 0d0
					enddo ! coord loops i
				enddo	! coord loops j
			enddo	! coord loops k

			call Grid_releaseBlkPtr(blockID,solnData)
		enddo ! block loop

		call Timers_stop("heat")
		return
	endif
	

  	! start the timer ticking
  	call Timers_start("heat")

	! loop over local blocks in domain, so the blocklist and apply heating and cooling as well as stellar winds
	! this loop should go into its own subroutine 
	do thisBlock = 1, blockCount
		blockID = blockList(thisBlock)

		! Get a pointer to solution data 
		call Grid_getBlkPtr(blockID,solnData)
		call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) !indices for the interior zones and all zones including guard zones
		call Grid_getDeltas(blockID, del) !grid spacing dx dz dy

		! Get a pointer to solution data 
		xSizeCoord = blkLimitsGC(HIGH,IAXIS)
		ySizeCoord = blkLimitsGC(HIGH,JAXIS)
		zSizeCoord = blkLimitsGC(HIGH,KAXIS)

		! allocate space for dimensions
		allocate(xCoord(xSizeCoord))
		allocate(yCoord(ySizeCoord))
		allocate(zCoord(zSizeCoord))
		allocate(dx(xSizeCoord))
		allocate(dy(ySizeCoord))
		allocate(dz(zSizeCoord))

	  	! actually only one dx would be needed, -SizeCoord might not be right size, for cubic zones irrelephant
		dx(:) = del(IAXIS)
		dy(:) = del(JAXIS)
		dz(:) = del(KAXIS)

	  	call Grid_getCellCoords(IAXIS,blockID,CENTER,getGuardCells,xCoord,xSizeCoord)
	  	call Grid_getCellCoords(JAXIS,blockID,CENTER,getGuardCells,yCoord,ySizeCoord)
	  	call Grid_getCellCoords(KAXIS,blockID,CENTER,getGuardCells,zCoord,zSizeCoord)     

	  	! loop over all zones in block
    	! for 2d k is just 1
		do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
			zz = zCoord(k)
			do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
				yy = yCoord(j)
				do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
					xx = xCoord(i)

					!  get zone variables
					tmp    = solnData(TEMP_VAR,i,j,k)
					rho    = solnData(DENS_VAR,i,j,k)
					ei     = solnData(EINT_VAR,i,j,k)
					eiold  = ei
					tmpold = tmp

					! ei to T conversion, if eos was called both should be proportional
					TtoEI  = ei/tmp

					! not sure ek is needed
					ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 + &
			     	    		& solnData(VELY_VAR,i,j,k)**2 + &
										& solnData(VELZ_VAR,i,j,k)**2)

					! get heating from feedback processes, i.e. radiation, SN
					! this is in [erg/(s cm^3)], as radiation heating only sees atomic Hydrogen
					eplus  = solnData(PHHE_VAR,i,j,k)

					! convert to erg /(g s)
					eplus  = eplus/rho

					! number density of the gas, all species, heating applied to all of them equally
					! heating rate should have taken specific species into account, i.e. radiation
					! this hould be hydrogen number density
					ndens  = rho/(he_abar*he_protonmass)

					! conversion factor for cooling rate 
					conf   = ndens*ndens/rho

					! add UV heating
					if ( (tmp < he_theatmax) .and. (tmp > he_theatmin) ) then
						
						! add any additional heating terms
						! photoelectric heating
						! stratify ?
						! dust shielding and cosmic ray heating.

						if (sim_constant_heating) then ! switch old - new heating rates
							
								if(he_stratifyHeat) then
									! converts to [erg/(g s)] 
									! he_peheat is in [erg/s]
									eplus = eplus +  he_peheat * exp(-abs(zz)/(he_h_UV)) / (rho/ndens)
								else
									eplus = eplus +  he_peheat / (rho/ndens)
								endif ! heat stratification
  						else
							
	   						call Parameterized_heating(zz, rho, sim_Z, he_peheat, sim_stratify_heating, sim_pe_h, sim_cr_h, myheat)						
  							eplus = eplus + myheat
  							
  						endif ! switch old - new heating rates
						
					endif	! allowed range of temperatures

					! no heating is done if below he_absTmin, so turn it to 0d0 
					!				  if (tmpold .gt. he_absTmin) then
					
					! checks for density and temperature range allowed to radiatively cool
					! heatmin/max should be consistent with tradmin/max
					! no heating if outside range
					if (  (tmp   <= he_tradmax .AND. tmp   >= he_tradmin) & 
					.AND. (ndens >= he_dradmin .AND. ndens <= he_dradmax) ) then
							
		 				ei  = eiold
				      	tmp = tmpold
				      	dt0 = dt
				      	dt1 = dt
				
						! absolute subcycle time goes from 0.0 to hydro dt
				      	timestep = 0.0
				
						! number of subcycles, upper limit could be implemented to switch to implicit solution 
						! but usually only few zones require a very large number (> 1e5) steps
				      	nstep = 0
				
						! Temperature values where old cooling curve is implemented ( t2 < t )
						! Interpolation between the old and new cooling curves      ( t1 < t < t2 )
						! Temperature below which the new cooling curve is used     (      t < t1 )
						t0 = 100.0
						t1 = 5000.0
						t2 = 10000.0

				     	do
				        nstep = nstep + 1

							! get da cooling rate [erg cm^3/s] (volumetric)
                			if(he_coolOff) then
                  				emin = 0d0
                			else
								! JC: call the standard cooling curve implemented by Joung & Mac Low 2006.
								if (tmp .gt. t1) then
	  			          			call radloss(tmp,emin)
									! convert to erg/(g s)
									emin = emin*conf
								else
									! JC: If the temperature is between 100 - 5000 Kelvin Interpolate the IMG and Ryan's
									! 	  cooling curves. Weak dependence in the density.
									if (tmp .gt. t0) then
										
										call radloss(tmp,coolold) 
						                coolold  = coolold  * conf

						                call radloss(t1,cooly1) 
						                cooly1 = cooly1 * conf

						                call radloss(t2,cooly2)
						                cooly2  = cooly2  * conf
										
										call low_temp_cooling(rho, t0, sim_Z, cooly0)
										!call low_temp_cooling(rho, tmp, sim_Z, coolnew)

						                call second_order_Interpolation(t0, t1, t2,cooly0, cooly1, cooly2, tmp, coolhere)
						                !coolhere = 10**coolhere

						                if (coolold .lt. coolhere) then
											emin = coolold
						                else
											emin = coolhere										
										end if
									else ! JC: Temperatures between 10 - 100
										
										call radloss(tmp,coolold) 
						                coolold  = coolold  * conf

										call low_temp_cooling(rho, tmp, sim_Z, coolhere)

						                if (coolold .lt. coolhere) then
											emin = coolold
						                else
											emin = coolhere
										end if
										
									end if ! 10 < T < 1000
								endif ! T > 1000
                			endif ! he_coolOff switch.
							
							! calculate energy change rate, at 1e4 K rapid change might be unstable if temperature is just above
							! with a much shallower slope
       						dei    = max(1e-50,abs(emin-eplus))

							! get timestep, large change -> small time step
							dt_dei = ei/abs(dei)
							
							! subcycle timestep size 
        					dt1    = he_subfactor*dt_dei

							! last step should land on hydro dt
        					dt0    = min(dt1,dt-timestep)
							
							! change internal energy by a fraction
							ei     = ei + dt0*(eplus-emin)
							
							! change temperature 
							tmp 	 = ei/TtoEI

							! over cooled, fix	
			        		if (tmp .lt. he_absTmin) then
								! fix internal energy to 10 K equivalent
								ei  = he_absTmin*TtoEI
			          			tmp = he_absTmin
			        		endif

							! update current absolute subcycle time
				        	timestep = timestep + dt0

							! reached hydro dt, step out
				        	if (abs(timestep-dt) .lt. 1.0e-6*dt) exit

						enddo
					else
						
						! just heat, might be useful if just uv heating is on
						! might lead to very high temperatures if if is just below T threshold
						! but large energy input
						ei  = ei + eplus
					endif ! end of checks for density and temperature in allowed ranges. 


				! print warning if temperature too high, maybe fix to upper limit

          		if (tmp > he_absTmax) then

					! write to file some info
					open(he_funit_log, file=trim(he_outfile), position='APPEND')
					write(he_funit_log,'(3(1X,ES16.9),6(1X,I16))') &
					dr_simTime, tmp, he_absTmax, dr_nStep, i, j, k, blockID, he_meshMe
					close(he_funit_log)

					! fix it and call eos with density temperature mode, afterwards gracefully exit
					ei  = he_absTmax*TtoEI
            		tmp = he_absTmax
					! adjust internal energy

#ifdef verbose
           			print*,'\\Correction to temp =', he_absTmax,' in Heat'
           			print*,'\\purely from shock heating'
           			print*,'\\core ID, block ID, zone ID xyz',he_meshMe,blockID,li,j,k
#endif
          		endif ! end if gas is too hot.

				! floor internal energy
				if(ei .lt. he_smallpres/solnData(DENS_VAR,i,j,k)) then
					ei =  he_smallpres/solnData(DENS_VAR,i,j,k)
				endif

				! update the global thermodynamic quantities due to the phen. heating
				solnData(ENER_VAR,i,j,k) = ei + ek
				solnData(EINT_VAR,i,j,k) = ei
!					solnData(TEMP_VAR,i,j,k) = tmp

				!clean up heating rate
				solnData(PHHE_VAR,i,j,k) = 0d0

				enddo ! coord loops
			enddo
		enddo

		!  crank changed state variables through EOS
		call Eos_wrapped(MODE_DENS_EI, blkLimitsGC, blockID)
!		call Eos_wrapped(MODE_DENS_TEMP, blkLimitsGC, blockID)

		!  clean up memory 
		call Grid_releaseBlkPtr(blockID,solnData)
		deallocate(xCoord)
		deallocate(yCoord)
		deallocate(zCoord)
		deallocate(dx)
		deallocate(dy)
		deallocate(dz) 	  
	enddo ! block loop

	call Timers_stop("heat")

	return
end subroutine Heat


subroutine Linear_Interpolation(x0, x1, y0, y1, xin, yout)
	
	implicit none
	
	real, intent(IN) :: x0, x1, y0, y1, xin
	real, intent(OUT) :: yout
	real :: a0, a1
	
	a0 = (xin - x1) / (x0 - x1)
	a1 = (xin - x0) / (x1 - x0)
	
    yout = y0 * a0 + y1 * a1

    return 
end subroutine Linear_Interpolation


subroutine second_order_Interpolation(x00, x11, x22, y00, y11, y22, xinin, youtout)
	
	implicit none
	
	real, intent(IN) :: x00, x11, x22, y00, y11, y22, xinin	
	real, intent(OUT) :: youtout
	
	real :: x0, x1, x2, y0, y1, y2, xin
	real :: a0, a1, a2, yout
	
	x0 = LOG10(x00)
	x1 = LOG10(x11)
	x2 = LOG10(x22)
	y0 = LOG10(y00)
	y1 = LOG10(y11)
	y2 = LOG10(y22)
	xin = LOG10(xinin)

    a0 = ((xin - x1) * (xin - x2)) / ((x0 - x1 ) *(x0 - x2))   
    a1 = ((xin - x0) * (xin - x2)) / ((x1 - x0) * (x1 - x2))   
    a2 = ((xin - x0) * (xin - x1)) / ((x2 - x0) * (x2 - x1))
        
    yout = y0*a0 + y1*a1 + y2*a2

	youtout = 10.0**yout
    
    return 
end subroutine second_order_Interpolation
