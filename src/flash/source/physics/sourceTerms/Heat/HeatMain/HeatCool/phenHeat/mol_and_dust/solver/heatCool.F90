!!****ih* source/physics/sourceTerms/Heat/HeatMain/HeatCool/phenHeat/mol_and_dust/solver/heatCool.F90
!!
!! NAME
!!
!!  heatCool
!!
!! SYNOPSIS
!!
!!  use heatCool
!!  call heating_and_cooling (
!!            real(INOUT) :: ei,
!!            real(IN)    :: timestep,
!!            real(INOUT) :: temp,
!!            real(OUT)   :: ncount)
!!
!! DESCRIPTION
!!
!!  Single subroutine call for all
!!  heating and cooling implemented
!!  in Flash for star formation and feedback.
!!
!!  Joshua Wall, Drexel University 2017
!!
!! ARGUMENTS
!!
!!  ei       : specific internal energy (erg/gram)
!!  timestep : time step over which to compute heating/cooling (sec)
!!  temp     : temperature (Kelvin)
!!  ncount   : number of sub-steps taken to compute heating/cooling
!!
!!***

module heatCool

implicit none

#include "Flash.h"

contains

subroutine heating_and_cooling(ei, timestep, temp, ncount)

!#define ONE_CELL_TESTING
!#define write_data
!#define timing

use OdeInt
!use heatCool
use cool_vars
use Heat_data
use Cool_data, only : he_int_method, mp, useDustCool
use Driver_interface, only : Driver_abortFlash

implicit none

!integer, parameter  :: dp = KIND(1d0)

real(dp), intent(in)             :: timestep
real(dp), intent(inout)          :: ei, temp

real(dp) :: eplus, emin, PeHeat, dust_heat, dust_cool, mol_cool

real(dp) :: eiold, tempold, ei_init, temp_init, eiolder, ei1, ei2
real(dp) :: dei, dei2, dt0, dt1, total_time, c1
real(dp) :: total_local_FUV, dust_cr_heat

logical  :: bracket_found, cycled, first_loop

!real(dp), parameter :: twoEV=2.0_dp*1.6022d-12
!real(dp), parameter :: min_temp = 10.0_dp
real(dp), parameter :: err_tol = 1d-6

#ifdef timing
real(dp) :: t1, t2
#endif

integer, intent(out)  :: ncount
integer, parameter :: max_count = 100000

#ifdef write_data
character(len=*), parameter    :: frmt='(es15.3,t20,es15.3)'
#endif

!character(len=*), parameter :: method='euler'


!interface
!    function dei_dt(del_t, ei_in)
!    implicit none
!    integer, parameter    :: dp=KIND(1d0)
!    real(dp), intent(in) :: del_t, ei_in
!    real(dp)             :: dei_dt
!    end function dei_dt
!end interface

total_time = 0.0d0

ndens   = rho / (mu_mol * mp)
TtoEI   = ei / temp
conf    = ndens*ndens/rho



#ifdef write_data
write(*,'(A,ES13.3E3)') "Initial ndens =", ndens
write(*,'(A,ES13.3E3)') "Initial conf =", conf
write(*,'(A,ES13.3E3)') "Initial ephen =", ephen
write(*,'(A,ES13.3E3)') "Initial TtoEI =", TtoEI

write(*,'(A,ES13.3E3)') "Initial energy =", ei
write(*,'(A,ES13.3E3)') "Initial temp =", temp

open(unit=10, file='euler.dat')

write(10,'(2ES13.3E3)') total_time, temp
#endif

! Euler's Method
if (he_int_method .eq. 'Euler') then

emin = 0d0
eplus = 0d0
PeHeat = 0d0
dust_heat = 0d0
dust_cool = 0d0
mol_cool = 0d0

ncount = 0

#ifdef timing
    call cpu_time(t1)
#endif

    do 
        eiolder = eiold
        eiold   = ei
        tempold = temp

        !call heating(temp, ndens, rho, ephen, zpos, eplus)
        !call cooling(temp, ndens, conf, emin)
        

        dei  = dei_dt(timestep, eiold)
        dt1  = he_subfactor*ei / max(abs(dei),1d-50)
        dt0  = min( (timestep - total_time), dt1)

#ifdef debug    
        write(*,'(A,ES13.3E3)') "Euler emin =", emin
        write(*,'(A,ES13.3E3)') "Euler dt1 =", dt1
#endif    
        ei   = ei + dei*dt0

! Attempt to dampen oscillations by averaging with the past solution.
! Note this makes the solution lag behind the true solution, so no good.
        !ei = 0.5*(eiold + ei)
    !ei = 1.0d0/3.0d0*(eiolder + eiold + ei)

        temp = ei / TtoEI

        total_time = total_time + dt0
            
#ifdef write_data
        write(10,frmt) total_time, temp
#endif

#ifdef debug    
        write(*,'(A,ES13.3E3)') "Euler energy =", ei
        write(*,'(A,ES13.3E3)') "Euler temp =", temp

        write(*,'(A,ES13.3E3)') "Euler dt =", dt0
        write(*,'(A,ES13.3E3)') "Euler total time =", total_time
#endif

        ncount = ncount + 1
        
        if ( ncount > max_count ) then
            call Driver_abortFlash("[heatCool]: Euler exceeded max steps!")
        end if
        
        if ( (timestep - total_time) < 1d-6) exit
        
    end do

! For stability, this should always be outside the do loop. Let the 
! integration converge, then set the min temp.
    if (temp .lt. he_absTmin) then
        temp = he_absTmin
        ei   = he_absTmin*TtoEI
    end if

        ! If we are tracking both heating TO dust from PE heating and cooling OF dust
        ! from collisions with hydrogen, we can also track the dust temperature
        ! consistently.
        !
        ! OK to stop tracking dust temperature when T > T_sputter because
        ! calculation assumes thermal equilibrium; history is not needed.
        ! -AT, 2019 May 15
    if (useDustCool .and. temp .le. he_dust_sputter_temp) then
        tdust = get_dust_temperature(tdust, temp, ndens, 1.0d0/mu_mol, dust_heat)
    end if

#ifdef timing
    call cpu_time(t2)
#endif

#ifdef write_data
    close(10)

    write(*,'(A,ES13.3E3)') "Euler final energy =", ei
    write(*,'(A,ES13.3E3)') "Euler final temp =", temp
    write(*,'(A,I8)') "Euler num steps =", ncount
    write(*,'(A,ES13.3E3)') "Euler total time =", total_time
    write(*,'(A,ES13.3E3,A)') "Euler took =", t2-t1, "secs"
#endif

else if (he_int_method .eq. "Implicit") then

!#define debug
!#define debug_bracket

emin = 0d0
eplus = 0d0
PeHeat = 0d0
dust_heat = 0d0
dust_cool = 0d0
mol_cool = 0d0

ncount = 0
cycled = .false.
first_loop = .true.

#ifdef write_data
    open(unit=10, file='beuler.dat')
    write(10,'(2ES13.3E3)') total_time, temp
#endif

#ifdef ONE_CELL_TESTING
        write(*,'(A,ES13.3E3)') "Start energy =", ei
        write(*,'(A,ES13.3E3)') "Start temp =", temp
        write(*,'(A,ES13.3E3)') "Start dust temp =", tdust
        write(*,'(A,ES13.3E3)') "[heatcool]: Start ion =", xHp
        call flush(6)
#endif

#ifdef timing
    call cpu_time(t1)
#endif
    do while ( (timestep - total_time) .gt. (1d-6*timestep))

        if (cycled) print*, "We cycled already!!!"

        ncount = ncount + 1
        
        eiolder = ei
        eiold   = ei
        tempold = temp
        

        dei  = dei_dt(timestep, eiold)
        
        if (first_loop) then
          !dt1  = he_subfactor*eiold / max(abs(dei),1d-50)
          !dt0  = min( (timestep - total_time), dt1)
          !first_loop = .false.
          dt0 = timestep
        else
          print*, "Why are you here?!"
        end if
    ! Take a trial step with Euler.

#ifdef ONE_CELL_TESTING
        write(*,'(A,ES13.3E3)') "[heatCool]: dei =", dei
        write(*,'(A,ES13.3E3)') "[heatCool]: dt0 =", dt0
        write(*,'(A,ES13.3E3)') "[heatCool]: dt1 =", dt1
#endif
        ! Okay, what we want here is to get just off the
        ! current solution eiold, without any chance of
        ! overshooting into the negatives. Note negative
        ! solutions are going to probably do strange things
        ! when the implicit solver calls all the heating and
        ! cooling fuctions with negative temperatures.
        
        ! zbrac will properly expand a slight difference
        ! in the solution to cover the zero in the solution.
        ! So we'll get a slight difference in the right direction
        ! with a tiny Euler step.

        eiold = eiolder + dei*dt0*1.0d-3

        
        !dei  = dei_dt(0.0d0, eiold)        
        !dt1  = he_subfactor*eiold / max(abs(dei),1d-50)
        !dt0  = min( (timestep - total_time), dt1, dt0)
        
    ! Now use the trial solution and the old solution as a first guess
    ! to bracket the real solution.
    
        ei1 = eiolder
        ei2 = eiold
        bracket_found = .false.
        ! If the temperature has bottomed out, dei may be zero.
        ! This means 1) the bracket method won't see a difference in the
        ! answers and 2) the root finding method is overkill anyway.
        ! Check and see if dei/dt ~ 0, and if so just take an Euler step.
        if (abs((ei1-ei2)/ei1) .lt. (err_tol*dt0/(timestep-total_time))) then
        
#ifdef debug_bracket
!                print*, "[heatCool]: Warning! ei1 ~ ei2."
!                write(*,'(A,ES13.3E3)') "Temp=", temp
!                write(*,'(A,ES13.3E3)') "dei=", dei
!                write(*,'(A,ES13.3E3)') 'dt=', dt0
!                call flush(6)
#endif
            dt0 = (timestep - total_time)
            ei = eiolder + dei*dt0
            temp = ei / TtoEI

            ! If we are tracking both heating TO dust from PE heating and cooling OF dust
            ! from collisions with hydrogen, we can also track the dust temperature
            ! consistently.
            !
            ! OK to stop tracking dust temperature when T > T_sputter because
            ! calculation assumes thermal equilibrium; history is not needed.
            ! -AT, 2019 May 15
            if (useDustCool .and. temp .le. he_dust_sputter_temp) then
              ! Total local FUV flux is stellar FUV + background Habing field.
              total_local_FUV = Gflux + he_Gzero*f_ext(rho, ndens, temp)*strat_factor
              ! The heating of dust due to FUV (the energy that DIDN'T eject PEs.
              call pe_heating(total_local_FUV, xHp, ndens, temp, rho, PeHeat, dust_heat)
              ! Now add heating from cosmic rays on the dust.
              call heating(temp, ndens, rho, 0.0d0, strat_factor, dust_cr_heat)
              ! convert from [ergs/g/s] to [ergs/cm^3/s]
              dust_cr_heat = dust_cr_heat * rho
              dust_heat = dust_heat + dust_cr_heat
              ! All the terms in the dust solver have this normalization.
              dust_heat = dust_heat / (ndens*1d-33) 
              tdust = get_dust_temperature(tdust, temp, ndens, 1.0d0/mu_mol, dust_heat)
            end if

            total_time = total_time + dt0
            cycled = .true.
            cycle
            !ei2 = ei1 + ei1*1.5d0
        end if

        call zbrac(dei_dt,ei1,ei2,eiolder,dt0,bracket_found)
        if (.not. bracket_found) then 
                write(*,*) "[heatcool]: Warning! Bracket not found in zbrac."

#ifdef debug_bracket 
                write(*,'(A,ES13.3E3)') "Temp after zbrac=", temp
                write(*,'(A,ES13.3E3)') "ei1 after zbrac=", ei1
                write(*,'(A,ES13.3E3)') "ei2 after zbrac=", ei2
                write(*,'(A,ES13.3E3)') "temp1 after zbrac=", ei1/TtoEI
                write(*,'(A,ES13.3E3)') "temp2 after zbrac=", ei2/TtoEI
                write(*,'(A,ES13.3E3)') "dei after zbrac=", dei
                write(*,'(A,ES13.3E3)') "dt0 after zbrac=", dt0
                write(*,'(A,ES13.3E3)') "dt1 after zbrac=", dt1
                call flush(6)
#endif
        end if

        if (ei1 < 0.0d0 .or. ei2 < 0.0d0) then
                ei1 = max(ei1, TtoEI*1.0d0)
                ei2 = max(ei2, TtoEI*1.0d0)
#ifdef debug_bracket 
                print*, "[heatCool]: Warning! ei1 or ei2 < 0 after bracket."
                !write(*,'(A,ES13.3E3)') "Temp=", temp
                write(*,'(A,ES13.3E3)') "ei1 after fix=", ei1
                write(*,'(A,ES13.3E3)') "ei2 after fix=", ei2
                !write(*,'(A,ES13.3E3)') "dei=", dei
                !write(*,'(A,ES13.3E3)') "dt0=", dt0
                !write(*,'(A,ES13.3E3)') "dt1=", dt1
                call flush(6)
                !stop
#endif
        end if

#ifdef debug
        ! Test bracketing
        !ei2 = ei1 + 0.1*ei1
        
         !write(*,'(A,2ES13.3E3,X,1L)') "Before bracketing", ei1, ei2, bracket_found
         !call zbrac(dei_dt,ei1,ei2,eiolder,dt0,bracket_found)
         write(*,'(A,2ES13.3E3)') "After bracketing and maxing to 1 K ei1, ei2=", ei1, ei2
        
         write(*,'(A,ES13.3E3)') "Implicit Euler eiolder =", eiolder
         write(*,'(A,ES13.3E3)') "Implicit Euler eiold =", eiold
         write(*,'(A,ES13.3E3)') "Implicit Euler temp =", temp
         !write(*,'(A,ES13.3E3)') "Implicit Euler eplus =", eplus
         !write(*,'(A,ES13.3E3)') "Implicit Euler PeHeat =", PeHeat
         !write(*,'(A,ES13.3E3)') "Implicit Euler emin =", emin
         !write(*,'(A,ES13.3E3)') "Implicit Euler mol_cool =", mol_cool
         !write(*,'(A,ES13.3E3)') "Implicit Euler dust_cool =", dust_cool
         write(*,'(A,ES13.3E3)') "Implicit Euler dei =", dei
         write(*,'(A,ES13.3E3)') "Implicit Euler dt1 =", dt1
         write(*,'(A,ES13.3E3)') "Implicit Euler dt0 =", dt0
#endif
        
    ! Now pass the old solution and the trial solution to a root finder.
    ! Accuracy requested is for an delta energy within 1 K. 
    ! now is: TtoEI*1.0d0; previously was: eiold*1d-4
        
        !call secant_method(dei_dt, dt0, eiolder, eiold, max(TtoEI*1.0d0,eiold*1d-6), ei)
        !call falsepos_method(dei_dt, dt0, eiolder, eiold, max(TtoEI*1.0d0,eiold*1d-6), ei)
        ! Ridder's method.
        !ei = ridder ( eiolder, eiold, max(TtoEI*1.0d0,eiold*1d-2), dei_dt, dt0 )

        ! Brent's method: zero(left_brac, right brac, y_old, mach eps, rel tol, dy/dx, dx)
        if (bracket_found) then
          ei = zero ( ei1, ei2, eiolder, 3d-16,  max(TtoEI*1.0d0,eiold*err_tol), dei_dt, dt0 )
        else ! Try Brent's method with absolute limits for temperature of 1 K and 1e9 K.
          ei = zero ( TtoEI*1.0d0, TtoEI*1.0d9, eiolder, 3d-16,  max(TtoEI*1.0d0,eiold*err_tol), dei_dt, dt0 )
        end if
        
#ifdef debug    
        write(*,'(A,ES13.3E3)') "Implicit Euler dei =", dei
        write(*,'(A,ES13.3E3)') "Implicit Euler dt1 =", dt1
#endif    

        temp = ei / TtoEI

        ! If we are tracking both heating TO dust from PE heating and cooling OF dust
        ! from collisions with hydrogen, we can also track the dust temperature
        ! consistently.
        !
        ! OK to stop tracking dust temperature when T > T_sputter because
        ! calculation assumes thermal equilibrium; history is not needed.
        ! -AT, 2019 May 15
        if (useDustCool .and. temp .le. he_dust_sputter_temp) then
            ! Total local FUV flux is stellar FUV + background Habing field.
            total_local_FUV = Gflux + he_Gzero*f_ext(rho, ndens, temp)*strat_factor
            ! The heating of dust due to FUV (the energy that DIDN'T eject PEs.
            call pe_heating(total_local_FUV, xHp, ndens, temp, rho, PeHeat, dust_heat)
            ! Now add heating from cosmic rays on the dust.
            call heating(temp, ndens, rho, 0.0d0, strat_factor, dust_cr_heat)
            ! convert from [ergs/g/s] to [ergs/cm^3/s]
            dust_cr_heat = dust_cr_heat * rho
            dust_heat = dust_heat + dust_cr_heat
            ! All the terms in the dust solver have this normalization.
            dust_heat = dust_heat / (ndens*1d-33) 
            tdust = get_dust_temperature(tdust, temp, ndens, 1.0d0/mu_mol, dust_heat)
        end if

        

#ifdef debug    
        write(*,'(A,ES13.3E3)') "Implicit Euler dei =", dei
        write(*,'(A,ES13.3E3)') "Implicit Euler dt1 =", dt1
#endif    
        
        ! If for some reason we were going to actually be accuracy based,
        ! we'd probably want to base timesteps on the rate
        ! of change in the energy. But currently, we just assume
        ! we are "close enough" with taking one single step, and
        ! this should actually be used. 
        dt1  = he_subfactor*eiolder / max(abs(ei-eiolder),1d-50)
        total_time = total_time + dt0
        dt0  = min( (timestep - total_time), dt1)
            
#ifdef write_data
        write(10,frmt) total_time, temp
#endif

#ifdef debug    
        write(*,'(A,ES13.3E3)') "Implicit Euler energy =", ei
        write(*,'(A,ES13.3E3)') "Implicit Euler temp =", temp

        write(*,'(A,ES13.3E3)') "Implicit Euler dt =", dt0
        write(*,'(A,ES13.3E3)') "Implicit Euler total time =", total_time
#endif

!        if ( ncount > max_count ) then
            !call Driver_abortFlash("[heatCool]: Euler exceeded max steps!")
            ! Just warn for now...
!            print*, "[heatCool]: WARNING WARNING WARNING! Implicit solver &
!                     exceeded max steps, temp probably wrong!"
!        end if
        
        
        !if ((timestep - total_time) < 1d-6) exit
    
    end do

! For stability, this should always be outside the do loop. Let the 
! integration converge, then set the min temp.
    if (temp .lt. he_absTmin) then
        temp = he_absTmin
        ei   = he_absTmin*TtoEI
    end if

#ifdef ONE_CELL_TESTING
        write(*,'(A,ES13.3E3)') "[heatcool]: End energy =", ei
        write(*,'(A,ES13.3E3)') "[heatcool]: End temp =", temp
        write(*,'(A,ES13.3E3)') "[heatcool]: End dust temp =", tdust
        write(*,'(A,ES13.3E3)') "[heatcool]: End ion =", xHp
        call flush(6)
        !stop
#endif

#ifdef timing
    call cpu_time(t2)
#endif
#ifdef write_data
    close(10)

    write(*,'(A,ES13.3E3)') "Backward Euler final energy =", ei
    write(*,'(A,ES13.3E3)') "BackwardEuler final temp =", temp
    write(*,'(A,I8)') "Backward Euler num steps =", ncount
    write(*,'(A,ES13.3E3)') "Backward Euler total time =", total_time
    write(*,'(A,ES13.3E3,A)') "Backward Euler took =", t2-t1, "secs"
#endif

return

else if (he_int_method .eq. 'trap') then

emin = 0d0
eplus = 0d0
PeHeat = 0d0
dust_heat = 0d0
dust_cool = 0d0
mol_cool = 0d0

ncount = 0

! Trapizoid rule
#ifdef write_data
    open(unit=10, file='trap.dat')
    write(10,'(2ES13.3E3)') total_time, temp
#endif

#ifdef timing
    call cpu_time(t1)
#endif

    do

        eiold   = ei
        tempold = temp

    ! First value of y' = f(y0)    
        dei  = dei_dt(timestep, eiold)
        dt1  = he_subfactor*ei / max(abs(dei),1d-50)
        dt0  = min( (timestep - total_time), dt1)
        dt0  = dt0
        
        ei1  = eiold + dei*dt0
        
    ! Second value of y' = f(y1)
        dei2 = dei_dt(timestep, ei1)
        
    ! Average of the two is the integration using trap rule.
        ei2  = eiold + 0.5*dt0*(dei + dei2)
        temp = ei2 / TtoEI

        total_time = total_time + dt0

#ifdef write_data
        write(10,frmt) total_time, temp
#endif

        ! If we are tracking both heating TO dust from PE heating and cooling OF dust
        ! from collisions with hydrogen, we can also track the dust temperature
        ! consistently.
        !
        ! OK to stop tracking dust temperature when T > T_sputter because
        ! calculation assumes thermal equilibrium; history is not needed.
        ! -AT, 2019 May 15
        if (useDustCool .and. temp .le. he_dust_sputter_temp) then
            ! Total local FUV flux is stellar FUV + background Habing field.
            total_local_FUV = Gflux + he_Gzero*f_ext(rho, ndens, temp)*strat_factor
            ! The heating of dust due to FUV (the energy that DIDN'T eject PEs.
            call pe_heating(total_local_FUV, xHp, ndens, temp, rho, PeHeat, dust_heat)
            ! Now add heating from cosmic rays on the dust.
            call heating(temp, ndens, rho, 0.0d0, strat_factor, dust_cr_heat)
            ! convert from [ergs/g/s] to [ergs/cm^3/s]
            dust_cr_heat = dust_cr_heat * rho
            dust_heat = dust_heat + dust_cr_heat
            ! All the terms in the dust solver have this normalization.
            dust_heat = dust_heat / (ndens*1e-33) 
            tdust = get_dust_temperature(tdust, temp, ndens, 1.0d0/mu_mol, dust_heat)
        end if

        ncount = ncount + 1
        
        if ( (timestep - total_time) < 1d-6 .or. ncount > max_count) exit

    end do

! For stability, this should always be outside the do loop. Let the 
! integration converge, then set the min temp.
    if (temp .lt. he_tradmin) then
        temp = he_tradmin
        ei   = he_tradmin*TtoEI
    end if

#ifdef timing
    call cpu_time(t2)
#endif

#ifdef write_data
    write(*,'(A,ES13.3E3)') "Trap final energy =", ei
    write(*,'(A,ES13.3E3)') "Trap final temp =", temp
    write(*,'(A,I8)') "Trap num steps =", ncount
    write(*,'(A,ES13.3E3)') "Trap total time =", total_time
    write(*,'(A,ES13.3E3,A)') "Trap took =", t2-t1, "secs"

    close(10)
#endif

else if (he_int_method .eq. "RK4") then

emin = 0d0
eplus = 0d0
PeHeat = 0d0
dust_heat = 0d0
dust_cool = 0d0
mol_cool = 0d0

ncount = 0

! RK5(4)

#ifdef write_data
    open(unit=10, file='rk4.dat')
    write(10,'(2ES13.3E3)') total_time, temp
#endif

#ifdef timing
    call cpu_time(t1)
#endif

! These settings seem to work well for everything from 10 K to 10^6 K.
! ei*1d-6, 1d-6, 1d2, timestep*c1

! Ah ha! Finally figured out what causes it to fail sometimes. Basically
! if you start at a very high temp (and energy) the absolute error will
! not be small enough with ei_init*1d-3 if you are cooling to a much 
! lower temp. For instance, if tmp_init = 1d7 K and tmp_final is ~ 4.0 K,
! then the final energy is smaller than the absolute tolerance to calculate
! it. So basically just make the abs error always prop to eps*(kB / mu / mp *1 K)
! to be safe.

! Abs error for calc:

    !call OdeInit(TtoEI*1.0d0, 1d-2, 1d0, timestep*he_subfactor)
    call OdeInit(ei*1d-6, 1d-6, 1d0, timestep*he_subfactor) ! abs err, rel err, min dt, max dt
        
    eiold = ei

    call OdeDriver(dei_dt, eiold, total_time, timestep, 1d3, ei)
        
    temp = ei/TtoEI

        ! If we are tracking both heating TO dust from PE heating and cooling OF dust
        ! from collisions with hydrogen, we can also track the dust temperature
        ! consistently.
        !
        ! OK to stop tracking dust temperature when T > T_sputter because
        ! calculation assumes thermal equilibrium; history is not needed.
        ! -AT, 2019 May 15
        if (useDustCool .and. temp .le. he_dust_sputter_temp) then
            ! Total local FUV flux is stellar FUV + background Habing field.
            total_local_FUV = Gflux + he_Gzero*f_ext(rho, ndens, temp)*strat_factor
            ! The heating of dust due to FUV (the energy that DIDN'T eject PEs.
            call pe_heating(total_local_FUV, xHp, ndens, temp, rho, PeHeat, dust_heat)
            ! Now add heating from cosmic rays on the dust.
            call heating(temp, ndens, rho, 0.0, strat_factor, dust_cr_heat)
            ! convert from [ergs/g/s] to [ergs/cm^3/s]
            dust_cr_heat = dust_cr_heat * rho
            dust_heat = dust_heat + dust_cr_heat
            ! All the terms in the dust solver have this normalization.
            dust_heat = dust_heat / (ndens*1e-33) 
            tdust = get_dust_temperature(tdust, temp, ndens, 1.0d0/mu_mol, dust_heat)
        end if

    if (temp .lt. he_absTmin) then
        temp = he_absTmin
        ei   = he_absTmin*TtoEI
    end if

#ifdef timing
call cpu_time(t2)
#endif

#ifdef write_data
    write(*,'(A,ES13.3E3)') "RK4 final energy =", ei
    write(*,'(A,ES13.3E3)') "RK4 final temp =", temp
    write(*,'(A,ES13.3E3,A)') "RK4 took =", t2-t1, "secs"
#endif

else

    call Driver_abortFlash("[heating_and_cooling]: Unsupported or unrecognized method called. Aborting.")

end if


end subroutine heating_and_cooling


function get_dust_temperature(tdust_in, temp_in, nden, muinv, phen_heat)

!#define DEBUG_DUSTSOLVE

      use Heat_data, only : he_absTmin
      use Cool_data, only : mp

      real, intent(IN) :: tdust_in, temp_in, nden, muinv, phen_heat ! heating from sources - JW
      real :: get_dust_temperature

      real, parameter :: tol = 1.e-4, pw6 = 1./6., smallt = .1

      real, parameter :: pi       = 3.141592653589793238E0
      real, parameter :: kB       = 1.38065812E-16
      real, parameter :: Newton   = 6.6725985E-8

      integer, parameter :: max_newton = 100
      integer :: l

      real :: rho, temp_gas

      real :: tdust, error, chi, lambda_J, kappa, lambda_dust, tau, eqfunc, det
      real :: tnew, coll, heating ! Assume density dependent heating from field stars -JW

! calculate new dust temperature
! use Goldsmith ApJ 557, 2001 values
! use "opacity" 1/Jeans_length in the optical thick regime
! cosmic ray heating: assume shelding factor chi = 1.d-4
! sigma_eff = 6.8533297d-05 (from L_dust = kappa * sigma_eff * tdust**4)
!
!       eqfunc = 3.9e-28*nden - kappa_eff*sigma_eff*tdust**4.0 +  2.e-33 * nden**2 *
!                 sqrt(temp/10.0) * (temp - tdust)
!
! normalize with 1.d-33/nden

#ifdef DEBUG_DUSTSOLVE
  write(*,*) "Calling dust temp solver."
#endif

      error = 10.*tol

      tdust = tdust_in
      
      if (tdust < he_absTmin) tdust = he_absTmin
      if (temp_in < he_absTmin) then
        temp_gas = he_absTmin
      else
        temp_gas = temp_in
      end if
      
      ! Should we force dust to be the gas temp? Might try switching this
      ! off. - JW 12/15/16
      !if(tdust .gt. temp_gas) tdust = temp_gas
      get_dust_temperature = tdust

      rho = nden*mp/muinv
      !lambda_J    = sqrt(pi*kB*temp_gas*muinv/(Newton*rho))
      
	  !heating = 3.9e5 + phen_heat
      heating =  phen_heat
         
! iterate until convergence
      l = 1

      do while ( (l <= max_newton) .AND. (error > tol) )

! Here we comment this out, following Draine ApJ 732:100 (12pp) 2011
! stating that the dust sublimation is not important in HII regions that
! we simulate, unlike the protoplanetary disks that Semenov 2008 did
! that led to this (note I don't expect we will ever reach densities
! above 10^6, maybe in extreme cases 10^8). - JW 

! Further, because we never get as dense as Banerjee's original work,
! we will just assume the dust is always optically thin to itself. - JW

!! chi is likely opacity per unit volume.

!         chi         = 3.3d-26 * (tdust*0.054833588)**2 * nden
!         if (tdust .GT. 200.) then
!            if (tdust .LT. 1500.) then
!! temperature independent opacity (see Semenov et al., A&A 410, 2003)
!               chi = 3.9688735d-24 * nden
!            else
!! dust melting kappa ~ T^-12
!               chi = 3.9688735d-24 * nden * (tdust*0.00066666667)**(-12.0)
!            endif
!         endif
!         kappa       = min(chi,1.0/lambda_J)*6.8533297e+28/nden
!         lambda_dust = kappa*tdust**4
!         tau         = chi * lambda_J

!         eqfunc  =  3.9E5 +  0.63245552*nden*(temp - tdust)*sqrt(temp) &
!              - lambda_dust
         
! Dust cooling in [ergs cm^-3 s^-1] from Goldsmith 2001 eqn 13. - JW

! Note all eqns are normed by dividing by ndens*1e-33. - JW

         lambda_dust = 6.8 * tdust**6.0
         
! Collisional term with gas from Goldsmith 2001. - JW

! Note this will actually cause the gas to cool too quickly in wind
! regions, so coll cooling of the gas should be switched off in gas 
! above ~ 10^5 K. Switch off before metal line cooling.
!         if (temp_gas .le. 1d5) then

! Goldsmith 2001.
           !coll = 0.63245552*nden*(temp_gas - tdust)*sqrt(temp_gas)
! Hollenbach and McKee 1989, as reported in Glover and Clark 2011.
           coll = 3.8d0 * nden *(1.0d0 - 0.8d0*exp(-75.0d0 / temp_gas)) *(temp_gas - tdust) * sqrt(temp_gas)
!         end if

         eqfunc = heating + coll - lambda_dust
              
!         eqfunc  =  heating +  0.63245552*nden*(temp - tdust)*sqrt(temp) &
!              - lambda_dust
                     
!         if (tau > 1.e0) then 
!            det = -4.0 * lambda_dust/tdust - 0.63245552*nden*sqrt(temp)
!         else
            !det = -6.0 * lambda_dust/tdust - 0.63245552*nden*sqrt(temp_gas)

! With the Hollenbach and McKee collisional bit:
             det = -6.0 * lambda_dust/tdust - 3.8d0*nden*(1.0d0-0.8d0*exp(-75.0d0/temp_gas))*sqrt(temp_gas)
            
!            if (tdust .GT. 200.) then
!!              if (tdust .LT. 1500.) then
!              det = -4.0 * lambda_dust/tdust - 0.63245552*nden*sqrt(temp)
!!              else
!!                det =  8.0 * lambda_dust/tdust - 0.63245552*nden*sqrt(temp)
!!              endif
!            endif
!         endif

         tnew = tdust - eqfunc/det
                  
         error = abs(tnew - tdust)/tdust
         
#ifdef DEBUG_DUSTSOLVE
         write(*,'(A,X,ES13.3)') 'tgas=', temp_gas
         write(*,'(A,X,ES13.3)') 'told=', tdust
         write(*,'(A,X,ES13.3)') 'tnew=', tnew
         write(*,'(A,X,3ES13.3)') 'coll, heating, lambda_d', coll, heating, lambda_dust
         call flush(6)
#endif
         
         tdust = tnew

         if (tdust < smallt) then
            tdust = smallt
            error = 0.1*tol
         endif

         l = l + 1

       enddo

#ifdef DEBUG_DUSTSOLVE
  write(*,'(A,I4)') "In dust solver, num of iteratons = ", l
#endif

! were there too many iterations??
       if (l > max_newton) then
          print *, 'Error converging in cool_block'
          print *, 'temp = ', temp_gas !, '(',thisBlock, i,j,k,')'
          print *, 'tdust = ', tdust
          call flush(6)
#ifdef DEBUG_DUSTSOLVE
          stop
#endif
!          call stamp_logfile ("cool_block:  iteration failed to converge")
!          call abort_flash   ("cool_block:  iteration failed to converge")
       endif
       
       if (tdust < he_absTmin) tdust = he_absTmin

       get_dust_temperature = tdust

       return

end function get_dust_temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Heating
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine heating(temp_in, ndens_in, rho_in, ephen_in, strat_factor_in, eplus_out)
    
use Heat_data ! Use Christian's heating data file.
#ifdef FERVENT
use rt_data, only : dust_gas_ratio
#else
real, parameter :: dust_gas_ratio = 0.01
#endif

integer, parameter      :: dp = KIND(1d0)

real(dp), intent(in)    :: temp_in, ndens_in, rho_in, ephen_in, strat_factor_in
real(dp), intent(inout) :: eplus_out
!                          lower_limit, upper_limit, 10% of ll, 10% of ul, the returned linear_interp factor
real(dp)                :: ll, ul, ll10, ul10, lin_factor


    eplus_out = 0.0d0
    
    !if (ephen_in > 0.0d0) print*, "[heating]: ephen = ", ephen_in

! add UV heating
    if ( (temp_in <= he_theatmax) .and. (temp_in >= he_theatmin) ) then

! photoelectric heating
! This is the BACKGROUND PE heating due to background starlight. - JW
! Note this should be less effective due to extinction in clouds,
! which still needs to be checked! - JW

! Here I add a local approximation for self-shielding based on the Jeans length
! Note a 2nd order correction to this would be the local density length scale:
! n / (dn/dr) = lambda_n, for which we could take N_H = n*min(lambda_J, lambda_n).

! See Bohlin+ 1978 for details and SILCC paper 1 for an example implementation of
! extinction without the approximation. See Banerjee, Pudritz and Anderson 2006 for
! example of this type of extinction approximation.
 
        if (he_use_cr_heating) then
          if (he_crIonNH == 1.0 .and. he_crIonExp == 1.0) then ! uniform
            eplus_out = strat_factor_in*he_crIonRate*he_crIonEnergy / (rho_in/ndens_in)
          else ! Use eqn 27 from Padovani et. al. 2009 
            eplus_out = strat_factor_in*( he_crIonRate*he_crIonEnergy &
                *(max(he_crIonNH,approx_column_dens(rho_in, ndens_in, temp_in))/he_crIonNH)**(-he_crIonExp)) &
                / (rho_in/ndens_in)
          end if
        end if

        ! Make a gentle linear transistion from heating on to heating off
        ! from 90% to 100% of temp at which user specified to switch off heating.
        ll   = he_theatmin
        ll10 = ll+0.1*ll
        ul   = he_theatmax
        ul10 = ul-0.1*ul
        
        if (temp_in >= ul10) then
            lin_factor = lin_interp(ul10, 1.0, ul, 0.0, temp_in)
        else if (temp_in <= ll10) then
            lin_factor = lin_interp(ll, 0.0, ll10, 1.0, temp_in)
        else
            lin_factor = 1.0d0
        end if
        eplus_out = eplus_out*lin_factor
    else
        eplus_out = 0.0
    endif
#ifdef ONE_CELL_TESTING
    write(*,'(A,ES13.3)') "CR heating=", eplus_out*rho_in
    write(*,'(A,ES13.3)') "CR lin_factor=", lin_factor
    write(*,'(A,ES13.3)') "CR temp=", temp_in
    write(*,'(A,ES13.3)') "CR he_theatmax=", he_theatmax
    write(*,'(A,ES13.3)') "CR he_theatmin=", he_theatmin
#endif

return
end subroutine heating

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Photoelectric heating - JW (yay I added new things...)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pe_heating(Gflux_in, xHp_in, ndens_in, temp_in, rho_in, PeHeat, dust_heating)
    ! returns PeHeat [erg g^-1 s^-1], gas heating rate from PE effect,
    ! and dust_heating [erg cm^-3 s^-1], dust heating rate from FUV photons
    ! that DIDN'T eject photo-electrons.

    use cool_vars
    use Heat_data

    real(dp), intent(in)  :: Gflux_in, xHp_in, ndens_in, temp_in, rho_in
    real(dp), intent(out) :: PeHeat, dust_heating
    ! PE heating flux in normalized Habing units, efficiency
    real              :: Gfactor, eff, xHp_local
    real(dp)          :: ll, ul, ll10, ul10, lin_factor

    ! W&D heating constants (Rv=3.1, bc=6x10^-5, case A, B0 field):
    real, parameter   :: C0=7.64d0, C1=4.52d0, C2=0.04371d0, C3=0.00557d0, &
                         C4=0.132d0, C5=0.452d0, C6=0.675d0

    PeHeat       = 0.0d0
    dust_heating = 0.0d0

    ! Limit temp_in to values below 10^5 K (reference: Simon Glover).
    ! This is because you can have PE heating in very hot gas where the
    ! dust should be destroyed, i.e. in stellar winds. - JW

    ! Note based on the plot in Draine's book (fig. 25.4) I think
    ! this number should be between 10^5.5 and 10^6 K since it appears
    ! dust sputters at a rate of ~ 1 cm^3 / 100 yrs at 10^6 K.
    ! Further note that this assumes the code can properly find the
    ! equilibrium temperature between atomic cooling and dust heating at
    ! around 2x10^4 K as in Bakes and Tielens 1994 fig 16. If not, then
    ! this should be set manually to turn off heating around that
    ! temperature, but cooling may still need to be effective up to
    ! 10^6 K. 3/26/18 -JW
    if ( (temp_in <= min(he_theatmax, he_dust_sputter_temp)) .and. &
         (temp_in >= he_theatmin) ) then

        ! Use current ionization as approx electron density
        ! since we only follow hydrogen in the code. - JW

        ! No NaNs please!
        if (xHp_in .eq. 0.0d0) then
            xHp_local = 1d-8
        else
            xHp_local = xHp_in
        end if

#ifdef ONE_CELL_TESTING
        write(*,'(A,X,4ES13.3)') "Gflux_in, temp_in, xHp_local, ndens_in", &
                                  Gflux_in, temp_in, xHp_local, ndens_in
#endif

        if (he_pe_recipe .eq. "BT94") then
        ! Wolfire+ 2003 PE heating, following Bakes/Tielens 1994
        ! We assume phi_PAH = 0.5 as in Wolfire.

            Gfactor = Gflux_in*sqrt(temp_in) / (xHp_local*ndens_in*0.5d0)

            eff = 4.9d-2 / (1.0 + 4.0d-3 * Gfactor**0.73) &
                + (3.7d-2 * (temp_in / 1.0d4)**0.7) / (1.0 + 2.0d-4 * Gfactor)

            ! eqn 19, specific energy units [ergs g^-1 s^-1].
            PeHeat = he_pe_norm*eff*Gflux_in / (rho_in/ndens_in)

            ! Heating of the dust itself from Goldsmith 2001 eqn 4 in [ergs cm^-3 s^-1].
                        !dust_heating = 3.9e-24 * (ndens / 2.0) * Gflux
            ! Convert to [ergs g^-1 s^-1].
                        !dust_heating = dust_heating / (rho/ndens)
            ! or if you like, its (1-eff) if eff is the efficiency of electrons to
            ! be liberated from the grains. This seems correct since Goldsmith and
            ! Wolfire are so close on their coefficients here if we note that the
            ! efficiency should be about 10% as stated in Stahler and Palla pg 54.

            ! Note that this means that no photons above 13.6 eV inside the HII
            ! region are being captured here, even though some of the >13.6 eV are
            ! hitting dust when they are very close to the star.

            ! Volumetric units [ergs cm^-3 s^-1].
            dust_heating = he_pe_norm*max(0.0,(1.0-eff))*Gflux_in * ndens_in

        else if (he_pe_recipe .eq. "WD01") then
        ! Weingarter and Draine 2001 PE Heating using R_v=3.1, b=6d-5 and B0
        ! field from table 2.

            Gfactor = Gflux_in*sqrt(temp_in) / (xHp_local*ndens_in)

            eff    = (C0 + C1*temp_in**C4) / &
                     (1.0d0 + (C2*Gfactor**C5)*(1.0d0 + C3*Gfactor**C6))

            ! Specific energy units [ergs g^-1 s^-1].
            PeHeat = 1d-26*eff*Gflux_in / (rho_in/ndens_in)

            ! Volumetric units [ergs cm^-3 s^-1].
            dust_heating = (he_pe_norm - 1d-26*eff)*Gflux_in * ndens_in
            dust_heating = max(dust_heating, 0.0)

        else if (he_pe_recipe .eq. "NOGF") then
        ! No G factor stuff, set FUV to a simple knob following Joung & M-MML
        ! 2006 and Hill+ (incl M-MML, JCIM) 2018. -AT, 2019 June 19

            eff = 0.05

            ! Specific energy units [ergs g^-1 s^-1].
            PeHeat = he_pe_norm*eff*Gflux_in / (rho_in/ndens_in)

            ! Volumetric units [ergs cm^-3 s^-1].
            dust_heating = he_pe_norm*(1.0-eff)*Gflux_in * ndens_in

        else

            call Driver_abortFlash("[heatCool] Invalid he_pe_recipe")

        end if

        ! Make a gentle linear transistion from heating on to heating off
        ! from 90% to 100% of temp at which user specified to switch off heating.
        ll   = he_theatmin
        ll10 = ll+0.1*ll
        ul   = min(he_theatmax, he_dust_sputter_temp)
        ul10 = ul-0.1*ul

        if (temp_in >= ul10) then
            lin_factor = lin_interp(ul10, 1.0, ul, 0.0, temp_in)
        else if (temp_in <= ll10) then
            lin_factor = lin_interp(ll, 0.0, ll10, 1.0, temp_in)
        else
            lin_factor = 1.0d0
        end if

        PeHeat = PeHeat*lin_factor

    end if  ! (temp_in <= ... .and. temp_in >= ...)

#ifdef ONE_CELL_TESTING
    write(*,'(A,ES13.3)') "PE heating=", PeHeat*rho_in
    write(*,'(A,ES13.3)') "PE eff=", eff
    write(*,'(A,ES13.3)') "PE on dust heating=", dust_heating
#endif

return
end subroutine pe_heating

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Atomic cooling from Dalgarno and McCray - JW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cooling(temp_in, ndens_in, conf_in, emin_out)

use Heat_data ! Use Christian's heating data file.
!use Cool_data ! Use Robi's cooling data file.
use cool_vars, only : xHp 

integer, parameter  :: dp = KIND(1d0)

real(dp), intent(in)  :: temp_in, ndens_in, conf_in
real(dp), intent(inout) :: emin_out
real(dp)          :: ll, ul, ll10, ul10, lin_factor
logical :: DGionFix = .true. !.false.
real(dp) :: DGfactor


    emin_out = 0.0d0
    ! checks for density and temperature range allowed to radiatively cool
    ! heatmin/max should be consistent with tradmin/max
    ! no heating if outside range
    if (  (temp_in   <= he_tradmax .AND. temp_in   >= he_tradmin) & 
     .AND. (ndens_in >= he_dradmin .AND. ndens_in <= he_dradmax) ) then
        
        ! get da cooling rate [erg cm^3/s] (volumetric)
        call Radloss(temp_in,emin_out)

        if (DGionFix) then

    ! NOTE: D&G curve this is based on has lower cooling for
    !       temps < 1e4 K and ionization fractions smaller than 0.01.
    !       Here I'm going to "fix" this curve, which only has the
    !       curve for xHp = 0.01 ionization, by dividing out the
    !       difference in the curves (see fig 2 in Dalgarno and McCray 72).
    !       I'll at least try and be conservative here. - JW
    
          if (xHp > 1e-2 .and. temp_in < 1.2e4) then
              ! Roughly 2-10x diff b/t xHp=1e-1 and xHp=1e-2 curve.
              DGfactor = lin_interp(1e-2, 1.0, 1e-1, 5.0, min(1e-1,xHp))
              emin_out = emin_out * DGfactor
              ! piecewise-continuous, monotonic transition to CIE curve
              if (temp_in > 1e4) then
                  emin_out = emin_out / lin_interp(1e4, 1.0, 1.2e4, DGfactor**0.5, temp_in)**2
              end if
          else if (xHp < 1e-2 .and. temp_in < 1e4) then
              ! Roughly 2-5x diff b/t xHp=1e-2 and xHp=1e-4 curve.
              DGfactor = lin_interp(1e-4, 0.35, 1e-2, 1.0, max(1e-4,xHp))
              emin_out = emin_out * DGfactor
              ! piecewise-continuous, monotonic transition to CIE curve
              if (temp_in > 9e3) then
                  emin_out = emin_out / lin_interp(9e3, 1.0, 1e4, DGfactor, temp_in)
              end if
          end if

        end if
        
        ! Make a gentle linear transistion from cooling on to cooling off
        ! from 90% to 100% of temp at which user specified to switch off cooling.
        ll   = he_tradmin
        ll10 = ll+0.1*ll
        ul   = he_tradmax
        ul10 = ul-0.1*ul
        
        if (temp_in >= ul10) then
            lin_factor = lin_interp(ul10, 1.0, ul, 0.0, temp_in)
        else if (temp_in <= ll10) then
            lin_factor = lin_interp(ll, 0.0, ll10, 1.0, temp_in)
        else
            lin_factor = 1.0d0
        end if
        
        emin_out = emin_out*lin_factor
        
#ifdef ONE_CELL_TESTING
    write(*,'(A,ES13.3)') "Atomic cooling=", emin_out*ndens_in*ndens_in
#endif

        ! convert to erg/(g s)
        emin_out = emin_out*conf_in

    else
        emin_out = 0.0d0
       
    endif

return
end subroutine cooling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Cooling method from Robi Banerjee and Daniel Seifried - JW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Molecular cooling (from Neufield et al 1995)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine molecular_cooling(ndens_in, temp_in, mu_mol_in, cool_pw)

use Cool_data
use Heat_data, only : he_theatmax, he_tradmin

integer, parameter  :: dp = KIND(1d0)

real(dp), intent(in)  :: ndens_in, temp_in, mu_mol_in
real(dp), intent(out) :: cool_pw
real(dp)          :: ll, ul, ll10, ul10, lin_factor

integer           :: ixtemp, ixdens 

! Molecular cooling from Neufield et al.
! Note that Banerjee et al 2006 states in appendix that dust cooling is only needed
! for densities greater than 10^7.5. - JW

    cool_pw = 0.0

! Switch off if Temp is too high (more than 2x10^4, so use theatmax instead of tradmax) (Simon says!).
! Could be improved by adding div (vel) for the Neufeld data (Simon says!).

! The min density should be n_H ~ 1000 cm^-3 in the flash par file (Simon says!).

    if ((ndens_in.ge.nd_cool_min).and.(ndens_in.le.nd_cool_max) .and. &
        (temp_in.ge.T_cool_min) .and. (temp_in .le. he_theatmax)) then  !temp_in.ge.he_tradmin
#ifdef timing                  
!                  t1 = MPI_WTime()
#endif                  
        ixtemp = find(cool_dat(1,:,1),max(temp_in,T_min),TEMP_PTS)
        ixdens = find(cool_dat(:,1,2),max(ndens_in,nd_min),DENS_PTS)

#ifdef timing
!                  t2 = MPI_WTime()
!                  mol_time = mol_time + t2 - t1
#endif


        if ( ixtemp .gt. 0 .and. ixtemp .le. TEMP_PTS .and. &
            ixdens .gt. 0 .and. ixdens .le. DENS_PTS ) then

            ! cooling power per H2 molecule
            cool_pw = cool_dat(ixdens,ixtemp,3) / (mu_mol_in * mp)

        else
             cool_pw = 0.0
             write(*,*) 'No cooling value found, T = ', temp_in   &
                  ,' ndens = ',ndens_in
        endif

        ! Make a gentle linear transistion from cooling on to cooling off
        ! from 90% to 100% of temp at which user specified to switch off cooling.
        ll   = T_cool_min
        ll10 = ll+0.1*ll
        ul   = he_theatmax
        ul10 = ul-0.1*ul
        
        if (temp_in >= ul10) then
            lin_factor = lin_interp(ul10, 1.0, ul, 0.0, temp_in)
        else if (temp_in <= ll10) then
            lin_factor = lin_interp(ll, 0.0, ll10, 1.0, temp_in)
        else
            lin_factor = 1.0d0
        end if

        cool_pw = cool_pw*lin_factor
        
        ! Do the same except now for the gas density.
        ll   = nd_cool_min
        ll10 = ll+0.1*ll
        ul   = nd_cool_max
        ul10 = ul-0.1*ul
        
        if (ndens_in >= ul10) then
            lin_factor = lin_interp(ul10, 1.0, ul, 0.0, ndens_in)
        else if (ndens_in <= ll10) then
            lin_factor = lin_interp(ll, 0.0, ll10, 1.0, ndens_in)
        else
            lin_factor = 1.0d0
        end if

        cool_pw = cool_pw*lin_factor

    end if
! advected internal energy
!                   ei_adv = solnData(ENER_VAR,i,j,k) - ek


#ifdef ONE_CELL_TESTING
  write(*,'(A,ES13.3)') "Molecular cooling", cool_pw*ndens_in*mu_mol_in*mp
#endif

end subroutine molecular_cooling

subroutine dust_cooling(ndens_in, temp_in, tdust_in, rho_in, cool_dg)

use Heat_data, only : he_dust_sputter_temp

integer, parameter  :: dp = KIND(1d0)

real, intent(in)  :: ndens_in, temp_in, tdust_in, rho_in
real, intent(out) :: cool_dg
real              :: lambda_gd
real              :: ll, ul, ll10, ul10, lin_factor


cool_dg = 0.0

! cooling power by dust (Goldsmith ApJ 557, Eq.(15)
! cooling of the gas by dust in erg cm^-3 s^-1 - JW

! Note this will actually cause the gas to cool too quickly in wind
! regions, so coll cooling of the gas should be switched off in gas 
! above ~ 10^5 K. Switch off before metal line cooling!!!

  if (temp_in .le. he_dust_sputter_temp) then

	! This Goldsmith coeff is too small, Hollebach and McKee is better. See
	! Glover and Klessen 2014 review for more.
	! Goldsmith 2001
	!lambda_gd   = 2.e-33 * ndens_in**2 *  (temp_in - tdust_in) * sqrt(temp_in*0.1)
	! Hollenbach and McKee 1989, as reported in Glover and Clark 2011
	lambda_gd  = 3.8d-33 * ndens_in**2.0 *(1.0d0 - 0.8d0*exp(-75.0d0 / temp_in)) *(temp_in - tdust_in) * sqrt(temp_in)
	cool_dg = lambda_gd / rho_in  ! convert to ergs g^-1 s^-1
	
	! Make a gentle linear transistion from cooling on to cooling off
	! from 90% to 100% of temp at which user specified to switch off cooling.
	ul   = he_dust_sputter_temp
	ul10 = ul-0.1*ul
	if (temp_in >= ul10) then
	    lin_factor = lin_interp(ul10, 1.0, ul, 0.0, temp_in)
	else
	    lin_factor = 1.0d0
	end if

	cool_dg = cool_dg*lin_factor

#ifdef ONE_CELL_TESTING
    write(*,'(A,ES13.3)') "dust cooling=", lambda_gd
#endif
  end if

! Here we comment this out, following Draine ApJ 732:100 (12pp) 2011
! stating that the dust sublimation is not important in HII regions that
! we simulate, unlike the protoplanetary disks that Semenov 2008 did
! that led to this (note I don't expect we will ever reach densities
! above 10^6, maybe in extreme cases 10^8). - JW 

!! optical depth Goldsmith 2001 Eq.10 (tau = kappa * rho * lambda_J
!! 3.8d11 Hz = 18.236997 K; 1/18.236997 = 0.054833588
!                  if( rho .ge. 1.e-12) then
!                  ! Here they use the Jean's length as a proxy for radius of the cloud. - JW
!                    lambda_J = sqrt(pi*kB*temp/(mu_mol*mp)/(Newton*rho))
!                    chi      = 3.3d-26 * (tdust*0.054833588)**2 * ndens
!                    if (tdust .GT. 200.) then
!                      if (tdust .LT. 1500.) then
!! temperature independent opacity (see Semenov et al., A&A 410, 2003)
!                        chi = 3.9688735d-24 * ndens
!                      else
!! dust melting kappa ~ T^-12
!                        chi = 3.9688735d-24 * ndens * (tdust*0.00066666667)**(-12.0)
!                      endif
!                    endif

!                    tau      = chi * lambda_J
!                  else
!                    tau = 0.0
!                  endif

!! Lets try switching the dust cooling solver back on... - JW 11/3/16
!                  if (tdust_solver .and. dust_cool) &
!                    tdust = get_dust_temperature(tdust, temp, ndens, 1./(mu_mol*mp), dust_heating)
                  
!                  if (tdust <= 0.0) then
!                     print*, "tdust = 0 after second solver call." !, he_meshMe, blockID, i, j, k
!                     tdust = 5.0
!                  end if

return
end subroutine dust_cooling


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ODE for rate of change in internal energy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dei_dt(t_in, ei_in)

!#define DEBUG_RADHEAT1
!#define DEBUG_RADHEAT2

use cool_vars
use Cool_data, only : useDustCool, T_cool_min
use Heat_data, only : he_Gzero, he_absTmin, he_tradmin, he_theatmax

!integer, parameter  :: dp = KIND(1d0)

real(dp), intent(in) :: t_in, ei_in
real(dp)             :: dei_dt, tmp_temp
real(dp)             :: emin, cool_pw, cool_dg        ! Cooling things.
real(dp)             :: eplus, PeHeat, dust_heating   ! Heating things.
real(dp)             :: cooling_floor
real(dp)             :: total_local_FUV

integer          :: nprint

          
          dei_dt   = 0.0d0
          tmp_temp = 0.0d0
          emin=0.0d0; cool_pw=0.0d0; cool_dg=0.0d0
          eplus=0.0d0; PeHeat=0.0d0; dust_heating=0.0d0
          cooling_floor = 0.0d0

          tmp_temp = ei_in / TtoEI

          ! Enforce absolute minimums on temperatures. Care that energy and temp are still related correctly.

          !if (temp < 10.0) temp  = 10.0
          !if (tdust < 5.0) tdust = 5.0
#ifdef DEBUG_RADHEAT2
          !  nprint = 0
          !if (nprint < 100) then
          write(*,*)              "On entry."
          write(*,'(A,ES13.3E3)') "ei_in =", ei_in
          write(*,'(A,ES13.3E3)') "rho =", rho
          write(*,'(A,ES13.3E3)') "ndens =", ndens
          write(*,'(A,ES13.3E3)') "mu_mol =", mu_mol
          write(*,'(A,ES13.3E3)') "eplus =", eplus
          write(*,'(A,ES13.3E3)') "emin =", emin
          write(*,'(A,ES13.3E3)') "PeHeat =", PeHeat
          write(*,'(A,ES13.3E3)') "temp =", tmp_temp
          write(*,'(A,ES13.3E3)') "dei_dt =", dei_dt
          write(*,'(A,ES13.3E3)') "t =",  t_in
          call flush(6)
#endif
          !nprint = nprint + 1
          if (tmp_temp < 0.0d0) then
            tmp_temp = 1.0
          end if


! Get all the heating from different stuff (radiation, SNe, Photoelectric on gas, rad on dust...)
          !call heating(tmp_temp, Gflux, ndens, rho, ephen, eplus, PeHeat, dust_heating)
          ! I would like to eventually move away from calling heating at all for
          ! FUV (only using it for CR) and instead
          ! just add in a base flux to all cells here from background UV
          ! ^^^^ This is now done! - JW
          
#ifdef ONE_CELL_TESTING
    write(*,'(A,ES13.3)') "stellar heating=", ephen*rho
    write(*,'(A,ES13.3)') "FUV background heating=", he_Gzero*f_ext(rho, ndens, tmp_temp)*strat_factor
#endif
          call heating(tmp_temp, ndens, rho, 0d0, strat_factor, eplus)

          ! Add in background FUV field from stars (note Gflux and Gzero
          ! are already in Habing field units.)
          total_local_FUV = Gflux + he_Gzero*f_ext(rho, ndens, tmp_temp)*strat_factor
          !call pe_heating(Gflux, xHp, ndens, tmp_temp, rho, PeHeat, dust_heating)
          call pe_heating(total_local_FUV, xHp, ndens, tmp_temp, rho, PeHeat, dust_heating)

          !if (tmp_temp < he_theatmax) &
             eplus = eplus + ephen

#ifdef DEBUG_RADHEAT2
          !  nprint = 0

          !if (nprint < 100) then
          !if (PeHeat > 0.0d0 .and. ephen > 0.0d0 .or. t_in .lt. 1.0) then
          write(*,*)              "After heating."
          write(*,'(A,ES13.3E3)') "ei_in =", ei_in
          write(*,'(A,ES13.3E3)') "eplus =", eplus
          write(*,'(A,ES13.3E3)') "emin =", emin
          write(*,'(A,ES13.3E3)') "PeHeat =", PeHeat
          write(*,'(A,ES13.3E3)') "temp =", tmp_temp
          write(*,'(A,ES13.3E3)') "dei_dt =", dei_dt
          write(*,'(A,ES13.3E3)') "t =",     t_in
          call flush(6)
#endif
          !if (t_in .lt. 1.0) stop
		  if (tmp_temp < 0.0d0) then
            write(*,*) "WARNING! Temperature is now negative after applying heating!"
          end if
          !nprint = nprint + 1
          
          !end if
#ifdef ONE_CELL_TESTING
    write(*,'(A)') "Getting cooling floor."
#endif
          
! Calculate the cooling at the temperature floor and subtract it
! away to smoothly approach the temperature floor.
! Okay, this cooling min must reflect the proper individual cooling
! minimums for each function, i.e. T_cool_min for molecular, he_tradmin
! for atomic cooling, and he_absTmin for dust.
! Indeed, setting the floor this way in theory should allow for us to
! totally remove the min temp checks on the individual functions, but
! we leave them in for safety. Note also, the mins in the functions MUST
! be >= checks that INCLUDE these floors so we get something back here!
          call cooling(he_tradmin, ndens, conf, emin)
          if (useDustCool) then
            call molecular_cooling(ndens, T_cool_min, mu_mol, cool_pw)
            ! Don't include dust cooling in the floor, since its based on the
            ! difference b/t the dust temp and the gas temp and the dust
            ! temp is already floored at the gas temp, so this would be
            ! zero there.
            !call dust_cooling(ndens, T_cool_min, tdust, rho, cool_dg)
          end if
            
          cooling_floor = emin + cool_pw ! + cool_dg

#ifdef ONE_CELL_TESTING
    write(*,'(A)') "Getting real cooling."
#endif
! Get all the cooling from stuff (atomic cooling, dust cooling, molecular cooling)
          call cooling(tmp_temp, ndens, conf, emin)
          if (useDustCool) then
            call molecular_cooling(ndens, tmp_temp, mu_mol, cool_pw)
            call dust_cooling(ndens, tmp_temp, tdust, rho, cool_dg)

! This should be done at the end of the loop... don't change any variables during the root finding process.
!            ! be warned, dust heating is in ergs cm^-3 s^-1
!            dust_heating = dust_heating / ndens / 1e-33 ! convert for dust solver
!            tdust = get_dust_temperature(tdust, tmp_temp, ndens, 1.0d0/mu_mol, dust_heating)
          end if
          !call cooling(tmp_temp, ndens, rho, conf, tdust, dust_heating, emin)
          
          dei_dt = eplus + PeHeat - emin - cool_pw - cool_dg + cooling_floor
          
! Note that if we would cool below 10 K, set delta e = 0 so that we don't
! do that (and so the RK4 converges rapidly). Note also this requires
! that Radloss actually report what the cooling should be at and below
! 10 K otherwise this won't work!

! Also, this should be moved outside the integrator so that it can converge
! normally to some temperature, then we can apply a floor after integration.
!          if (tmp_temp .le. 10.0_dp .and. dei_dt .lt. 0.0_dp) then
!            !print*, "Yo!"
!            dei_dt = 0.0_dp
!          end if
          
#ifdef DEBUG_RADHEAT2
          !  nprint = 0

          !if (PeHeat > 0.0d0 .and. ephen > 0.0d0 .or. t_in .lt. 1.0) then
          write(*,*)              "After cooling."
          write(*,'(A,ES13.3E3)') "ei_in =", ei_in
          write(*,'(A,ES13.3E3)') "eplus =", eplus
          write(*,'(A,ES13.3E3)') "ephen =", ephen
          write(*,'(A,ES13.3E3)') "PeHeat =", PeHeat
		  write(*,'(A,ES13.3E3)') "atomic cooling =", emin
          write(*,'(A,ES13.3E3)') "mol cooling =", cool_pw
          write(*,'(A,ES13.3E3)') "dust cooling =", cool_dg
		  write(*,'(A,ES13.3E3)') "cooling floor=", cooling_floor
          write(*,'(A,ES13.3E3)') "Heating TO dust =", dust_heating
          write(*,'(A,ES13.3E3)') "temp =", tmp_temp
          write(*,'(A,ES13.3E3)') "dust temp =", tdust
          write(*,'(A,ES13.3E3)') "dei_dt =", dei_dt
          write(*,'(A,ES13.3E3)') "t =",     t_in
          call flush(6)
          
          !nprint = nprint + 1
#endif
          
		  if (tmp_temp < 0.0d0) then
            write(*,*) "WARNING! Temperature is now negative after applying cooling!"
          end if

#ifdef debug
print*, "emin after mol/dust cooling =", emin
print*, "cool_pw =", cool_pw
print*, "cool_dg =", cool_dg
call flush(6)
#endif

          return
          

end function dei_dt

subroutine Radloss(T,radia)
!
    implicit none
    
    integer, parameter  :: dp = KIND(1d0)

    real(dp), intent(IN) :: T
    real(dp), intent(OUT) :: radia

!
!  The factor of 5.647e-101 comes from changing units from cgs to kpc, Myr, & Msun
!
if (T .ge. 1.70e4) then
  if (T .ge. 5.62e5) then
    if (T .ge. 2.75e6) then
      if (T .ge. 3.16e7) then
    radia = 3.090e-27*sqrt(T) !*5.647d-101
      else 
    radia = 5.188e-21/T**0.33 !*5.647d-101
      endif
    else 
      if (T .ge. 1.78e6) then
    radia = 3.890e-4/T**2.95 !*5.647d-101
      else 
    radia = 1.3e-22*(T/1.5e5)**0.01 !*5.647d-101
      endif
    endif
  else 
    if (T .ge. 7.94e4) then
      if (T .ge. 2.51e5) then
    radia = 3.98e-11/T**2 !*5.647d-101
      else 
    radia = 6.31e-22*(T/1.e6)**0.01 !*5.647d-101
      endif
    else 
      if (T .ge. 3.98e4) then
    radia = 1.e-31*T**2 !*5.647d-101
      else 
    radia = 1.479e-21/(T**0.216d0) !*5.647d-101
      endif
    endif
  endif
else
  if (T .ge. 1.e3) then
    if (T .ge. 6.31e3) then
      if (T .ge. 1.0e4) then
    radia = 7.63d-81*(T**13.8d0) !*5.647d-101
      else 
    radia = 3.13e-27*(T**0.396d0) !*5.647d-101
      endif
    else 
      if (T .ge. 3.16e3) then
    radia = 2.64e-26*(T**0.152d0) !*5.647d-101
      else 
    radia = 5.28e-27*(T**0.352d0) !*5.647d-101
      endif
    endif
  else 
    if (T .ge. 3.98e1) then
      if (T .ge. 2.00e2) then
    radia = 3.06e-27*(T**0.431d0) !*5.647d-101
      else 
    radia = 1.52e-28*(T**0.997d0) !*5.647d-101
      endif
    else 
      if (T .ge. 2.51e1) then
    radia = 2.39e-29*(T**1.50d0) !*5.647d-101
      else 
    radia = 1.095e-32*(T**3.885d0) !*5.647d-101
      endif
    endif
  endif
endif
! don't cool below 10 K
!if (T .le. 1.e1) radia = 0.e-28 !*5.647d-101  !below T=0, no cooling

    return
end subroutine Radloss

function approx_column_dens(rho_in, ndens_in, temp_in)

use Heat_data, only : he_boltz
#ifdef FERVENT
use rt_data, only : dust_gas_ratio
implicit none
#else
implicit none
real, parameter :: dust_gas_ratio = 0.01
#endif

integer, parameter      :: dp=KIND(1d0)
real(dp), intent(in)    :: temp_in, ndens_in, rho_in
real(dp)                :: lam_j, approx_column_dens
real(dp), parameter     :: Newton = 6.6725985d-8, pi = 3.1415926535897932384d0

! Here I add a local approximation for column density.

    lam_j              = 1.0/(rho_in/ndens_in) * sqrt(5./3.*pi*he_boltz*temp_in/ndens_in / Newton)
    approx_column_dens = lam_j * ndens_in

return
end function

function f_ext(rho_in, ndens_in, temp_in)

use Heat_data, only : he_boltz
#ifdef FERVENT
use rt_data, only : dust_gas_ratio
implicit none
#else
implicit none
real, parameter :: dust_gas_ratio = 0.01
#endif

integer, parameter      :: dp=KIND(1d0)
real(dp), intent(in)    :: temp_in, ndens_in, rho_in
real(dp)                :: lam_j, Av, f_ext
real(dp), parameter     :: Newton = 6.6725985d-8, pi = 3.1415926535897932384d0, N_H0=1.87e21

! Here I add a local approximation for self-shielding based on the Jeans length
! Note a 2nd order correction to this would be the local density length scale:
! n / (dn/dr) = lambda_n, for which we could take N_H = n*min(lambda_J, lambda_n).

! See Bohlin+ 1978 for details and SILCC paper 1 for an example implementation of
! extinction without the approximation. See Banerjee, Pudritz and Anderson 2006 for
! example of this type of extinction approximation.

    lam_j = 1.0/(rho_in/ndens_in) * sqrt(5./3.*pi*he_boltz*temp_in/ndens_in / Newton)
    Av    = lam_j * ndens_in * dust_gas_ratio / N_H0 ! 1% dust to gas ratio. - JW
    f_ext = exp(-3.5*Av)

return
end function

subroutine secant_method(f, dt, xlow, xhigh, xacc, root)

implicit none
integer, parameter   :: dp=KIND(1d0)
real(dp), intent(in)  :: dt       ! Solution at this step, timestep.
real(dp), intent(in)  :: xlow, xhigh ! initial brackets for the root.
real(dp), intent(in)  :: xacc        ! Accuracy between this root and previous root.
real(dp), intent(out) :: root        ! The answer!

real(dp)              :: x0, x1, del_x ! Old x, new x and step between.
real(dp)              :: y0, y1, tmp   ! Solutions in order from the iteration.
integer, parameter    :: max_iter=1000 ! Maximum allowed tries.
integer               :: n_tries

real(dp), parameter   :: kB = 1.3807d-16, mp = 1.6726d-24

interface
    function f(x_par, y_par)
    implicit none
    integer, parameter   :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface


n_tries = 0

! Make sure the bracket is in the correct order.
!if (xlow < xhigh) then
!    x0 = xlow
!    x1 = xhigh
!else
!    x0 = xhigh
!    x1 = xlow
!end if
x0 = xlow
x1 = xhigh

! Assume higher is closer to zero...
y0      = one_step(dt, x0, xlow, f) ! This just returns dt*dy/dt, so it shows which way the derivative is going.
!y0      = y0 + x0                   ! This makes it an Explicit Euler step.
y1      = one_step(dt, x1, xlow, f) ! Here look at the root from a trial Implicit Euler step from the old value.

!if (y0*y1 .gt. 0.0) &
!        y0 = 

!root    = x0
!oldroot = x1
!print*, "Temp =", y0*mp/kB*1.3d0, "ei =", y0
!print*, "Temp =", y1*mp/kB*1.3d0, "ei =", y1
! Check this assumption, and if false swap things around.
if (abs(y0) < abs(y1)) then
!    root    = x0
!    oldroot = x1 
    tmp     = y0
    y0      = y1
    y1      = tmp
    tmp     = x0
    x0      = x1
    x1      = tmp
end if

del_x = x1-x0

!print*, "Temp =", y1*mp/kB*1.3d0

do while ((abs(del_x) > xacc) .and. (y1 /= 0.0d0)) ! Check for convergence.

    n_tries = n_tries + 1

!    oldroot = root
    !x1      = (y1*xl - y0*xh) / (y1 - y0)
    del_x   = y1*(x0 - x1) / (y1 - y0)
    x0      = x1
    x1      = x1 + del_x
    y0      = y1
    y1      = one_step(dt, x1, xlow, f)

!    print*, "Temp =", x1*mp/kB*1.3d0* (5.0_dp/3.0_dp - 1.0_dp) 

!    if (x1*y1 >= 0.0) then
!        xl = root

    if (n_tries > max_iter) then
        print*, "[secant_method]: Max tries exceeded.", n_tries
        print*, "[secant_method]: temp =", y1*mp/kB*1.3d0
        print*, "[secant_method]: old bracket =", y0
        print*, "[secant_method]: new bracket =", y1
        print*, "[secant_method]: delta x =", del_x
        call flush(6)
        stop
    end if

end do

!print*, "n_tries =", n_tries

root = x1

end subroutine secant_method

subroutine falsepos_method(f, dt, xlow, xhigh, err_tol, root)

implicit none
integer, parameter    :: dp=KIND(1d0)
real(dp), intent(in)  :: dt       ! Solution at this step, timestep.
real(dp), intent(in)  :: xlow, xhigh ! initial brackets for the root.
real(dp), intent(in)  :: err_tol     ! Accuracy required in answer.
real(dp), intent(out) :: root        ! The answer!

real(dp)              :: x0, x1, del_x ! Old x, new x and step between.
real(dp)              :: y0, y1, tmp   ! Solutions in order from the iteration.
real(dp)              :: xacc          ! Accuracy between old root and new root.
integer, parameter    :: max_iter=1000 ! Maximum allowed tries.
integer               :: n_tries

real(dp), parameter   :: kB = 1.3807d-16, mp = 1.6726d-24

interface
    function f(x_par, y_par)
    implicit none
    integer, parameter   :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface


n_tries = 0

! Make sure the values bracket zero.

x0 = xlow
x1 = xhigh

y0      = one_step(dt, x0, xlow, f)
y1      = one_step(dt, x1, xlow, f)

if (y1*y0 > 0.0d0) then
    print*, "False position method requires that zero be bracketed!"
    call flush(6)
    stop
end if

! We assume y0 is neg and y1 pos.
! Check this assumption, and if false swap things around.
if (y1 < 0.0d0) then
    tmp     = y0
    y0      = y1
    y1      = tmp
    tmp     = x0
    x0      = x1
    x1      = tmp
end if

del_x = x1-x0

xacc = del_x

!print*, "Temp =", y1*mp/kB*1.3d0

do while ((abs(xacc) > err_tol) .and. (y1 /= 0.0d0)) ! Check for convergence.

    n_tries = n_tries + 1

    root    = x0 + del_x*y0 / (y0-y1)
    tmp     = one_step(dt, root, xlow, f)

    if (tmp < 0.0d0) then
        xacc  = x0 - root
        x0    = root
        y0    = tmp
    else
        xacc  = x1 - root
        x1    = root
        y1    = tmp
    end if
    
    del_x = x1-x0

!    print*, "Temp =", x1*mp/kB*1.3d0* (5.0_dp/3.0_dp - 1.0_dp) 

!    if (x1*y1 >= 0.0) then
!        xl = root

    if (n_tries > max_iter) then
        print*, "[secant_method]: Max tries exceeded.", n_tries
        print*, "[secant_method]: temp =", y1*mp/kB*1.3d0
        print*, "[secant_method]: old bracket =", y0
        print*, "[secant_method]: new bracket =", y1
        print*, "[secant_method]: delta x =", del_x
        call flush(6)
        stop
    end if

end do

!print*, "n_tries =", n_tries
return
!root = x1

end subroutine falsepos_method

function one_step(dx, y, yold, dydx)

! Take one backward Euler step along the solution curve.
! For use with implicit solution methods like midpoint, secant or
! approximate Newton's methods to find the y that makes y-f(y) = 0.

implicit none
integer, parameter   :: dp=KIND(1d0)
real(dp), intent(in)  :: dx, y, yold
real(dp)              :: one_step

interface
    function dydx(x_par, y_par)
    implicit none
    integer, parameter    :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: dydx
    end function dydx
end interface

one_step = y - dx*dydx(dx, y) - yold

end function one_step


function better_step(dx, y, yold, dydx)

! This implicit eqn uses linear approx (like Newtons Method) to
! set up the equation to be solve via a root finding method.
! See eqn 17.5.28 in Numerical Reciepes.
! For use with implicit solution methods like midpoint, secant or
! approximate Newton's methods to find the y that makes y-f(y) = 0.

implicit none
integer, parameter   :: dp=KIND(1d0)
real(dp), intent(in)  :: dx, y, yold
real(dp)              :: better_step

interface
    function dydx(x_par, y_par)
    implicit none
    integer, parameter    :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: dydx
    end function dydx
end interface

!better_step = y - dx*dydx(dx,yold) * (1.0-dx*(dydx(dx,y)-dydx(dx,yold))/max(1e-50,(y-yold)))**(-1.0) - yold
better_step = y - dx*dydx(dx,y) * (1.0-dx*(dydx(dx,y)-dydx(dx,yold))/max(1d-50,(y-yold)))**(-1.0) - yold


end function better_step

function zero ( a, b, a_old, machep, t, f, dt )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The interval [A,B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6 * MACHEPS * abs ( C ) + 2 * T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the change of 
!    sign interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real ( kind = 8 ) ZERO, the estimated value of a zero of
!    the function F.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) e
  !real ( kind = 8 ) f
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fc
  real ( kind = 8 ) m
  real ( kind = 8 ) machep
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) sa
  real ( kind = 8 ) sb
  real ( kind = 8 ) t
  real ( kind = 8 ) tol
  real ( kind = 8 ) zero
  
  real (kind = 8), intent(in) :: dt, a_old
  integer, parameter    :: max_iter=10000 ! Maximum allowed tries.
  integer               :: n_tries
  
interface
    function f(x_par, y_par)
    implicit none
    integer, parameter   :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface

  
!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = one_step(dt, sa, a_old, f) ! f ( sa )
  fb = one_step(dt, sb, a_old, f) ! f ( sb )
!  fa = better_step(dt, sa, a_old, f) ! f ( sa )
!  fb = better_step(dt, sb, a_old, f) ! f ( sb )

  c = sa
  fc = fa
  e = sb - sa
  d = e
  
  n_tries = 0
    
  if ( (fa > 0.0 .and. fb > 0.0) .or. &  
       (fa < 0.0 .and. fb < 0.0)) then
    print*, "Zero must be bracketed in zbrent!."
    write(*,'(A,X,2ES12.3)') "f(a), f(b) =", fa, fb
    call flush(6)
    stop
  endif

  do

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00 * machep * abs ( sb ) + t
    m = 0.5D+00 * ( c - sb )

    if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
      exit
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb / fa

      if ( sa == c ) then

        p = 2.0D+00 * m * s
        q = 1.0D+00 - s

      else

        q = fa / fc
        r = fb / fc
        p = s * ( 2.0D+00 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0D+00 ) )
        q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. &
        p < abs ( 0.5D+00 * s * q ) ) then
        d = p / q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb + d
    else if ( 0.0D+00 < m ) then
      sb = sb + tol
    else
      sb = sb - tol
    end if

    fb = one_step(dt, sb, a_old, f) !f ( sb )
!    fb = better_step(dt, sb, a_old, f) !f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb - sa
      d = e
    end if
    
    n_tries = n_tries + 1
    
    if (n_tries .gt. max_iter) then
        print*, "[zero]: Max tries exceeded. Exiting."
        call flush(6)
        stop
    end if

  end do

  zero = sb

  return
end function


function ridder( a, b, t, f, dt )

  implicit none
  integer, parameter    :: dp=KIND(1d0)
  real(dp), intent(in) :: a, b, t, dt

  integer               :: n_tries
  integer, parameter    :: max_iter=10000 ! Maximum allowed tries.

  real(dp) :: ridder, xa, xb, fa, fb, xm, fm, sq_root, xnew, fnew
  
interface
    function f(x_par, y_par)
    implicit none
    integer, parameter   :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: f
    end function f
end interface

  n_tries = 0
  fa = one_step(dt, a, a, f)
  fb = one_step(dt, b, a, f)
  
  if ( (fa > 0.0 .and. fb < 0.0) .or. ( fa < 0.0 .and. fb > 0.0) ) then
  
    xa = a
    xb = b
    ridder = -9d99
    
    do
    
      xm = 0.5d0*(xa+xb)
      fm = one_step(dt, xm, a, f)
      sq_root = sqrt(fm*fm-fa*fb)
      
      if (sq_root == 0.0d0) then
        ridder = xm
        exit
      end if
      
      if (fa >= fb) then
        xnew = xm + (xm-xa)*(fm/sq_root)
      else
        xnew = xm - (xm-xa)*(fm/sq_root)
      end if
      
      if (abs(xnew - ridder) <= t) exit
      
      ridder = xnew
      fnew = one_step(dt, xnew, a, f)
      if (fnew == 0.0d0) exit
      
      
      if (sign(fm,fnew) /= fm) then
        xa = xm
        fa = fm
        xb = xnew
        fb = fnew
      else if (sign(fa,fnew) /= fa) then
        xb = xnew
        fb = fnew
      else if (sign(fb,fnew) /= fb) then
        xa = xnew
        fa = fnew
      else
        print*, "[ridder]: Something definitely broke!"
      end if
      
      if (abs(xb-xa) < t) exit
      
      n_tries = n_tries + 1
      
      if (n_tries >= max_iter) then
        print*, "[ridder]: Max tries exceeded!"
        exit
      end if
    
      
    end do
    

  else if (fa == 0.0d0) then
    ridder = a
  else if (fb == 0.0d0) then
    ridder = b
  else
    print*, "[ridder]: Initial values must bracket zero!"
    print*, "[ridder]: fa, fb =", fa, fb
  end if
    
    
end function

subroutine zbrac(func,x1,x2,xold,dt,success)

! Bracket finding algorithm from Numerical Reciepes
! Press et. al. 2007.

use cool_vars, only : TtoEI

    implicit none
    
    integer, parameter      :: dp=KIND(1d0)
    real(dp), intent(inout) :: x1,x2
    real(dp), intent(in)    :: xold,dt
    logical, intent(inout)  :: success
    real(dp), parameter     :: FACTOR=1.6d0
    integer, parameter      :: NTRY=5000
    integer                 :: j
    real(dp)                :: f1,f2
      
interface
    function func(x_par, y_par)
    implicit none
    integer, parameter   :: dp=KIND(1d0)
    real(dp), intent(in) :: x_par, y_par
    real(dp)             :: func
    end function func
end interface
      
      if (x1.eq.x2) then
        print*, '[zbrac]: Bad initial range in zbrac!'
      end if
      f1=one_step(dt, x1, xold, func)
      !write(*,'(A,X,ES13.3)') 'f1=', f1
      f2=one_step(dt, x2, xold, func)
      !write(*,'(A,X,ES13.3)') 'f2=', f2
      success=.true.
      do j=1, NTRY
        if (f1*f2 .lt. 0.0d0) & 
          return
        if (abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          !write(*,'(A,X,ES13.3)') 'x1=', x1
          if (x1 < 0.0d0) x1 = TtoEI*1.0d0
          f1=one_step(dt, x1, xold, func)
          !write(*,'(A,X,ES13.3)') 'f1=', f1
        else
          x2=x2+FACTOR*(x2-x1)
          !write(*,'(A,X,ES13.3)') 'x2=', x2
          if (x2 < 0.0d0) x2 = TtoEI*1.0d0
          f2=one_step(dt, x2, xold, func)
          !write(*,'(A,X,ES13.3)') 'f2=', f2
        endif
      end do
      success=.false.
      return
end subroutine zbrac

function get_mu(xHp_in, ndens_in)
real :: get_mu
real :: ll, ul
real :: log_ndens
real, parameter  :: mp = 1.6726e-24, lower_ion_limit = 1d-3
real, intent(in) :: xHp_in, ndens_in

    if (xHp_in .gt. lower_ion_limit) then ! Probably no H2 if more than some % ionized.
        ll   = lower_ion_limit
        ul   = 1.0 !  Linearly interp on X
        get_mu = lin_interp(ll, 1.3d0, ul, 0.61d0, xHp_in)
    else ! This regime made to look like Micic et al 2012 Fig 5, where H2 frac looks log linear in ndens.
        if (ndens_in .le. 1e4) then
            if (ndens_in .ge. 1e1) then
                ll   = 1 ! log lower lim
                ul   = 4 ! log upper lim
                log_ndens = log10(ndens_in)
                get_mu = lin_interp(ll, 1.3d0, ul, 24.0d0/11.0d0, log_ndens)
            else
                get_mu = 1.3d0
            end if
        else
            get_mu = 24.0d0/11.0d0
        end if
    end if
    
    return
end function get_mu

function lin_interp(x1, y1, x2, y2, x)
  real*8, intent(in) :: x1, x2, y1, y2, x
  
  real*8 lin_interp
  real*8 ms

  ms = (y2-y1)/(x2-x1)
  lin_interp = ms*(x-x1) + y1
  
  if (x .lt. x1 .or. x .gt. x2) print*, "[heatCool:lin_interp:] WARNING! &
                                          x is out of bounds and you're now &
                                          extrapolating (not interpolating!)"
  return
end function lin_interp

end module heatCool
