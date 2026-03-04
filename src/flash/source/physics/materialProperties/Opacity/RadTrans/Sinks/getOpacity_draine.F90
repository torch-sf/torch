!!****f* source/physics/materialProperties/Opacity/RadTrans/Sink
!!
!! NAME
!!
!!  Opacity_sink
!!
!!
!! SYNOPSIS
!!
!!  call getOpacity_draine()
!!
!! DESCRIPTION
!!
!! Read opacities for direct radition from sink particles using the Draine
!! 2003 opacity table. This assumes that the stored table is present in the
!! same directory under the name 'DraineStellarOpacities.dat', and contains 
!! linearly separated values of Tstar and its corresponding Planck mean opacity
!! obtained from the Draine 2003 opacity table. 
!!
!! ARGUMENTS
!!
!!  sink_Tstar (in) - Colour temperature of sink particle
!!  opac_sink (out) - Planck mean opacity appropriate to this temperature
!!
!! PARAMETERS
!!
!! AUTHOR
!!  Shyam Menon (2021)
!!
!!***

subroutine getOpacity_draine(sink_Tstar,opac_sink)
#include "Flash.h"
#include "constants.h"
#ifdef RAYTRACE_3DRT
  use Particles_sinkData, ONLY : useSinkParticles
#endif
  implicit none
  real, intent(in) :: sink_Tstar
  real, intent(out) :: opac_sink
  real, save, allocatable, dimension(:) :: Tstar_table, kappa_table
  real, allocatable, dimension(:,:) :: fileoutput
  real, save :: Tstar_min, Tstar_max, dTstar
  real :: Tstar
  logical, save :: first_call = .true.
  character(len=50) :: filename = 'DraineStellarOpacities.dat'
  integer :: i,iostat, index, istat
  integer, save :: nlines = 0
  !TODO: Dust-to-Gas ratio hard-coded. Does this need to change?
  real, parameter :: dustToGasRatio = 0.01

! Return if Raytrace not present or if useSinkParticles is set to False
#ifdef RAYTRACE_3DRT
  if(.not. useSinkParticles) return
#else 
  return
#endif

  ! First call
  if (first_call) then

    ! Find number of lines in the file
    OPEN (unit=07,file=trim(filename),status='old',access='sequential')
    do 
      read(07,*,END=10)
      nlines = nlines+1
    end do
    10 CLOSE(unit=07)

    ! Allocate memory for variables
    allocate(fileoutput(2,nlines),stat=istat)
    if(istat .ne. 0) call Driver_abortFlash("[Sink_Opac]: Could not allocate fileoutput")
    allocate(Tstar_table(nlines),stat=istat)
    if(istat .ne. 0) call Driver_abortFlash("[Sink_Opac]: Could not allocate Tstar_table")
    allocate(kappa_table(nlines),stat=istat)
    if(istat .ne. 0) call Driver_abortFlash("[Sink_Opac]: Could not allocate kappa_table")

    !Read table
    OPEN (unit=07,file=trim(filename),status='old',access='sequential')

    do i = 1, nlines
      read(07,*) fileoutput(:,i)
    end do
    CLOSE(unit=07)

    Tstar_table = fileoutput(1,:)
    kappa_table = fileoutput(2,:)

    ! Save min/max T and spacing in Tstar table
    Tstar_min = Tstar_table(1)
    Tstar_max = Tstar_table(nlines)
    dTstar = Tstar_table(2)-Tstar_table(1)

    first_call = .false.
  endif

  ! limit Tstar to be within table margins
  Tstar = sink_Tstar
  if (Tstar .le. Tstar_min) Tstar = Tstar_min
  if (Tstar .ge. Tstar_max) Tstar = Tstar_max - dTstar ! so we don't run into an index problem below

  ! Assume linearly spaced bins from Tstar_min with spacing dTstar
  index = (Tstar - Tstar_min) / dTstar + 1

  ! Linearly interpolate from 2 closest stored values in table
  opac_sink = (kappa_table(index)*(Tstar_table(index+1)-Tstar) + &
              kappa_table(index+1)*(Tstar-Tstar_table(index))) / &
              (Tstar_table(index+1)-Tstar_table(index))

  opac_sink = opac_sink * dustToGasRatio

end subroutine
