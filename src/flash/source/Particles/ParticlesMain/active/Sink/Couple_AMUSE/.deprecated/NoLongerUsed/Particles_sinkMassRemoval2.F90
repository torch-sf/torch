!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Particles_sinkMassRemoval
!!
!! Take the mass in a sink out of it.
!! Only when the free fall time has passed.
!!
!! K thanks. - Josh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Particles_sinkMassRemoval(dt, mass, x, y, z, vx, vy, vz, &
                                     ang_x, ang_y, ang_z)

use Particles_sinkData
use Driver_data, ONLY: dr_simTime
use RuntimeParameters_interface, ONLY : RuntimeParameters_get
use PhysicalConstants_interface, ONLY : PhysicalConstants_get
use Grid_interface, ONLY : Grid_getMinCellSize

#include "Flash.h"
#include "constants.h"

implicit none

real, intent(in) :: dt
real, dimension(10), intent(out) :: mass, x, y, z, vx, vy, vz, ang_x, ang_y, ang_z

integer :: p
real, save :: free_fall_time, delx, crit_density, newton
real :: time, sink_time

logical :: first_call = .true.

if (first_call) then

    call RuntimeParameters_get("sink_density_thresh", crit_density)
    call Grid_getMinCellSize(delx)
    call PhysicalConstants_get("newton", newton)
    
    free_fall_time = sqrt(3.0 * PI / (32.0*newton*crit_density))

	first_call=.false.

end if

mass = 0.0; x = 0.0; y = 0.0; z = 0.0
vx = 0.0; vy = 0.0; vz = 0.0
ang_x = 0.0; ang_y = 0.0; ang_z = 0.0

time = dr_simTime

do p=1, localnp

    sink_time = particles_local(CREATION_TIME_PART_PROP, p)
    
    if ((time - sink_time) .lt. free_fall_time) then
    
		print*, "We're still cooking. Timer says ", time-sink_time
		cycle

    else
		
		print*, "Now we're taking mass from this cell!"
		mass = particles_local(MASS_PART_PROP, p)
		x = particles_local(POSX_PART_PROP, p)
			y = particles_local(POSY_PART_PROP, p)
			z = particles_local(POSZ_PART_PROP, p)
		vx = particles_local(VELX_PART_PROP, p)
		vy = particles_local(VELY_PART_PROP, p)
		vz = particles_local(VELZ_PART_PROP, p)
		ang_x = particles_local(X_ANG_PART_PROP, p)
		ang_y = particles_local(Y_ANG_PART_PROP, p)
		ang_z = particles_local(Z_ANG_PART_PROP, p)
	
		particles_local(MASS_PART_PROP, p) = 0.0
		particles_local(POSX_PART_PROP, p) = 0.0
			particles_local(POSY_PART_PROP, p) = 0.0
			particles_local(POSZ_PART_PROP, p) = 0.0
		particles_local(VELX_PART_PROP, p) = 0.0
		particles_local(VELY_PART_PROP, p) = 0.0
		particles_local(VELZ_PART_PROP, p) = 0.0
		particles_local(X_ANG_PART_PROP, p) = 0.0
		particles_local(Y_ANG_PART_PROP, p) = 0.0
		particles_local(Z_ANG_PART_PROP, p) = 0.0
	
    end if
    
end do

end subroutine Particles_sinkMassRemoval
