!!****if* source/Particles/ParticlesMain/active/Sink/pt_sinkInterface
!!
!! NAME
!!
!!  pt_sinkInterface
!!
!!***

module pt_sinkInterface
  implicit none
  
	interface
      subroutine pt_sinkAccelSinksOnSinks(local_min_radius, local_max_accel)
        real, intent(out) :: local_min_radius, local_max_accel
      end subroutine pt_sinkAccelSinksOnSinks
	end interface

	interface
    	subroutine pt_sinkAccelGasOnSinks()
    	end subroutine pt_sinkAccelGasOnSinks
  	end interface

	interface
 		subroutine pt_sinkCorrectForPeriodicBCs(distx, disty, distz)
   		real, intent(INOUT)  :: distx, disty, distz
 		end subroutine pt_sinkCorrectForPeriodicBCs
	end interface

  	interface
    	function pt_sinkCreateParticle(x, y, z, pt, block_no, MyPE)
      		real, intent(IN)    :: x, y, z, pt
      		integer, intent(IN) :: block_no, MyPE
      		integer :: pt_sinkCreateParticle
    	end function pt_sinkCreateParticle
  	end interface
  
  	interface
		subroutine pt_sinkDumpParticles(simtime, outfilename)
	  		real, intent(IN) :: simtime
	  		character(len=*), intent(IN) :: outfilename	
		end subroutine pt_sinkDumpParticles
			! JC - Old Code version, outfilename was static and created within the routine.
			!subroutine pt_sinkDumpParticles(simtime)
			!	real, intent(IN)   :: simtime
			!end subroutine pt_sinkDumpParticles
  	end interface

  	interface
    	subroutine pt_sinkFindList(x, y, z, rad, create_part, pindex_found, np_found)
      		use Particles_sinkData, only: maxsinks
      		real, intent(IN) :: x, y, z, rad
      		logical, intent(IN) :: create_part
      		integer, dimension(maxsinks), intent(OUT) :: pindex_found
      		integer, intent(OUT) :: np_found
    	end subroutine pt_sinkFindList
  	end interface
  
  	interface
	  	subroutine pt_sinkGatherGlobal(propinds, nprops)
	    	integer, dimension(:), intent(in), optional :: propinds
	    	integer, intent(in), optional :: nprops
		end subroutine pt_sinkGatherGlobal
			! JC - Old code version: new array called propinds used to know which properties to collect nad the number of entries.
			!subroutine pt_sinkGatherGlobal()
			!end subroutine pt_sinkGatherGlobal
  	end interface
  
  	interface
    	subroutine pt_sinkGetSubCycleTimeStep(dt, dt_global, local_min_radius, local_max_accel)
      		real, intent(out)  :: dt
      		real, intent(in)   :: dt_global, local_min_radius, local_max_accel
    	end subroutine pt_sinkGetSubCycleTimeStep
  	end interface
  
  	interface
    	subroutine pt_sinkMergingAfterCreation(delta_at_lrefmax)
      		real, intent(IN)    :: delta_at_lrefmax
    	end subroutine pt_sinkMergingAfterCreation
  	end interface
  
  	interface
    	subroutine pt_sinkParticleMerging(dt)
      		real, intent(IN)    :: dt
    	end subroutine pt_sinkParticleMerging
  	end interface
  
	interface
      subroutine pt_sinkEwaldCorrection(x, y, z, xcorr, ycorr, zcorr)
        real, intent(IN)  :: x, y, z
        real, intent(OUT) :: xcorr, ycorr, zcorr
      end subroutine pt_sinkEwaldCorrection
	end interface

	interface
      subroutine pt_sinkPrepareEwald()
      end subroutine pt_sinkPrepareEwald
	end interface

end module pt_sinkInterface
