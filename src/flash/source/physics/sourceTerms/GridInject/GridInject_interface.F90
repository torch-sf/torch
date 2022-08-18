!!****h* source/physics/sourceTerms/GridInject/GridInject_interface
!!
!! NAME
!!
!!  GridInject_interface
!!
!! SYNOPSIS
!!
!!   use GridInject_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Energy Deposition unit that defines its
!!  public interfaces.
!!
!!***

Module GridInject_interface

  interface
    subroutine GridInject_init ()
    end subroutine GridInject_init
  end interface

  interface
    subroutine GridInject_thermalSN (xloc, yloc, zloc, energy, mass, &
      r_min, r_max, nms, r_exp, m_exp, rho_avg)
      real, intent(IN)    :: xloc, yloc, zloc, energy, mass, r_min, r_max
      integer, intent(IN) :: nms
      real, intent(OUT)   :: r_exp, m_exp, rho_avg
    end subroutine GridInject_thermalSN
  end interface

  interface
    subroutine GridInject_kineticSN (xloc, yloc, zloc, energy, mass, &
        snap_to_grid)
      real, intent(IN)              :: xloc, yloc, zloc, energy, mass
      logical, intent(IN), optional :: snap_to_grid
    end subroutine GridInject_kineticSN
  end interface

  interface
    subroutine GridInject_wind (xloc, yloc, zloc, injectMassIn, &
        injectVelocityIn, twind, dt, bgDens, snap_to_grid)
      real, intent(IN)    :: xloc, yloc, zloc
      real, intent(IN)    :: injectMassIn, injectVelocityIn, twind, dt
      real, intent(INOUT) :: bgDens
      logical, optional, intent(IN) :: snap_to_grid
    end subroutine GridInject_wind
  end interface

  interface
    subroutine GridInject_getInjBlks (xloc, yloc, zloc, radius, injBlks, InjBlkNum)
      real, intent(IN)      :: xloc, yloc, zloc, radius
      integer, intent(OUT)  :: injBlks(MAXBLOCKS)
      integer, intent(OUT)  :: injBlkNum
    end subroutine GridInject_getInjBlks
  end interface

end Module GridInject_interface
