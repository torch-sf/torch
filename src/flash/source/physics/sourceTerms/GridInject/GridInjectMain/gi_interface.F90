!!****ih* source/physics/sourceTerms/GridInject/GridInjectMain/gi_interface
!!
!! This is an interface module for the GridInject unit that defines some
!! private interfaces
!!***

Module gi_interface

  interface gi_distance
    function gi_distance (xx, yy, zz, xloc, yloc, zloc)
      implicit none
      real, intent(IN) :: xx, yy, zz, xloc, yloc, zloc
      real :: gi_distance
    end function gi_distance
  end interface

  interface
    subroutine gi_distanceVector (xx, yy, zz, xloc, yloc, zloc, dx, dy, dz)
      implicit none
      real, intent(IN) :: xx, yy, zz, xloc, yloc, zloc
      real, intent(OUT) :: dx, dy, dz
    end subroutine gi_distanceVector
  end interface

  interface
    subroutine gi_normal_rand (mean, std_dev, randnum)
      implicit none
      real, intent(IN) :: mean, std_dev
      real, intent(OUT) :: randnum
    end subroutine gi_normal_rand
  end interface

  interface
    subroutine gi_overlap (ishp, rad, center, cell_bot, cell_top, nsteps, overlap_vol)
      implicit none
      integer, intent(IN) :: ishp, nsteps
      real, intent(IN)    :: rad, center(3), cell_bot(3), cell_top(3)
      real, intent(OUT)   :: overlap_vol
    end subroutine gi_overlap
  end interface

end Module gi_interface
