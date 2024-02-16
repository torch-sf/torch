Module pt_windInterface

  interface
    subroutine inject_direct(loc_in, injectMassIn, injectVelocityIn, twind, dt, bgDens)
    ! Updated 20240215 to accurately reflect inject_direct call -SA

      implicit none
      real, intent(in) :: loc_in(3)
      real, intent(in) :: injectMassIn, injectVelocityIn, twind
      real, intent(inout) :: dt
      real, intent(inout) :: bgDens
    end subroutine inject_direct
  end interface

  interface
    subroutine overlap(ishp, rad, center, cell_bot, cell_top, nsteps, overlap_vol)

      implicit none
        integer, intent(in)     :: ishp, nsteps
        real, intent(in)        :: rad, center(3), cell_bot(3), cell_top(3)
        real, intent(out)       :: overlap_vol
    end subroutine overlap
  end interface

  interface
    subroutine sphere_and_cell_frac(vol4,R4,x4,y4,z4,cellsize4)
        implicit none
        real,intent(in) :: R4,x4,y4,z4,cellsize4
        real(kind=8) :: R,x,y,z,cellsize,vol
        real,intent(out) :: vol4
    end subroutine sphere_and_cell_frac
  end interface

end Module pt_windInterface
