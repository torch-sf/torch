Module pt_windInterface

  interface
    subroutine inject_direct(loc, injectMass, injectVelocity, twind, dt)

      implicit none
      real, intent(in) :: loc(3)
      real, intent(in) :: injectMass, injectVelocity, twind
      real, intent(inout) :: dt
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

end Module pt_windInterface
