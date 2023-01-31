!! In Particles/localAPI/
subroutine pt_initVoronoiPositions (success,updateRefine)


  implicit none

  integer, INTENT(OUT) :: updateRefine
  logical,intent(INOUT) :: success

  updateRefine = .false.
  success = .true. ! DEV: returns true because this stub creates no particles,
                   ! therefore all of those zero particles were created successfully
  return
  !----------------------------------------------------------------------
end subroutine pt_initVoronoiPositions
