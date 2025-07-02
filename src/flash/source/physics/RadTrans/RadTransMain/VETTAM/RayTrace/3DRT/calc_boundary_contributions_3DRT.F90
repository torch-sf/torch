!===============================================================================
!
! Subroutine: calc_boundary_contributions_3DRT
!
!===============================================================================
!
!
!===============================================================================

subroutine calc_boundary_contributions_3DRT ( &
     myPe, cutList, nrOfCuts,             &
     regDestX, regDestY, regDestZ,        &
     xMin, yMin, zMin,                    &
     xMax, yMax, zMax,                    &
     regFactX, regFactY, regFactZ,        &
     dirX, dirY, dirZ,                    &
     maxLevel,                            &
     blockCoordAll,                       &
     faceValueAllMean,                    &
     refLevelAll, nrOfPe,                 &
     leafListAll, b, maxNrOfLeafBlocks,   &
     x_periodic,y_periodic,z_periodic,    &
     irradiation                          )
!
  use tree,  only: maxblocks, lrefine
  use dBase, only: &
       nxb, nyb, nzb, nguard,    &
       maxcells, k2d, k3d, ndim, &
       dBasePropertyInteger,     &
       dBaseKeyNumber,           &
       dBaseGetData
!
  implicit none
!
  integer, intent(in), dimension(maxblocks)         :: cutList
  integer, intent(in)                               :: nrOfCuts
  real,    intent(in)                               :: regDestX
  real,    intent(in)                               :: regDestY
  real,    intent(in)                               :: regDestZ
  real,    intent(in)                               :: regFactX
  real,    intent(in)                               :: regFactY
  real,    intent(in)                               :: regFactZ
  real,    intent(in)                               :: dirX
  real,    intent(in)                               :: dirY
  real,    intent(in)                               :: dirZ
  real,    intent(in)                               :: xMin,yMin,zMin
  real,    intent(in)                               :: xMax,yMax,zMax
  integer, intent(in)                               :: maxLevel
  real,    intent(in), dimension(:,:)               :: blockCoordAll
  integer, intent(in), dimension(:)                 :: refLevelAll
  real,    intent(in), dimension(:,:)               :: faceValueAllMean
  integer, intent(in)                               :: nrOfPe
  integer, intent(in)                               :: myPe
  integer, intent(in)                               :: b
  integer, intent(in)                               :: maxNrOfLeafBlocks
  logical, intent(in)                               :: x_periodic
  logical, intent(in)                               :: y_periodic
  logical, intent(in)                               :: z_periodic
  integer, intent(in), dimension(maxblocks*nrOfPe)  :: leafListAll
  real,    intent(inout)                            :: irradiation 
!
    !
  return
  !
end subroutine calc_boundary_contributions_3DRT
