!!____    ____ ____ _   _ ____ 
!!|___ __ |__/ |__|  \_/  [__  
!!|       |  \ |  |   |   ___]                                                    '    
!!written by C. Baczynski, 2012-2014


!! Description:
!!   initialisation routine for ray tracing
!!
!! Input: 
!!   restart: simulation restart switch

subroutine Particles_rayInit(restart)

  use Particles_rayData
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm, Driver_getNumProcs

  use HEALPixModule , ONLY :  mk_pix2xy, mk_xy2pix, mk_xy2pix1, pix2vec_nest
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_data, ONLY :  gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax, &
                         gr_domainBC

  use pt_rayAsyncComm, ONLY : ph_globalMe, ph_globalComm, ph_globalNumProcs, & 
                              ph_meshMe, ph_meshComm, ph_meshNumProcs, ph_numNeigh, & 
                              ph_size, ph_rank, ph_raysToBundle, ph_CommCheckInterval
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"
#include "Flash_mpi.h"

  logical, intent(in) :: restart
  integer :: ierr
  integer :: commoffset 

!! HEALPIx specific initialisation
  call mk_pix2xy()
  call mk_xy2pix()
  call mk_xy2pix1()

!! raytracing options
  call RuntimeParameters_get ("ph_sampling", ph_sampling)
  call RuntimeParameters_get ("ph_initHPlevel", ph_initHPlevel)
  call RuntimeParameters_get ("ph_inBlockSplit", ph_inBlockSplit)
  call RuntimeParameters_get ("ph_rotRays", ph_rotRays)
  call RuntimeParameters_get ("ph_periodicBoxL", ph_periodicBoxL)
  call RuntimeParameters_get ("ph_maxNRays", ph_maxNRays)
  call RuntimeParameters_get ("ph_raysToBundle", ph_raysToBundle)
  call RuntimeParameters_get ("ph_CommCheckInterval", ph_CommCheckInterval)
  call RuntimeParameters_get ("ph_xPeriodic", xperiodic)
  call RuntimeParameters_get ("ph_yPeriodic", yperiodic)
  call RuntimeParameters_get ("ph_zPeriodic", zperiodic)
  call RuntimeParameters_get ("useRadTransfer", useRadTransfer)
  call RuntimeParameters_get ("ph_radPressure", ph_radPressure)
  call RuntimeParameters_get ("early_term_FUV", early_term_FUV)
  call RuntimeParameters_get ("ph_EUVonDust", ph_EUVonDust)
  call PhysicalConstants_get("speed of light",speedoflight)

  ph_locsampling = ph_sampling

! assign indexes in rays array
! direct MPI
  inion  = 1  ! real
  ieion  = 2  ! real
! from global source list and healpix + 1 intersection calc.
  irad   = 3  ! real
! direct MPI
  ihnum  = 4  ! double integer = real?
! from healpix reconstruction
  ivelx  = 5  ! real
  ively  = 6  ! real
  ivelz  = 7  ! real
! from global source list 
  iposx  = 8 ! real
  iposy  = 9 ! real
  iposz  = 10 ! real
! cross sections 
  isigh  = 11 ! real

! packed 2 bytes
  iblk   = 1  ! integer
! from MPI
  iproc  = 2  ! integer
! packed 1 byte
  ihlev  = 3  ! integer
! not transported used as packed info buffer
  istpd  = 4  ! integer
! source ID 1 byte
  isid   = 5 
! radiation bin (for now either >13.6 eV (ionizing) or 5.6-13.6 eV (PE) - JW
  itype  = 6 ! integer

! MPI buffer array positions
! for H ionisation
  itnion  = 1
  ithnum  = 2
! healpix number
  itinfo  = 3
! radius
  itrad   = 4
! source id
  itid    = 5
! current block
  itblk   = 6

  glDX    = gr_imax - gr_imin
  glDY    = gr_jmax - gr_jmin
  glDZ    = gr_kmax - gr_kmin

! define additional transport property indices 
  iteion  = 7

  itvelx  = 8
  itvely  = 9
  itvelz  = 10

  itposx  = 11
  itposy  = 12
  itposz  = 13

  itsigh  = 14
  itstpd  = 15
  ittype  = 16 ! radiation bin info - JW

! processor id always known so no last property
! all properties, -1 for processor id
  ph_transProp = ph_numIntProp + ph_numRealProp - 1

  ph_periodic = .false.

  ! are any boundaries periodic?
  if(xperiodic .or. yperiodic .or. zperiodic) then
    ph_periodic = .true.    

! only x periodic
    if(xperiodic .and. .not. yperiodic .and. .not. zperiodic) then
      ph_BCcase = 1
    endif

! only y periodic
    if(yperiodic .and. .not. xperiodic .and. .not. zperiodic) then
      ph_BCcase = 2
    endif

! only z periodic
    if(zperiodic .and. .not. yperiodic .and. .not. xperiodic) then
      ph_BCcase = 3 
    endif

! x and y 
    if(xperiodic .and. yperiodic .and. .not. zperiodic) then
      ph_BCcase = 4
    endif

! x and z
    if(xperiodic .and. zperiodic .and. .not. yperiodic) then
      ph_BCcase = 5
    endif

! y and z
    if(yperiodic .and. zperiodic .and. .not. xperiodic) then
      ph_BCcase = 6
    endif

! x and y and z
    if(xperiodic .and. yperiodic .and. zperiodic) then
      ph_BCcase = 7
    endif
  endif

  xhydroper = .false.
  yhydroper = .false.
  zhydroper = .false.

! check if hydro bounds are different from ray bounds
  if(gr_domainBC(LOW,IAXIS) .eq. PERIODIC .and. .not. xperiodic) then
    xhydroper = .true.
  endif

  if(gr_domainBC(LOW,JAXIS) .eq. PERIODIC .and. .not. yperiodic) then
    yhydroper = .true.
  endif

  if(gr_domainBC(LOW,KAXIS) .eq. PERIODIC .and. .not. zperiodic) then
    zhydroper = .true.
  endif

! assorted MPI init
! setup communication 
  call Driver_getMype        (GLOBAL_COMM, ph_globalMe       )
  call Driver_getComm        (GLOBAL_COMM, ph_globalComm     )
  call Driver_getNumProcs    (GLOBAL_COMM, ph_globalNumProcs )
  call Driver_getMype        (  MESH_COMM, ph_meshMe         )
  call Driver_getComm        (  MESH_COMM, ph_meshComm       )
  call Driver_getNumProcs    (  MESH_COMM, ph_meshNumProcs   )

! get rank
  call MPI_Comm_rank(ph_meshComm, ph_rank, ierr)

! get size
  call MPI_Comm_size(ph_meshComm, ph_size, ierr)

  if(ph_meshMe==MASTER_PE) print*,'pew pew rays in space'

! allocate all data at start
  allocate(raysIntProp (ph_numIntProp,  ph_maxNRays))
  allocate(raysRealProp(ph_numRealProp, ph_maxNRays))

  return
end subroutine Particles_rayInit
