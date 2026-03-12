!!****if* source/physics/RadTrans/RadTransMain/useHybridChar3DRTFLD/RayTrace/3DRT/rt_init_raytrace_3drt
!!
!!  NAME 
!!
!!  rt_init
!!
!!  SYNOPSIS
!!
!!  call rt_init_raytrace_3drt()
!!
!!  DESCRIPTION 
!!    Initialize local data for each radiative transfer model
!!
!!***

#undef DEBUG_RT

subroutine rt_init_raytrace_3drt

  use RadTrans_hybridCharModule, ONLY : io, io_sour, io_conv, getAngles
  use raytrace_data, ONLY : rt_maxNrOfBoundIter, rt_nrOfAngleGroups, rt_dirX, rt_dirY, rt_dirZ, &
                      rt_healpix_nSide, rt_nPhi, rt_nTheta, allocGeneration
  use rt_data_raytrace_3drt
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_data, ONLY : gr_globalDomain, gr_domainBC
  use Driver_data, ONLY : dr_globalMe
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"
#include "constants.h"

  !
  !    Prepare log files
  !
  call prepare_log_file
  !
  !    Prepare the solid angle grid
  !
  call prepare_solid_angle_grid
  !
  !    Set the begin and end indices of the internal zones of a block.
  !    NOTE: since we use the zone corners the starting indices are set to
  !          NGUARD instead of the usual NGUARD+1
  !
  !    Check for periodic boundary conditions
  !
  if(gr_domainBC(LOW,IAXIS).eq.PERIODIC.and.gr_domainBC(HIGH,IAXIS).eq.PERIODIC) then
     x_periodic = .true.
  elseif(gr_domainBC(LOW,IAXIS).ne.PERIODIC.and.gr_domainBC(HIGH,IAXIS).ne.PERIODIC) then
     x_periodic = .false.
  else
     call Driver_abortFlash("radiation: no consistent boundary conditions in x")
  endif
  if(NDIM>=2) then
     if(gr_domainBC(LOW,JAXIS).eq.PERIODIC.and.gr_domainBC(HIGH,JAXIS).eq.PERIODIC) then
        y_periodic = .true.
     elseif(gr_domainBC(LOW,JAXIS).ne.PERIODIC.and.gr_domainBC(HIGH,JAXIS).ne.PERIODIC) then
        y_periodic = .false.
     else
        call Driver_abortFlash("radiation: no consistent boundary conditions in y")
     endif
  endif
  if(NDIM==3) then
     if(gr_domainBC(LOW,KAXIS).eq.PERIODIC.and.gr_domainBC(HIGH,KAXIS).eq.PERIODIC) then
        z_periodic = .true.
     elseif(gr_domainBC(LOW,KAXIS).ne.PERIODIC.and.gr_domainBC(HIGH,KAXIS).ne.PERIODIC) then
        z_periodic = .false.
     else
        call Driver_abortFlash("radiation: no consistent boundary conditions in z")
     endif
  endif
  !
  if(x_periodic.or.y_periodic.or.z_periodic) then
     nrOfBoundIter = rt_maxNrOfBoundIter
  else
     nrOfBoundIter = 0 
  endif
  if(nrOfBoundIter>0) then
     if(dr_globalMe.eq.MASTER_PE) then
        write(io,*) 'WARNING: rt_maxNrOfBoundIter > 0'
        write(io,*) 'Periodic Boundary conditions for the radiation transfer'
        write(io,*) 'are not properly implemented yet. We ignore it for the '
        write(io,*) 'moment and force rt_maxNrOfBoundIter to be 0.'
     endif
     rt_maxNrOfBoundIter = 0
  endif
  !
  !    Calculate the physical domain size and the total number of grid points on
  !    a regular grid that corresponds to a fully refined domain. Use these to
  !    calculate conversion factors.
  !
  domainSizeX = gr_globalDomain(HIGH,IAXIS) - gr_globalDomain(LOW,IAXIS)
  domainSizeY = gr_globalDomain(HIGH,JAXIS) - gr_globalDomain(LOW,JAXIS)
  domainSizeZ = gr_globalDomain(HIGH,KAXIS) - gr_globalDomain(LOW,KAXIS)
  !
  !    Output of some general information.
  !
  if(dr_globalMe.eq.MASTER_PE) then
     !
     write(io,*) 'Number of Angle Groups:',rt_nrOfAngleGroups
     write(io,*) 'Angles per Angle Group:',nrOfAnglesPerGroup
     write(io,*) 'Total Number of Angles:', nrOfAngles
     write(io,*) '-------------------------------------------------------------------'
     !
  endif

  allocGeneration = -1

  return ! subroutine rt_init_raytrace_3drt

  !
  ! ======================================================================
  !

contains

  !
  ! ======================================================================
  !
  subroutine prepare_log_file
    !
    implicit none
    
    character (len=80)  :: myPeChar

    !
    !    Prepare raytracer logfiles
    !
#ifndef DEBUG_RT
    if(dr_globalMe.eq.MASTER_PE) then
#endif
       write(myPeChar,'(I5.5)') dr_globalMe
       open(io,file='raytrace_'//trim(myPeChar)//'.log')
       !
       write(io,*) '======================'
       write(io,*) 'Raytrace 3DRT Logfile ' 
       write(io,*) '======================'
       !
       !    Give some geometrical informations, just to make sure.
       !
       write(io,*) '-------------------------------------------------------------------'
       write(io,*) 'Some general geometrical Information:'
       write(io,*) 'Number of Dimensions (NDIM):',NDIM
       write(io,*) 'Linear Number of Cells (NXB,NYB,NZB):',NXB,NYB,NZB
       write(io,*) 'Dimension Indictation (K2D,K3D):',K2D,K3D
       write(io,*) 'Number of Guard Cells (NGUARD):',NGUARD
       !
#ifndef DEBUG_RT
    endif
#endif
    !

    if(dr_globalMe.eq.MASTER_PE) then
       open(io_sour,file='source_function_convergence_'            &
            //trim(myPeChar)//'.out')
       !
       write(io_sour,*) '========================================='
       write(io_sour,*) 'Source Function Convergence  ' 
       write(io_sour,*) '' 
       write(io_sour,*) 'In this file, you can see the maximum    '
       write(io_sour,*) 'change in the source function.           '
       write(io_sour,*) 'This file has no effect and is for your  '
       write(io_sour,*) 'convenience only.'
       write(io_sour,*) '========================================='
       write(io_sour,*) '' 
       !
       open(io_conv,file='mean_intensity_convergence_'             &
            //trim(myPeChar)//'.out')
       !
       write(io_conv,*) '========================================='
       write(io_conv,*) 'Mean Intensity Convergence  ' 
       write(io_conv,*) '' 
       write(io_conv,*) 'In this file, you can see the maximum    '
       write(io_conv,*) 'change in the mean intensity.            '
       write(io_conv,*) 'This file has no effect and is for your' 
       write(io_conv,*) 'convenience only.'
       write(io_conv,*) '========================================='
       write(io_conv,*) ''       
    endif

    
  end subroutine prepare_log_file
   
  !
  ! ======================================================================
  !
  
  subroutine prepare_solid_angle_grid
     !
     implicit none

     integer :: iAngle
     real    :: theta, phi, dirFact(3)
     
     !
     !  Total number of directions
     !
     if(NDIM==3) then
        if(rt_healpix_nSide>0) then
           !
           ! The healpix case
           !
           nrOfAngles = 12*rt_healpix_nSide*rt_healpix_nSide
           !
           if(dr_globalMe.eq.MASTER_PE) then 
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'Using Healpix Tesselation!'
              write(io,*) 'Resolution Parameter (rt_healpix_nSide)',rt_healpix_nSide
              write(io,*) 'Total Number of Angles (nrOfAngles)',nrOfAngles
              if (iand(rt_healpix_nSide-1,rt_healpix_nSide) /= 0) then
                 write(io,*) "ERROR: rt_healpix_nSide=",rt_healpix_nSide," is not a power of 2."
                 call Driver_abortFlash('Error in raytrace_3DRT, rt_healpix_nSide is not a power of 2')
              endif
              write(io,*) '-------------------------------------------------------------------'
           endif
           !
        elseif(rt_healpix_nSide==-1) then
           !
           ! only one direction
           !
           dirFact(1) = rt_dirX
           dirFact(2) = rt_dirY
           dirFact(3) = rt_dirZ
           if(dr_globalMe.eq.MASTER_PE) then 
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'Ignoring Healpix Tesselation'
              write(io,*) 'because rt_healpix_nSide==0'
              write(io,*) 'Resolution Parameter (rt_healpix_nSide)',rt_healpix_nSide
              write(io,*) 'Using direction vector from the parameter context instead'
              write(io,*) 'dirFact', dirFact
              write(io,*) '-------------------------------------------------------------------'
           endif
           dirFact = dirFact / SQRT(SUM(dirFact**2))
           nrOfAngles = 1
           rt_nrOfAngleGroups = 1
        elseif(rt_healpix_nSide==-2) then
           !
           ! old discretization 
           !
           if(dr_globalMe.eq.MASTER_PE) then
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'Ignoring Healpix Tesselation'
              write(io,*) 'invalid Resolution Parameter (rt_healpix_nSide):',rt_healpix_nSide
              write(io,*) 'rt_healpix_nSide has to be a power of 2 (1,2,4,8,...) '
              write(io,*) 'for the Healpix Tesselation to work'
              write(io,*) 'using rt_nPhi and rt_nTheta for angular discretization the old way'
              write(io,*) 'rt_nPhi',rt_nPhi
              write(io,*) 'rt_nTheta',rt_nTheta
              write(io,*) '-------------------------------------------------------------------'
           endif
           nrOfAngles = rt_nPhi*rt_nTheta
        else
           if(dr_globalMe.eq.MASTER_PE) then
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'No parallel rays at all!'
              write(io,*) 'nrOfAngles = 0.'
              write(io,*) '-------------------------------------------------------------------'
           end if
           nrOfAngles = 0
        endif
     endif
     if(NDIM==2) then
        if(dr_globalMe.eq.MASTER_PE) then
           write(io,*) '-------------------------------------------------------------------'
           write(io,*) 'ATTENTION:'
           write(io,*) 'In the 2D case, we have rt_nPhi valid directions,' 
           write(io,*) 'forcing runtime parameter rt_nTheta=1.'
           write(io,*) 'ignoring rt_healpix_nSide and rt_nTheta in 2D (no healpix)'
           write(io,*) 'using only rt_nPhi:',rt_nPhi
           write(io,*) 'discretization is done in the xy-plane'
           write(io,*) '-------------------------------------------------------------------'
        endif
        rt_nTheta     = 1 
        nrOfAngles = rt_nPhi
     endif
     if(NDIM==1) then
        rt_nTheta     = 1 
        nrOfAngles = rt_nPhi*rt_nTheta
        if(nrOfAngles==2) then
           if(dr_globalMe.eq.MASTER_PE) then
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'ATTENTION:'
              write(io,*) 'In the 1D case, we only have 2 valid directions,' 
              write(io,*) 'The valid directions are along the x- and -x-direction.'
              write(io,*) 'rt_nTheta is ignored and forced to be rt_nTheta=1.'
              write(io,*) '-------------------------------------------------------------------'
           endif
        else
           dirFact = dirFact / SQRT(SUM(dirFact**2))
           if(dr_globalMe.eq.MASTER_PE) then
              write(io,*) '-------------------------------------------------------------------'
              write(io,*) 'ATTENTION:'
              write(io,*) 'In the 1D case, we only have 2 valid directions,' 
              write(io,*) 'Your chosen rt_nPhi is invalid, rt_nPhi:',rt_nPhi
              write(io,*) 'rt_nTheta is ignored and forced to be rt_nTheta=1.'
              write(io,*) 'I assume you only want one direction, so I choose the'
              write(io,*) 'normalized direction from the parameter context, dirFact(1):',dirFact(1)
              write(io,*) '-------------------------------------------------------------------'
           endif
        endif
        
     endif
     !
     if(nrOfAngles>0) then 
        dOmega = 4.d0*PI/real(nrOfAngles)
     else if(nrOfAngles.eq.0) then
        dOmega = 0.
     else 
        if(dr_globalMe.eq.MASTER_PE) then
           write(io,*) 'Something went wrong with nrOfAngles:',nrOfAngles
           write(io,*) 'nrOfAngles must not be negative!'
        endif
        call Driver_abortFlash('Error in rt_init_raytrace_3drt: nrOfAngles is negative')
     endif
     !
     if(rt_nrOfAngleGroups<1) then
        if(dr_globalMe.eq.MASTER_PE) then
           write(io,*) 'WARNING: invalid rt_nrOfAngleGroups:',rt_nrOfAngleGroups
           write(io,*) 'rt_nrOfAngleGroups must be at least 1!'
           write(io,*) 'forcing nrOfAngleGroups=1'
        endif
        rt_nrOfAngleGroups=1
     endif
     !
     if(nrOfAngles.gt.0) then
       nrOfAnglesPerGroup = (nrOfAngles/rt_nrOfAngleGroups)
     else
       ! If the parallel rays are deactivated, still have nrOfAnglesPerGroup != 0,
       ! so we have memory allocated for point source communications.
       nrOfAnglesPerGroup = rt_nrOfAngleGroups
     end if
     !
     if(dr_globalMe.eq.MASTER_PE) then
        open(8,file='solid_angle_grid.out')
        write(8,*) '====================================================='
        write(8,*) 'This file contains the angular discretization used' 
        write(8,*) 'in the raytrace_3DRT subroutine for the calculation'
        write(8,*) 'of the angle averaged mean intensity.'
        write(8,*) 'We print it for your convenience...'
        write(8,*) 'nrOfAngles',nrOfAngles 
        write(8,*) 'dOmega',dOmega
        write(8,*) '====================================================='
        write(8,*) 'iAngle,theta,phi,dirX,dirY,dirZ'
        if(nrOfAngles>1) then
           do iAngle=1,nrOfAngles
              call getAngles(iAngle,0.,0.,0.,theta,phi,dirFact)
              write(8,'(1I3,5F13.8)') iAngle,theta,phi,dirFact
           enddo
        else
           write(8,*) theta,phi,dirFact
        endif
        close(8)
     endif
     !
  end subroutine prepare_solid_angle_grid

end subroutine rt_init_raytrace_3drt
