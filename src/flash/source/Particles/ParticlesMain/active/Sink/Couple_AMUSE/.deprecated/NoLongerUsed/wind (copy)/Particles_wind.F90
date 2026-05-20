! Wind injection onto the grid. Winds are injected in pairs and in random
! orientations to avoid grid effects.

! Joshua Wall and Andrew Pellegrino, Drexel University, 05-2016


subroutine Particles_wind(starMass, mass, xloc, yloc, zloc, dt)

!#define debug

use Grid_data, ONLY: gr_meshComm, gr_meshMe

use Grid_interface, ONLY: Grid_getBlkIDFromPos, Grid_getBlkPhysicalSize, &
    Grid_getBlkNeighBlkIDFromPos, Grid_getBlkNeighLevels, Grid_getDeltas, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getMinCellSize, &
    Grid_getBlkIndexLimits, Grid_getCellIndexFromPosition, &
    Grid_notifySolnDataUpdate, Grid_fillGuardCells
    
use Particles_data, ONLY : pt_MaxPerProc

use Eos_interface, ONLY : Eos_wrapped

use RuntimeParameters_interface, ONLY : RuntimeParameters_get

implicit none

#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"
#include "Particles.h"
#include "Eos.h"

real, intent(in) :: starMass, mass, xloc, yloc, zloc
real, intent(in) :: dt ! Return the timestep from the velocity of the gas at max resolution
                       ! May do something with dt later

integer :: nInjectionCells
integer, parameter :: numgroups = 300
integer, parameter :: relNumProp = 5 ! Number of relevant properties
integer, parameter :: injectNumProp=NPART_PROPS ! Number of physical properties of the cloud-in-cell
                           ! injection "particle."

! Injection CIC properties (for each cell): Position x,y,z, mass and total energy 
! Note to use the particle to mesh mapping scheme, the location property
! of our "particles" need to be the same as regular particles.
integer, parameter :: IPOSX=POSX_PART_PROP, IPOSY=POSY_PART_PROP, IPOSZ=POSZ_PART_PROP, &
                      IMASS=MASS_PART_PROP, ITHRE=THRE_PART_PROP, &
                      IVELX=VELX_PART_PROP, IVELY=VELY_PART_PROP, IVELZ=VELZ_PART_PROP, &
                      IBLK=BLK_PART_PROP  , IPROC=PROC_PART_PROP

! Injection cube properties to map to the grid.
integer, dimension(relNumProp) :: partProp = [IMASS, ITHRE, IVELX, IVELY, IVELZ]
!integer, dimension(relNumProp) :: partProp = [IMASS, IVELX, IVELY, IVELZ]
! Variables on the grid to map to.
integer, dimension(relNumProp) :: gridProp = [DENS_VAR, EINT_VAR, VELX_VAR, &
!integer, dimension(relNumProp) :: gridProp = [DENS_VAR, VELX_VAR, &
                                     VELY_VAR, VELZ_VAR] 

real :: pi = 4*atan(1.0)
!real    :: injectBlkSize(3), injectBlkCenter, dVol, 
real    :: MinDelta
real    :: injectMass
!integer :: lb
integer :: injectProc, injectBlk, neighProc, neighBlk
logical :: iHaveInjectBlk, dupBlk, snap_to_grid
logical, save :: first_call=.true.
logical :: use_cube
integer :: numInjectParticles, numBlks
integer, allocatable, dimension(:) :: affectedBlks
integer :: i, j, k, n, mode1, prop, ii, jj, kk, ierr, l
integer :: blkLimits(2,NDIM), blkLimitsGC(2,NDIM)
real    :: background_density, density, cs, vs
real    :: xloc1, yloc1, zloc1
real, save :: gamma

real, allocatable, dimension(:,:) :: rands

real, dimension(:,:,:,:), pointer :: solndata

real, allocatable, dimension(:,:) :: injectGroup

! Rotation matrices
real :: zrot(3,3), i3mat(3,3), rotationMatrix(3,3), rv(3), vvtranspose(3,3)
! Rotation angles
real :: axisRot, axisPhi, axisTheta

real :: vterm, rStar
real :: injectRadius

! Note cell sizes are the same as the deltas for the
! highest level of refinement.

! We will probably have to reshape this to be (num_props, num_cells)
! to be more like the particles array.

! Reshape into this array.
real, allocatable, dimension(:,:) :: particles_injection

! Dummy particle for mapping from grid to injection cell.
real, allocatable, dimension(:,:) :: dummy_particle

! Attribute array for mapping from the grid to the injection "particles."
integer :: attrib(2,relNumProp)

#ifdef debug
if (gr_meshMe == MASTER_PE) then
        print *, "-------------------- Injecting wind ---------------------"
end if
#endif

snap_to_grid = .false. !.true.
iHaveInjectBlk = .false.
density = 0.0
! If true, inject winds in a 3x3x3 cube on top of the grid around the star.
! If false, inject winds in random directions at an arbitrary distance
use_cube = .true.
nInjectionCells = 2 !should be even if use_cube is false


if (use_cube) then
    nInjectionCells = 27
else
    !make sure it's even to insert pairs
    nInjectionCells = nInjectionCells + mod(nInjectionCells,2)
end if

#ifdef debug
if (gr_meshMe == MASTER_PE) then
        print*, "nInjectionCells = ", nInjectionCells
        print*, "numgroups = ", numgroups
end if
#endif

allocate (affectedBlks( nInjectionCells*numgroups ))
allocate (injectGroup( nInjectionCells*numgroups, injectNumProp ))
allocate (particles_injection( NPART_PROPS, nInjectionCells*numgroups ))
allocate (dummy_particle( NPART_PROPS, nInjectionCells*numgroups ))

if (first_call) then
  call RuntimeParameters_get('gamma', gamma)
#ifdef debug  
  print*, "Gamma equals", gamma
#endif
  first_call = .false.
end if

#ifdef debug
if (gr_meshMe == MASTER_PE) then
  write(*,'(A,3ES12.3E3,A,L)') "[Particles_wind]: About to inject winds at", &
           xloc, yloc, zloc, " and snap to grid is ", snap_to_grid
end if
#endif

call Grid_getMinCellSize(MinDelta)

#ifdef debug
  write(*,'(A3ES12.3)') "Min Delta =", MinDelta
#endif

! Snap to grid?
if (snap_to_grid) then
  xloc1 = xloc + 0.5*MinDelta
  yloc1 = yloc + 0.5*MinDelta
  zloc1 = zloc + 0.5*MinDelta
else
  xloc1 = xloc
  yloc1 = yloc
  zloc1 = zloc
end if

injectGroup(:,:) = 0.0
dummy_particle(:,:) = 0.0
attrib(1,:) = partProp
attrib(2,:) = gridProp


! Fill in the injection blocks
! Maybe also try to map the heating to the PHHE variable on the grid...

if (use_cube) then
    allocate (rands( 3, numgroups ))

    !Place 27 imaginary cells in [numgroups] cubes
    do j=1, numgroups
        do i=1, nInjectionCells
            injectGroup((j-1)*nInjectionCells+i,IPOSX) = MinDelta * (mod( (i-1),   3) - 1)
            injectGroup((j-1)*nInjectionCells+i,IPOSY) = MinDelta * (mod( (i-1)/3, 3) - 1)
            injectGroup((j-1)*nInjectionCells+i,IPOSZ) = MinDelta * (     (i-1)/9     - 1)
        end do

        ! Choose random orientation
        ! Method from Arvo 1995
        call random_number(rands)
        axisRot = rands(1,j)
        axisTheta = rands(2,j) * 2 * pi
        axisPhi = rands(3,j) * 2 * pi
        ! reflection vector
        rv = (/ cos(axisPhi)*sqrt(axisRot), sin(axisPhi)* &
                    sqrt(axisRot), sqrt(1-axisRot) /)
        ! reflection vector times its transpose
        vvtranspose(:,1) = rv(1) * rv
        vvtranspose(:,2) = rv(2) * rv
        vvtranspose(:,3) = rv(3) * rv

        !identity matrix of size 3. I can't believe fortran doesn't have a builtin
        !for this
        i3mat(1,:) = (/ 1, 0, 0 /) 
        i3mat(2,:) = (/ 0, 1, 0 /) 
        i3mat(3,:) = (/ 0, 0, 1 /) 
           
        !zrot = transpose(reshape( &
        !            (/  cos(axisTheta), sin(axisTheta), 0.0, &
        !               -sin(axisTheta), cos(axisTheta), 0.0, &
        !                         0.0,          0.0, 1.0 /), &
        !        shape(zrot)))

        ! rotation matrix around the z axis
        zrot(1,:) = (/  cos(axisTheta), sin(axisTheta), 0.0 /)
        zrot(2,:) = (/ -sin(axisTheta), cos(axisTheta), 0.0 /)
        zrot(3,:) = (/             0.0,            0.0, 1.0 /)

        rotationMatrix = matmul(2*vvtranspose-i3mat,zrot)

        do i=1, nInjectionCells
            injectGroup((j-1)*nInjectionCells+i,IPOSX:IPOSZ) = &
            matmul(injectGroup((j-1)*nInjectionCells+i,IPOSX:IPOSZ),rotationMatrix)
            !matmul(rotationMatrix,injectGroup((j-1)*nInjectionCells+i,IPOSX:IPOSZ))
        end do
    end do
else
    allocate (rands( nInjectionCells, 2 ))
    !injectRadius = MinDelta * sqrt(3.0)
    injectRadius = MinDelta

    ! Place [nInjectionCells] cells at distance R randomly
    call random_number(rands)
    !injectRadius = MinDelta * sqrt(3.0)
    injectRadius = MinDelta
    do i=1, nInjectionCells, 2
        ! choose random spherical coordinates
        ! axisTheta and axisPhi describe a random point on a sphere in spherical
        ! polar coordinates with a uniform distribution. axisRot defines the
        ! rotation in 3 dimensions around the axis defined by that point on a
        ! sphere and the center.

        axisTheta = acos(2 * rands(i,1) - 1)
        axisPhi = rands(i,2) * 2 * pi

        ! give cell position at distance injectRadius at these coordinates
        injectGroup(i,IPOSX) = injectRadius * sin(axisTheta) * cos(axisPhi)
        injectGroup(i,IPOSY) = injectRadius * sin(axisTheta) * sin(axisPhi)
        injectGroup(i,IPOSZ) = injectRadius * cos(axisTheta)

        ! place next cell in opposing position to conserve momentum
        injectGroup(i+1,IPOSX) = -injectGroup(i,IPOSX)
        injectGroup(i+1,IPOSY) = -injectGroup(i,IPOSY)
        injectGroup(i+1,IPOSZ) = -injectGroup(i,IPOSZ)

    
        if (gr_meshMe == MASTER_PE) then
            print*, "phi = ", axisPhi
            print*, "theta = ", axisTheta
            print*, "X = ", injectGroup(i,IPOSX), injectGroup(i+1,IPOSX)
            print*, "Y = ", injectGroup(i,IPOSY), injectGroup(i+1,IPOSY)
            print*, "Z = ", injectGroup(i,IPOSZ), injectGroup(i+1,IPOSZ)
        end if
    end do

end if

! velocity of each block

rStar = 4e11 ! cm, R ~ M**0.5, M ~ 30 solar masses
! Fit of wind velocity at r = inf vs. star mass, from Dale et al. 2013
vterm = 1.019430 * ( starMass - 3.579183e34 )**0.24 + 6e7

! rStar << injectRadius in most cases, vs ~ vterm
! vs = vterm * sqrt(1-rStar/injectRadius)

do i=1, nInjectionCells*numgroups
    ! Wind speed solution from Stahler, Palla p. 540

    ! Divide total velocity into cartesian directions s.t. velocity is 
    ! radially outwards
    if (abs(injectGroup(i,IPOSX)) < epsilon(injectGroup(i,IPOSX)) .and. &
        abs(injectGroup(i,IPOSY)) < epsilon(injectGroup(i,IPOSY)) .and. &
        abs(injectGroup(i,IPOSZ)) < epsilon(injectGroup(i,IPOSZ)) ) then
        ! center cell. v = 0
        injectGroup(i,IVELX:IVELZ) = 0.0
    else
        injectGroup(i,IVELX:IVELZ) = vterm * injectGroup(i,IPOSX:IPOSZ) / &
                        sqrt(sum(injectGroup(i,IPOSX:IPOSZ)*injectGroup(i,IPOSX:IPOSZ)))
    end if
end do

! mass of each block
injectMass = mass / (nInjectionCells * numgroups)

do i=1, nInjectionCells*numgroups
        injectGroup(i,IMASS) = injectMass ! WARNING! The particle->grid map divides by volume.
end do

#ifdef debug
print *, "------ injectGroup ------"
print *, "mindelta ", MinDelta
print *, "loc ", xloc1, yloc1, xloc1
print*, "Injection positions are :"
do i=1, nInjectionCells*numgroups
     print*, injectGroup(i,IPOSX:IPOSZ)
end do

print*, "Injection velocities are :"
do i=1, nInjectionCells*numgroups
    print*, injectGroup(i,IVELX:IVELZ)
end do

print*, "Injection masses are :"
do i=1, nInjectionCells*numgroups
    print*, injectGroup(i,IMASS)
end do

#endif
! Now center the cube at position xloc1, yloc1, zloc1.
! Note this must be done before we map the grid up to the
! dummy particle. - JW

injectGroup(:,IPOSX) = injectGroup(:,IPOSX) + xloc1
injectGroup(:,IPOSY) = injectGroup(:,IPOSY) + yloc1
injectGroup(:,IPOSZ) = injectGroup(:,IPOSZ) + zloc1

! Set up dummy particle array

do i=1, nInjectionCells*numgroups 
    ! Map the current grid values to a dummy particle at this location.
    call Grid_getBlkIDFromPos([injectGroup(i,IPOSX),injectGroup(i,IPOSY), &
                                   injectGroup(i,IPOSZ)], injectBlk, injectProc, &
                                   gr_meshComm)

!         if (gr_meshMe .eq. injectProc) then
            dummy_particle(IBLK, i) = injectBlk
            dummy_particle(IPROC,i) = injectProc
            dummy_particle(IPOSX,i) = injectGroup(i,IPOSX)
            dummy_particle(IPOSY,i) = injectGroup(i,IPOSY)
            dummy_particle(IPOSZ,i) = injectGroup(i,IPOSZ)

#ifdef debug
        write(*,'(A,3ES12.3)') "Dummy part blk  =", dummy_particle(IBLK, i)
        write(*,'(A,3ES12.3)') "Dummy part prc  =", dummy_particle(IPROC,i)
        write(*,'(A,3ES12.3)') "Dummy part posx =", dummy_particle(IPOSX,i)
        write(*,'(A,3ES12.3)') "Dummy part posy =", dummy_particle(IPOSY,i)
        write(*,'(A,3ES12.3)') "Dummy part posz =", dummy_particle(IPOSZ,i)
#endif
end do

! Each dummy particle must be mapped separately, because there is no
! way to know in advance which processor the block containing the dummy
! particle will be on. If you try and map from the wrong processor,
! it will act like it worked except it will use a parent block and likely
! the map will occur outside the domain. You'll get errors about hitting
! an index in solnvec outside its allow index range. - JW

do i=1, nInjectionCells*numgroups
    if (dummy_particle(IPROC,i) == gr_meshMe) then
        call Grid_mapMeshToParticles(dummy_particle(:,i), NPART_PROPS, IBLK, &
                                     1, [IPOSX, IPOSY, IPOSZ], relNumProp, attrib, &
                                     WEIGHTED, CENTER)
    end if
end do

#ifdef debug
  print*, "Dummy part mass =", dummy_particle(IMASS,:)
  print*, "Dummy part vel  =", dummy_particle(IVELX:IVELZ,:)
#endif

do i=1, relNumProp ! Loop over the dummy particle properties, reducing each in place.
    call MPI_ALLREDUCE(MPI_IN_PLACE, dummy_particle(partProp(i),:), nInjectionCells*numgroups, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, gr_meshComm, ierr)
end do

! Now average the background density by dividing by the # of injection cells
background_density = sum(dummy_particle(IMASS,:)) / (nInjectionCells * numgroups)
#ifdef debug
write(*,'(A,ES12.3E3)') "Background density is", background_density
#endif


! Calculate the timestep for this gas using method in Simpson et al. 2015
density = background_density + injectMass/(MinDelta**3.0)
! dt = 0.1 * MinDelta / vs ! / (cs + vs)

#ifdef debug
if (gr_meshME == MASTER_PE) then
  write(*,'(A,ES12.3E3)') "[Particles_wind]: dt is =", dt
  write(*,'(A,ES12.3E3)') "[Particles_wind]: inject mass per cell is =", injectMass
  write(*,'(A,ES12.3E3)') "[Particles_wind]: inject vel per cell is =", vterm
end if
#endif

! Put the injection cube into an array that the mapping functions
! will understand (something that looks like the Flash particles array).

! Also check to see if this cell is local on this processor, if not 
! leave that cell filled with zeroes on this processor.

particles_injection(:,:) = 0.0
numBlks = 0
affectedBlks = 0
numInjectParticles = 0
! New stuff, should actually conserve momentum when adding the injection stuff.

do i=1, nInjectionCells*numgroups
      call Grid_getBlkIDFromPos([injectGroup(i,IPOSX),injectGroup(i,IPOSY),injectGroup(i,IPOSZ)], &
                                injectBlk, injectProc, gr_meshComm)
      ! If you're on the local proc, great, copy over the values to the array.
      ! We also should check here for max refinement!
      if (gr_meshMe .eq. injectProc) then
        iHaveInjectBlk = .true.
        numInjectParticles   = numInjectParticles + 1
#ifdef debug         
        print*, "Injection block and proc are:", injectBlk, injectProc
#endif

      ! Store affected blocks.
      ! Check if this block is already stored.
      dupBlk = .false.
      do n=1, numBlks
        if(affectedBlks(n) == injectBlk) then
          dupBlk = .true.
          exit
        endif
      end do
      ! If not, store it now.
      if (dupBlk .eqv. .false.) then
        numBlks = numBlks + 1
        affectedBlks(numBlks) = injectBlk
      end if

      ! Map the current grid values to a dummy particle at this location.  
        !dummy_particle(:,1) = 0.0
        !dummy_particle(IBLK,1)  = injectBlk
        !dummy_particle(IPROC,1) = injectProc
        !dummy_particle(IPOSX:IPOSZ,1) = injectGroup(i,IPOSX:IPOSZ)

        !call Grid_mapMeshToParticles(dummy_particle, NPART_PROPS, IBLK, &
        !                    1, [IPOSX, IPOSY, IPOSZ], 4, attrib, &
        !                    WEIGHTED, CENTER)
#ifdef debug
      ! print*, "Dummy part mass =", dummy_particle(IMASS,1)
      ! print*, "Dummy part vel  =", dummy_particle(IVELX:IVELZ,1)
#endif

      ! Calculate the average background density for the timestep calculation.
      ! This will have to be MPI summed for parallel calculations!
      ! background_density = background_density + dummy_particle(IMASS,1)
       
      ! Now use this to calculate the new velocity in the injection cell, where
      ! we actually use momentum conservation this time. Yay.
      ! Note that we subtract the cell values for velocity because when we map back
      ! we actually add to whatever is in the cell currently. This
      ! lets us map zeros safely for cells that are off proc.
    
        !!! Note here!!!!
        ! When mapping particle quantities to the grid, the routine
        ! itself divides the quantity by the volume of the cell to
        ! make it a "density" type thing. This means we just map mass,
        ! and we map velocity times the cell volume so that the proper
        ! velocity gets mapped to the cell. We also multiply thermal
        ! energy by the cell volume and divide by mass, since eint is
        ! specific internal energy. Note this does not happen when
        ! mapping from the GRID to the PARTICLE (for the dummy).

        particles_injection( IPOSX:IPOSZ , numInjectParticles ) &
              = injectGroup(i,IPOSX:IPOSZ)

        particles_injection( IVELX:IVELZ , numInjectParticles ) & 
              = ( ( injectMass * injectGroup(i,IVELX:IVELZ) &
              + (MinDelta**3.0) * dummy_particle(IMASS,i) / numgroups &
              * dummy_particle(IVELX:IVELZ,i) ) &
              / (numgroups * injectMass + (MinDelta**3.0) &
              * dummy_particle(IMASS,i) ) &
              - dummy_particle(IVELX:IVELZ, i) ) * MinDelta**3.0 / numgroups

!        particles_injection( IVELX:IVELZ , numInjectParticles ) & 
!              = injectMass / ( dummy_particle(IMASS,i) + injectMass )  &
!              * (injectGroup(i,IVELX:IVELZ) - dummy_particle(IVELX:IVELZ,i)) &
!              * MinDelta**3.0

        particles_injection(ITHRE, numInjectParticles ) &
              = injectGroup(i, ITHRE)

        particles_injection(IMASS, numInjectParticles ) &
              = injectGroup(i,IMASS)

        particles_injection(IBLK, numInjectParticles )  = injectBlk

        particles_injection(IPROC, numInjectParticles ) = injectProc

      ! Now we only add particle_injection cells if that location is local to
      ! this processor. This fixes the broken mapping when mapping these
      ! injection particles back to the grid, since otherwise this processor
      ! would attempt to map a particle outside its grid domain. As above,
      ! this results in getting a parent block and finally with solnvec
      ! complaining that your calling outside its index range. - JW 

      end if
end do

#ifdef debug
print *, "------ particles_injection ------"
print*, "Injection positions are :"
do i=1, nInjectionCells*numgroups
  print*, particles_injection(IPOSX:IPOSZ,i)
end do
print*, "Injection velocities are :"
do i=1, nInjectionCells*numgroups
  print*, particles_injection(IVELX:IVELZ,i)
end do
print*, "Injection masses are :"
do i=1, nInjectionCells*numgroups
  print*, particles_injection(IMASS,i)
end do

if (numInjectParticles==0) print*, "I'm proc ", gr_MeshMe, " and I've not no &
                                   injection particles to speak of."
#endif

! Now all we have to do is map each injection cube cell to the grid
! using the built in functions that map particles to the grid.

! Note that all processors have to participate in this call, even the
! ones with no particles. This is because they might border a block that
! does have particles and therefore gets involved in MPI communications
! with neighbor processors.

!if (iHaveInjectBlk) then

! Might need to loop over the cube properties here, as well as over i and j.

mode1 = 1 ! Don't zero the grid variables, just add the particle_injection
          ! amounts to them.

do prop=1, relNumProp
#ifdef debug  
    print*, "I'm proc", gr_meshMe, "mapping part prop", partProp(prop), "to grid prop", gridProp(prop)
#endif

    call Grid_mapParticlesToMesh(particles_injection, NPART_PROPS, &
                               numInjectParticles, pt_maxPerProc, partProp(prop), gridProp(prop), &
                                 mode1)
end do

! Update the EOS variables if we have an affected block for only those blocks.

!  call Grid_notifySolnDataUpdate()
!end if

if (iHaveInjectBlk) then
  do n=1, numBlks
#ifdef debug
    print*, "[Particles_wind]: Calling Eos_wrapped on blk, proc ", affectedBlks(n), gr_meshMe
#endif
    call Grid_getBlkIndexLimits(affectedBlks(n),blkLimits,blkLimitsGC)
    call Eos_wrapped(MODE_DENS_EI, blkLimits, affectedBlks(n))
  end do
!end if
else
#ifdef debug
   print*, "I'm proc", gr_meshMe, "and I have no injection blocks."
#endif
!  return
end if

call MPI_Barrier(gr_meshComm, ierr)
!print*,"I'm proc", gr_meshMe, "between barriers."
call Grid_fillGuardCells(CENTER, ALLDIR) !, eosMode=MODE_DENS_EI, doEos=.true.)

!call MPI_Barrier(gr_meshComm, ierr)
#ifdef debug
print*, "I'm proc", gr_meshMe, "returning now."
#endif
return

end subroutine Particles_wind
