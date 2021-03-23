module base_grid_interface

implicit none


#include "Flash.h"
#include "constants.h"

#ifndef MPI_INCLUDED
#include "Flash_mpi.h"
#define MPI_INCLUDED
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Grid variable get/set operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Sets the internal energy of a block/grid.
! Note here I assumed you set the density properly first!
FUNCTION set_grid_energy_density(i, j, k, index_of_grid, nproc, enrho, n)

  INTEGER :: n, m, myProc
! cell indices, block index on local proc, proc #
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: en, rho, enrho(n)
  INTEGER :: set_grid_energy_density

  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n
    if (myProc == nproc(m)) then

        !call Grid_putBlkData(index_of_grid, CENTER, EINT_VAR, INTERIOR, [i,j,k],en)
        call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)

        en = enrho(m) / rho

        call Grid_putPointData(index_of_grid(m), CENTER, EINT_VAR, INTERIOR, [i(m),j(m),k(m)], en)

    end if
  end do

  set_grid_energy_density=0
END FUNCTION

! Gets the internal energy of a block/grid.
FUNCTION get_grid_energy_density(i, j, k, index_of_grid, nproc, enrho, n)

  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: en, rho, enrho(n)
  INTEGER :: get_grid_energy_density

  en=0.0
  rho=0.0
  enrho=0.0

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n
    if (myProc == nproc(m)) then

      !call Grid_getBlkData(index_of_grid, CENTER, EINT_VAR, INTERIOR, [i,j,k], en)
      call Grid_getPointData(index_of_grid(m), CENTER, EINT_VAR, INTERIOR, [i(m),j(m),k(m)], en)
      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)

      enrho(m)=en*rho

    end if
  end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, enrho, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  else

    call MPI_Reduce(enrho, enrho, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  end if

  get_grid_energy_density=0
END FUNCTION


FUNCTION get_grid_momentum_density(i, j, k, index_of_grid, nproc, &
                                   rhovx, rhovy, rhovz, n)

  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: rhovx(n), rhovy(n), rhovz(n), vx, vy, vz, rho
  INTEGER :: get_grid_momentum_density

!!! Note in Flash velocities are a  "per mass" variable, unlike most grid codes.
!!! Note this means units of velocity, instead of momentum -MW

  rhovx=0.0; rhovy=0.0; rhovz=0.0
  vx=0.0; vy=0.0; vz=0.0; rho=0.0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, [i(m),j(m),k(m)], vy)
      call Grid_getPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, [i(m),j(m),k(m)], vz)
      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)

      !print*, "vx =", vx, "rho =", rho
      !print*, "vy =", vy, "rho =", rho

      rhovx(m) = rho*vx
      rhovy(m) = rho*vy
      rhovz(m) = rho*vz

    end if

  end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, rhovx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  else

    call MPI_Reduce(rhovx, rhovx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(rhovy, rhovy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(rhovz, rhovz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  end if



  get_grid_momentum_density=0
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Note this function currently sets the velocity field NOT
!!! the momentum density. Flash stores velocity not momentum.
!!! Work in progress!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Modified to assume you set the density correctly already.

FUNCTION set_grid_momentum_density(i, j, k, index_of_grid, nproc, &
                                   rhovx, rhovy, rhovz, n)

  INTEGER :: n, m, myProc
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rhovx, rhovy, rhovz
  DOUBLE PRECISION :: vx, vy, vz, rho
  INTEGER :: set_grid_momentum_density

  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n

  if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho)

      vx = rhovx(m) / rho
      vy = rhovy(m) / rho
      vz = rhovz(m) / rho

      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, [i(m),j(m),k(m)], vx)
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, [i(m),j(m),k(m)], vy)
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, [i(m),j(m),k(m)], vz)

    end if

  end do

  set_grid_momentum_density=0
END FUNCTION

FUNCTION get_grid_velocity(i, j, k, index_of_grid, nproc, &
                           vx, vy, vz, n)

  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: vx, vy, vz
  INTEGER :: get_grid_velocity, myProc, communicator, ierr

  vx = 0.0
  vy = 0.0
  vz = 0.0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vx(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELY_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vy(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELZ_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vz(m))

    end if

  end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, vx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, vy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, vz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  else

    call MPI_Reduce(vx, vx, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(vy, vy, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(vz, vz, n, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
  end if

get_grid_velocity=0
END FUNCTION


FUNCTION set_grid_velocity(i, j, k, index_of_grid, nproc, vx, vy, vz, n)

  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: vx, vy, vz
  INTEGER :: set_grid_velocity, myProc

  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vx(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vy(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, &
                             INTERIOR, [i(m),j(m),k(m)], vz(m))

    end if

  end do

set_grid_velocity=0
END FUNCTION


FUNCTION set_grid_density(i, j, k, index_of_grid, nproc, rho, n)

  INTEGER :: n, m, myProc
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: rho(n)
  INTEGER :: set_grid_density

  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_putPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, [i(m),j(m),k(m)], rho(m))

    end if

  end do

  set_grid_density=0
END FUNCTION


FUNCTION get_grid_density(i, j, k, index_of_grid, nproc, rho, n)

  INTEGER :: get_grid_density, n, m
  INTEGER :: myProc, communicator, ierr
  DOUBLE PRECISION :: rho(n)
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc

  rho=0.0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m),CENTER,DENS_VAR,INTERIOR,[i(m),j(m),k(m)],rho(m))

    end if

  end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, rho, n, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierr)
  else

    call MPI_Reduce(rho, rho, n, MPI_DOUBLE_PRECISION, &
                MPI_SUM, 0, communicator, ierr)
  end if

  get_grid_density=0
END FUNCTION


FUNCTION set_grid_state(i, j, k, index_of_grid, nproc, rho, rhovx, rhovy, rhovz, rhoen, n)

  INTEGER :: n, m, myProc
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rho, rhovx, rhovy, rhovz, rhoen
  DOUBLE PRECISION :: vx, vy, vz, en
  INTEGER :: set_grid_state

  call Driver_getMype(GLOBAL_COMM, myProc)

  do m=1, n

    if (myProc == nproc(m)) then

      vx = rhovx(m) / rho(m); vy = rhovy(m) / rho(m); vz = rhovz(m) / rho(m)
      en = rhoen(m) / rho(m)

      call Grid_putPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], rho(m))
      call Grid_putPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_putPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vy)
      call Grid_putPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vz)
      call Grid_putPointData(index_of_grid(m), CENTER, ENER_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], en)

    end if

  end do

  set_grid_state=0
END FUNCTION

!!! TODO TODO TODO
FUNCTION get_grid_state(i, j, k, index_of_grid, nproc, &
                        rho, rhovx, rhovy, rhovz, rhoen, n)

  INTEGER :: n, m
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION, dimension(n) :: rho, rhovx, rhovy, rhovz, rhoen
  DOUBLE PRECISION :: vx, vy, vz, en
  INTEGER :: get_grid_state, myProc, communicator, ierr

  rho = 0.0
  rhovx=0.0; rhovy=0.0; rhovz=0.0
  rhoen=0.0

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)

  do m=1, n

    if (myProc == nproc(m)) then

      call Grid_getPointData(index_of_grid(m), CENTER, DENS_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], rho(m))
      call Grid_getPointData(index_of_grid(m), CENTER, VELX_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vx)
      call Grid_getPointData(index_of_grid(m), CENTER, VELY_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vy)
      call Grid_getPointData(index_of_grid(m), CENTER, VELZ_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], vz)
      call Grid_getPointData(index_of_grid(m), CENTER, ENER_VAR, INTERIOR, &
                            [i(m),j(m),k(m)], en)

      rhovx(m)=vx*rho(m); rhovy(m)=rho(m)*vy; rhovz(m)=rho(m)*vz
      rhoen(m)=rho(m)*en

    end if
  end do


  if (myProc == 0) then

    !call MPI_Reduce(MPI_IN_PLACE, [rhovx, rhovy, rhovz, rhoen, rho], 5, &
    !                MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovx, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovy, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhovz, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rhoen, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, rho, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)

  else

    !call MPI_Reduce([rhovx, rhovy, rhovz, rhoen, rho], [0.0, 0.0, 0.0, 0.0, 0.0], 5, &
    !                MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(rhovx, rhovx, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(rhovy, rhovy, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(rhovz, rhovz, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(rhoen, rhoen, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)
    call MPI_Reduce(rho, rho, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
        communicator, ierr)

  end if

  get_grid_state=0
END FUNCTION


! Gets the photoelectric flux of a block/grid.
FUNCTION get_grid_flux_photoelectric(i, j, k, index_of_grid, nproc, flux_pe, n)

  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: flux_pe(n)
  INTEGER :: get_grid_flux_photoelectric


call Driver_getComm(GLOBAL_COMM, communicator)
call Driver_getMype(GLOBAL_COMM, myProc)

do m=1, n

  if (myProc == nproc(m)) then

    !call Grid_getBlkData(index_of_grid, CENTER, AFUF_VAR, INTERIOR, [i,j,k], flux_pe)
    call Grid_getPointData(index_of_grid(m), CENTER, AFUF_VAR, INTERIOR, [i(m),j(m),k(m)], flux_pe(m))

  end if

end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, flux_pe, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  else

    call MPI_Reduce(flux_pe, flux_pe, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  end if


  get_grid_flux_photoelectric=0
END FUNCTION

! Gets the photoionizing flux of a block/grid.
FUNCTION get_grid_flux_ionizing(i, j, k, index_of_grid, nproc, flux_ion, n)

  INTEGER :: n, m, myProc, communicator, ierr
  INTEGER, dimension(n) :: i, j, k, index_of_grid, nproc
  DOUBLE PRECISION :: flux_ion(n)
  INTEGER :: get_grid_flux_ionizing


call Driver_getComm(GLOBAL_COMM, communicator)
call Driver_getMype(GLOBAL_COMM, myProc)

do m=1, n

  if (myProc == nproc(m)) then

    !call Grid_getBlkData(index_of_grid, CENTER, AUVF_VAR, INTERIOR, [i,j,k], flux_ion)
    call Grid_getPointData(index_of_grid(m), CENTER, AUVF_VAR, INTERIOR, [i(m),j(m),k(m)], flux_ion(m))

  end if

end do

  if (myProc == 0) then

    call MPI_Reduce(MPI_IN_PLACE, flux_ion, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  else

    call MPI_Reduce(flux_ion, flux_ion, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    0, communicator, ierr)
  end if


  get_grid_flux_ionizing=0
END FUNCTION


FUNCTION get_grid_range(nx, ny, nz, index_of_grid, nproc)
  use Grid_interface, only : Grid_getBlkIndexLimits
  INTEGER :: nx, ny, nz, nproc, myProc, communicator, ierr
  INTEGER :: index_of_grid, blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
  INTEGER :: get_grid_range

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  nx = 0; ny = 0; nz = 0

  if (myProc == nproc) then

      call Grid_getBlkIndexLimits(index_of_grid, blkLimits, blkLimitsGC)
      nx = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) +1  !Note +1 here bc I am
      ny = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) +1 !assuming indexing
      nz = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) +1 !starts at 1.
    !  imin = 1 !(blkLimitsGC(HIGH,IAXIS) - blkLimits(HIGH,IAXIS))/2 !Same here. -Josh
    !  jmin = 1 !(blkLimitsGC(HIGH,JAXIS) - blkLimits(HIGH,JAXIS))/2
    !  kmin = 1 !(blkLimitsGC(HIGH,KAXIS) - blkLimits(HIGH,KAXIS))/2
  end if


  if (myProc == 0) then
      call MPI_Reduce(MPI_IN_PLACE, nx, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(MPI_IN_PLACE, ny, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(MPI_IN_PLACE, nz, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
  else
      call MPI_Reduce(nx, nx, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(ny, ny, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(nz, nz, 1, MPI_INT, MPI_SUM, 0, communicator, ierr)
  end if

  get_grid_range=0
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Notes on get_pos_of_index and get_index_of_pos!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Be aware that if you get the position of an index in a
!!! parent grid of Flash, then pass it back to get the index
!!! of the position, Flash is smart enough to give you the index
!!! of that location in the highest level refinement child grid
!!! instead of giving you the location in the parent grid you
!!! started with. This is a "feature" that took me a while to
!!! figure out. Working as intended! -Josh

FUNCTION get_index_of_position(x, y, z, i, j, k, index_of_grid, Proc_ID, n)
  use Grid_interface, only : Grid_getBlkIndexLimits
  INTEGER, intent(in) :: n
  DOUBLE PRECISION, intent(in), dimension(n) :: x, y, z
  INTEGER, intent(out), dimension(n) :: i, j, k, index_of_grid, Proc_ID
  INTEGER :: get_index_of_position

  DOUBLE PRECISION, DIMENSION(MDIM) :: loc, local_pos, blockCenter, &
                                       blockSize, delta
  INTEGER :: locBlkID, myProc, locProc, communicator
  INTEGER :: ii, nn, ierr

  i = 0  ! we'll take advantage of zero-valued default during MPI_reduce
  j = 0
  k = 0
  index_of_grid = 0
  Proc_ID = 0

  call Driver_getComm(GLOBAL_COMM, communicator)
  call Driver_getMype(GLOBAL_COMM, myProc)

  do nn=1, n

    loc(1) = x(nn)
    loc(2) = y(nn)
    loc(3) = z(nn)
    call Grid_getBlkIDFromPos(loc, locBlkID, locProc, communicator)

    if (myProc == locProc) then

      call Grid_getBlkCenterCoords(locBlkID, blockCenter)
      !call Grid_getBlkCornerID(locBlkID, cornerID, cornerIDMax)
      !call Grid_getBlkIndexLimits(locBlkID, blkLimits, blkLimitsGC)
      call Grid_getBlkPhysicalSize(locBlkID, blockSize)
      call Grid_getDeltas(locBlkID, delta)

      do ii=1, MDIM
        !delta(ii) = blockSize(ii)/(blkLimits(HIGH,ii) - blkLimits(LOW,ii))
        local_pos(ii) = loc(ii) - blockCenter(ii) + blockSize(ii)/2.0
      end do

      ! x,y,z at bottom-left blk faces may give i,j,k=0;
      ! max(...,1) ensures that x,y,z always maps to cells within blk
      ! 
      ! in principle, x,y,z at top-right blk faces may give i,j,k=nxb+1,...
      ! but Grid_getBlkIDFromPos should preclude that scenario, because it
      ! checks "onUpperBoundary"

      i(nn) = max(ceiling(local_pos(1)/delta(1)), 1)

      if (MDIM .gt. 1) then
        j(nn) = max(ceiling(local_pos(2)/delta(2)), 1)
      end if

      if (MDIM .gt. 2) then
        k(nn) = max(ceiling(local_pos(3)/delta(3)), 1)
      end if

      index_of_grid(nn) = locBlkID
      Proc_ID(nn) = locProc

    end if

  end do

  if (myProc == 0) then
    call MPI_Reduce(MPI_IN_PLACE, i, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, j, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, k, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, index_of_grid, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(MPI_IN_PLACE, Proc_ID, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
  else
    call MPI_Reduce(i, i, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(j, j, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(k, k, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(index_of_grid, index_of_grid, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
    call MPI_Reduce(Proc_ID, Proc_ID, n, MPI_INTEGER, MPI_SUM, 0, communicator, ierr)
  end if

  ! clean up for end-user: default values of 0 make MPI_REDUCE easier,
  ! but NONEXISTENT (-1) is a more obvious warning value, and follows
  ! convention from Grid_getBlkIDFromPos(...)
  if (myProc == 0) then  ! after MPI_Reduce, only rank=0 has complete data
    do nn=1,n
      if (index_of_grid(nn) == 0) then  ! valid FLASH blockIDs start from 1
        i(nn) = NONEXISTENT
        j(nn) = NONEXISTENT
        k(nn) = NONEXISTENT
        index_of_grid(nn) = NONEXISTENT
        Proc_ID(nn) = NONEXISTENT
      end if
    end do
  end if

  get_index_of_position=0
END FUNCTION

FUNCTION get_position_of_index(i, j, k, index_of_grid, Proc_ID, x, y, z, n)

  INTEGER, intent(in) :: n
  INTEGER, intent(in), dimension(n) :: i, j, k, index_of_grid, Proc_ID
  DOUBLE PRECISION, intent(out), dimension(n) :: x, y, z
  INTEGER :: get_position_of_index

  INTEGER :: indices(3), ii, myProc, ierr, communicator
  DOUBLE PRECISION :: loc(3)

  call Driver_getMype(GLOBAL_COMM, myProc)
  call Driver_getComm(GLOBAL_COMM, communicator)

  do ii = 1,n
    loc=0.0
    indices(1)=i(ii)
    indices(2)=j(ii)
    indices(3)=k(ii)

    if (myProc == Proc_ID(ii)) then
      call Grid_getSingleCellCoords(indices, index_of_grid(ii), CENTER, INTERIOR, loc)
    end if

    x(ii)=loc(1)
    y(ii)=loc(2)
    z(ii)=loc(3)

  end do

    if (MyProc == 0) then
      call MPI_Reduce(MPI_IN_PLACE, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(MPI_IN_PLACE, y, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(MPI_IN_PLACE, z, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    else
      call MPI_Reduce(x, x, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(y, y, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
      call MPI_Reduce(z, z, n, MPI_DOUBLE_PRECISION, MPI_SUM, 0, communicator, ierr)
    end if

    !x(ii)=loc(1)
    !y(ii)=loc(2)
    !z(ii)=loc(3)
  !end do
  get_position_of_index=0
END FUNCTION

FUNCTION get_leaf_indices(dummy,ind,ret_cnt,num_of_blks,nparts)

  use Driver_data, only : dr_globalNumProcs

  integer :: get_leaf_indices, nparts, num_of_blks, ret_num_blks
  integer :: myProc, ierr, comm, i
  integer, dimension(nparts) :: ind, ret_ind, ret_cnt, dummy
  integer, dimension(dr_globalNumProcs) :: disp, rec_count
  dummy = 0
  disp = 0
  rec_count = 0
  ret_cnt = 0
  ret_ind = -1
  ind = -1
  call Grid_getListOfBlocks(LEAF, ind, num_of_blks)
  call Driver_getComm(GLOBAL_COMM, comm)
  call Driver_getMype(comm, myProc)
  !print*, "ind =", ind, myProc
  !print*, "numblks = ", num_of_blks, myProc
! Gather the array on the root process. Note that we require the
  ! user to pass the proper length of the final array.

  ! Make an array of the # of leaf grids from each processor.
  call MPI_Gather(num_of_blks, 1, MPI_INTEGER, &
                  rec_count, 1, MPI_INTEGER, &
                  0, comm, ierr)
  ret_cnt(:dr_globalNumProcs) = rec_count
  ! Set the displacement for the incoming data based on how many
  ! particles are coming in from each processor. Note the displacement
  ! for the root process is zero, for rank 1 disp = num on root,
  ! for rank 2 disp = num on root + num on 1, etc etc.

  do i=1, dr_globalNumProcs-1

    disp(i+1) = disp(i) + rec_count(i)

  end do
  !print*, "dr_globalNumProcs =", dr_globalNumProcs
  !print*, "rec_count = ", rec_count
  !print*, "disp = ", disp

  ! Now actually gather the leaf grids using the variable length array
  ! gather command in MPI.

  call MPI_Gatherv(ind, num_of_blks, MPI_INTEGER, &
                   ret_ind, rec_count, disp, MPI_INTEGER, &
                   0, comm, ierr)

  ind = ret_ind

  !print*, "ind =", ind, myProc

  !call MPI_Reduce(ret_num_blks,num_of_blks, 1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)

  num_of_blks = sum(rec_count)
  !print*, "numblks =", num_of_blks, myProc

  get_leaf_indices=0
END FUNCTION

FUNCTION get_max_refinement(max_refine)

  use Grid_interface, only : Grid_getMaxRefinement

  INTEGER :: max_refine, get_max_refinement
!  INTEGER, PARAMETER :: mode=4 ! Mode 4 looks at the actual refinement level of the existing blocks.
!  INTEGER, PARAMETER :: scope=4 ! Makes comm=MPI_WORLD_COMM

  call Grid_getMaxRefinement(max_refine)

  get_max_refinement=0
END FUNCTION

end module base_grid_interface
