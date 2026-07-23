!*******************************************************************************

!  Routine:     get_cooling_data

!  Description: reads in the cooling table

        function get_cooling_data(cool_dat, dpts, tpts)

          use Driver_interface, ONLY : Driver_abortFlash

          character(len=512) :: filepath, datadir
          integer get_cooling_data, dpts, tpts
          integer i,j
          real cool_dat(dpts, tpts, 3)

          integer io_status, status
          real temp, dens, cool_pwr, dummy

          get_cooling_data = 0

          ! Get environment variable
          call get_environment_variable("TORCH_DIR", datadir, status=status)
          if (status /= 0 .or. len_trim(datadir) == 0) then
             write (*,*) 'get_cooling_data:  TORCH_DIR env variable not set'
             call Driver_abortFlash("Error: unable to open file cool.dat, TORCH_DIR not set")
             stop  
          end if

          filepath = trim(datadir)//"src/flash/source/physics/sourceTerms/Heat/HeatMain/"// &
                  "HeatCool/phenHeat/mol_and_dust/cool.dat"

          open( unit = 99, file = filepath, status = 'old' &
      &          ,iostat = io_status)


          if (io_status /= 0) then
             write (*,*) 'get_cooling_data:  unable to open file', &
     &                 ' "cool.dat"'
             call Driver_abortFlash("Error: unable to open file cool.dat")
             stop
          endif

          read(99, '(2(2X,I4.4))', iostat = io_status) i, j
          if (io_status /= 0 .OR. i .NE. tpts .OR. j .NE. dpts) then
             write (*,*) 'get_cooling_data:  wrong number of data points in cool.dat', dpts, tpts
             call Driver_abortFlash("Error: wrong number of data points in cool.dat")
          endif

          do j = 1, tpts
            do i = 1,dpts
              read(99, 5010, iostat=io_status) temp, dens, dummy &
     &                  ,dummy, dummy, cool_pwr
              if (io_status .eq. 0) then
                cool_dat(i,j,1) = 10.e0**temp
                cool_dat(i,j,2) = 10.e0**dens
                cool_dat(i,j,3) = 10.e0**cool_pwr
              else
                write (*,*) 'get_cooling_data: read error'
                call Driver_abortFlash("Error: reading cool.dat")
                stop
              endif
            enddo
          enddo

! number of records
          get_cooling_data = (i-1)*(j-1)
5010      FORMAT(6(2X,F8.4))

          return
        end function get_cooling_data

