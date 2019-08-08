!*******************************************************************************

!  Routine:     find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

        function find (x, x0, N)

        integer find, N, il, ir, im
        real    x(N), x0

! default value
        find = 0

        if (x0 .lt. x(1)) then

          find = 0

        elseif (x0 .gt. x(N)) then

          find = N

        else

          il = 1
          ir = N
10          if (ir .eq. il+1) goto 20
            im = (il + ir) / 2
            if (x(im) .gt. x0) then
              ir = im
            else
              il = im
            endif
            goto 10
20        find = il

        endif

        return
        end


