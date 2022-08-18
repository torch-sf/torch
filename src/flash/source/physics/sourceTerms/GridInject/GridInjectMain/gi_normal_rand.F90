! https://www.jstor.org/stable/2027592
! A Convenient Method for Generating Normal Variables
! G. Marsaglia and T. A. Bray
! SIAM Review
! Vol. 6, No. 3 (Jul., 1964), pp. 260-264

! This implementation from
! http://www.techchug.com/articles/tutorials/2017/09/Generating-Random-Numbers-With-a-Normal-Distribution-in-Fortran

subroutine gi_normal_rand(mean, std_dev, randnum)

  implicit none

  real, intent(IN) :: mean, std_dev
  real, intent(OUT) :: randnum
  real :: x, y, r
  real, save :: spare
  logical, save :: has_spare
  logical, save :: first_call = .true.

  if (first_call) then
    call random_seed()
    first_call = .false.
  end if

  ! use a spare saved from a previous run if one exists
  if (has_spare) then
    has_spare = .FALSE.
    randnum = mean + (std_dev * spare)
    return
  else
    r = 1.0
    do while ( r >= 1.0 )
      ! generate random number pair between 0 and 1
      call random_number(x)
      call random_number(y)
      ! normalise random numbers to be in square of side-length = R
      x = (x * 2.0) - 1.0
      y = (y * 2.0) - 1.0
      r = x*x + y*y
    end do

    ! calculate the co-efficient to multiply random numbers x and y
    ! by to achieve normal distribution
    r = sqrt((-2.0 * log(r)) / r)

    randnum = mean + (std_dev * x * r)
    spare = y * r
    has_spare = .TRUE.
    return
  end if
end subroutine gi_normal_rand
