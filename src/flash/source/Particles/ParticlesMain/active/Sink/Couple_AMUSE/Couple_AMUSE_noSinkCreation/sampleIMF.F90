module sampleIMF

contains

subroutine sample_IMF(total_mass, min_mass, stars, n_stars, initialize_seed)

implicit none

real, intent(INOUT)                          :: total_mass
real, intent(IN)                             :: min_mass
real, dimension(:), intent(INOUT)            :: stars
integer, intent(out)                         :: n_stars
logical, intent(in), optional                :: initialize_seed

!real, dimension(:), allocatable :: stars_temp
real :: a, rand_mass, prob, rand_num, initial_mass
real, parameter :: alpha=-2.35
integer :: arr_len

logical, save :: first_call=.true.
logical  :: allow_oversampling

allow_oversampling = .false.

if (present(initialize_seed) .and. first_call) &
  first_call = initialize_seed
if (first_call) then
  call init_random_seed()
  first_call = .false.
end if

! Note that the average mass of a star should be about twice the minimum mass.
! We use this to estimate the length of the array we need to hold the returned stars.

arr_len = ceiling(total_mass/(2.0*min_mass)*1.2) ! Assume average mass + 20%

!allocate(stars_temp(arr_len))

initial_mass = total_mass

do while (total_mass .gt. min_mass)

  !a = norm_factor(total_mass, min_mass, alpha)
  rand_mass = 0.0
  
  do while (rand_mass < min_mass)
    
    call random_number(rand_mass)
    
    ! Use the initial mass every time. This means that you generally
    ! will oversample by a small amount, and very rarely might oversample
    ! by a lot. But it keeps you from pulling too much from the low mass
    ! end of the IMF. We just "borrow" the extra mass from the surrounding
    ! gas.
    
    !rand_mass = rand_mass*total_mass
    rand_mass = rand_mass*initial_mass
    
  end do
  
  !prob = a * Salpeter(rand_mass)
  prob = Salpeter(rand_mass)
  call random_number(rand_num)
  
  if (prob >= rand_num) then
    n_stars = n_stars + 1
    stars(n_stars) = rand_mass
    !total_mass = total_mass - rand_mass
    
    if (allow_oversampling) then
    
        if (sum(stars(1:n_stars)) .ge. 1.1*initial_mass) &
            stars(n_stars) = 1.1*initial_mass - sum(stars(1:n_stars-1))
    
    else
    
        if (sum(stars(1:n_stars)) .ge. initial_mass) &
            stars(n_stars) = initial_mass - sum(stars(1:n_stars-1))
            
    end if
    
    total_mass = initial_mass - sum(stars(1:n_stars))
    !print*, "[sampleIMF]: total mass =", total_mass
  end if
  
end do

!allocate(stars(n_stars))
!stars(1:n_stars) = stars_temp(1:n_stars)

!deallocate(stars_temp)

end subroutine sample_IMF

subroutine IMF_prob(total_mass, min_mass, rand_mass, prob, initialize_seed)

implicit none

real, intent(IN)                             :: total_mass
real, intent(IN)                             :: min_mass
real, intent(out)                            :: prob, rand_mass
logical, intent(in), optional                :: initialize_seed

!real, dimension(:), allocatable :: stars_temp
real :: a, rand_num
real, parameter :: alpha=-2.35

logical, save :: first_call=.true.

if (present(initialize_seed) .and. first_call) &
  first_call = initialize_seed
if (first_call) then
  call init_random_seed()
  first_call = .false.
end if

a = norm_factor(total_mass, min_mass, alpha)

call random_number(rand_num)
rand_mass = rand_num*(total_mass - min_mass) + min_mass

prob = Salpeter(rand_mass)
return
end subroutine IMF_prob


real function Salpeter(mass)
implicit none
real, intent(in) :: mass

Salpeter = mass**(-2.35)

end function Salpeter


real function norm_factor(m_tot, m_min, alpha)
implicit none
real, intent(IN) :: m_tot, m_min, alpha
real             :: integral

integral = 1.0/(alpha + 1.0) * (m_tot**(alpha + 1.0) - m_min**(alpha + 1.0))
norm_factor = 1.0 / integral

end function norm_factor


! NOTE: This subroutine obtained from the g fortran man pages.
!SUBROUTINE init_random_seed() 
!IMPLICIT NONE
!INTEGER :: i, n, clock
!INTEGER, DIMENSION(:), ALLOCATABLE :: seed

!CALL RANDOM_SEED(size = n)
!ALLOCATE(seed(n))

!CALL SYSTEM_CLOCK(COUNT=clock)

!seed = clock + 37 * (/ (i - 1, i = 1, n) /)
!CALL RANDOM_SEED(PUT = seed)

!DEALLOCATE(seed)
!END SUBROUTINE

  subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t
  
    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed

end module sampleIMF
