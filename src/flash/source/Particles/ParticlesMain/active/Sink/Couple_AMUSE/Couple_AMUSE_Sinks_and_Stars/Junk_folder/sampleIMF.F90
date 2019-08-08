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

! Stubby

return

end subroutine

end module sampleIMF
