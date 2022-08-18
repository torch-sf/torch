module cool_vars

use OdeData, only : dp
! Added for doing RK5 integration of heating and cooling so I
! don't have to pass them to the d (ei) / dt function. - JW

real(dp), save     :: ndens, rho, conf, Gflux, ephen, &
                      tdust, TtoEI, xHp, strat_factor, mu_mol

!logical, save      ::  he_int_method

end module cool_vars
