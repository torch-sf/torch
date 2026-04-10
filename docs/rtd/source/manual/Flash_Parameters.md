# Flash Parameters

## Simulation
| Variable (=default)           | Type | Description |
|:------------------------------|:-----|:------------|
| basenm            = "turbsph_"           |      |             |
| log_file          = "turbsph.log"        |      |             |
| run_comment       = "Turbulent sphere"   |      |             |
| output_directory  = "./data"             |      |             |

## Simulation/SimulationMain/Cube
| Variable (=default)           | Type | Description |
|:------------------------------|:-----|:------------|
| sim_cubeFile      = "./cube128"          |      |             |
| sim_init_Hp       = 1e-4                 |      |             |
| sim_tdust         = 18.0                 |      |             |
| sim_A_n           =  1.2784              |      |             |
| sim_A_i           =  0.6698              |      |             |
| bx0               = 0.0e0                |      |             |
| by0               = 0.0e0                |      |             |
| bz0               = 3.0e-6               |      |             |

## physics/Eos/EosMain
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| gamma             = 1.66666666666666667  |      |             |

## Grid/GridMain
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| xmax            =  2.1602e+19             |      |             |
| xmin            = -2.1602e+19             |      |             |
| ymax            =  2.1602e+19             |      |             |
| ymin            = -2.1602e+19             |      |             |
| zmax            =  2.1602e+19             |      |             |
| zmin            = -2.1602e+19             |      |             |

| Variable           | Type | Description |
|:-------------------|:-----|:------------|
| xl_boundary_type      = "outflow" |      |             |
| xr_boundary_type      = "outflow" |      |             |
| yl_boundary_type      = "outflow" |      |             |
| yr_boundary_type      = "outflow" |      |             |
| zl_boundary_type      = "outflow" |      |             |
| zr_boundary_type      = "outflow" |      |             |

## Grid/GridMain/paramesh/paramesh4
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| nblockx         = 1                     |      |             |
| nblocky         = 1                     |      |             |
| nblockz         = 1                     |      |             |
| lrefine_max     = 2                     |      |             |
| lrefine_min     = 2                     |      |             |
| refine_var_1    = "pres"                 |      |             |
| refine_filter_1 = 1e-2                   |      |             |
| refine_cutoff_1 = 0.98                   |      |             |
| derefine_cutoff_1 = 0.6                  |      |             |
| gr_sanitizeDataMode  = 3                  |      |             |
| gr_sanitizeVerbosity = 4                   |      |             |
| gr_pmrpConserve          = .true.          |      |             |

## Driver/DriverMain, IO/IOMain
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| restart         = .false.                  |      |             |
| nend            = 999999999                |      |             |
| tmax            = 6.3113852e13             |      |             |
| dtinit          = 3.15e8                    |      |             |
| dtmin           = 3.15e7                    |      |             |
| dtmax           = 3.15e12                   |      |             |
| dr_shortenLastStepBeforeTMax = .true.       |      |             |
| wall_clock_checkpoint       = 0.432e+05         |      |             |
| wall_clock_time_limit       = 1.0e99            |      |             |
| checkpointFileNumber        = 0                 |      |             |
| checkpointFileIntervalStep  = 0                 |      |             |
| checkpointFileIntervalTime  = 3.1556926e13      |      |             |
| plotFileNumber              = 0                 |      |             |
| plotFileIntervalStep        = 0                 |      |             |
| plotFileIntervalTime        = 3.1556926e12      |      |             |
| particleFileNumber          = 0                 |      |             |
| particleFileIntervalStep    = 0                 |      |             |
| particleFileIntervalTime    = 3.1556926e12      |      |             |
| dr_usePosdefComputeDt   = .true.           |      |             |
| dr_numPosdefVars        = 4                  |      |             |
| dr_posdefDtFactor       = 0.9                 |      |             |
| dr_posdefVar_1          = "ener"               |      |             |
| dr_posdefVar_2          = "eint"               |      |             |
| dr_posdefVar_3          = "dens"               |      |             |
| dr_posdefVar_4          = "pres"               |      |             |
| plot_var_1 = "dens"                           |      |             |
| plot_var_2 = "pres"                           |      |             |
| plot_var_3 = "temp"                           |      |             |
| plot_var_4 = "velx"                           |      |             |
| plot_var_5 = "vely"                           |      |             |
| plot_var_6 = "velz"                           |      |             |
| plot_var_7 = "ihp"                            |      |             |
| plot_var_8 = "iha"                            |      |             |
| plot_var_9 = "eint"                           |      |             |
| plot_var_10 = "ener"                          |      |             |
| plot_var_11 = "tdus"                          |      |             |
| plot_var_12 = "magx"                          |      |             |
| plot_var_13 = "magy"                          |      |             |
| plot_var_14 = "magz"                          |      |             |
| plot_var_15 = "uvfl"                          |      |             |
| plot_var_16 = "fufl"                          |      |             |
| plot_var_17 = "magp"                          |      |             |
| plot_var_18 = "gpot"                          |      |             |
| plot_var_19 = "bgpt"                          |      |             |

## Hydro solver
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| useHydro       = T                      |      |             |
| cfl		= 0.8                           |      |             |
| eintSwitch      = 1.e-4                  |      |             |
| killdivb    = T                          |      |             |
| UnitSystem  = "CGS"                      |      |             |
| smallt = 1.0d-99                        |      |             |
| smallp = 1.3807e-19                     |      |             |
| smallu = 1.0d-99                        |      |             |
| smallx = 1.0d-99                        |      |             |
| smalle = 1.0d-99                        |      |             |
| smlrho = 1.67e-28                       |      |             |

## split/ppm solver
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| hydrocomputedtoption = -1                  |      |             |
| order		    = 3                           |      |             |
| slopeLimiter    = "minmod"                  |      |             |
| LimitedSlopeBeta= 1.                          |      |             |
| charLimiting	= .true.                       |      |             |
| use_avisc	    = .true.                       |      |             |
| cvisc		    = 0.1                          |      |             |
| use_flattening	= .false.                      |      |             |
| use_steepening	= .false.                      |      |             |
| use_upwindTVD	= .true.                       |      |             |
| use_gravhalfupdate = .true.                  |      |             |
| use_hybridOrder     = .true.                  |      |             |
| hy_fallbackLowerCFL = .true.                  |      |             |
| hy_cflFallbackFactor = 0.9                     |      |             |
| E_modification	= .true.                      |      |             |
| E_upwind        = .true.                      |      |             |
| energyFix	    = .true.                      |      |             |
| ForceHydroLimit	= .false.                     |      |             |
| prolMethod      = "injection_prol"            |      |             |

##	RIEMANN SOLVERS
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| RiemannSolver	= "HLLD"      |      |             |
| entropy         = .false.     |      |             |
| EOSforRiemann	= .false.     |      |             |

##	STRONG SHOCK HANDELING SCHEME:
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| shockDetect	    = .true.      |      |             |
| shockLowerCFL   = .true.       |      |             |
| useSTS                  = .false. |      |             |
| nstepTotalSTS           = 5        |      |             |
| nuSTS                   = 0.2      |      |             |

## Particles/ParticlesMain
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| useParticles        = T                |      |             |
| useSinkParticles    = .true.           |      |             |
| pt_maxPerProc   = 10000                |      |             |

## Particles/ParticlesMain/active/SinkNoAdvance
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| jeans_ncells_ref    = 12.0               |      |             |
| jeans_ncells_deref  = 24.0               |      |             |
| sink_density_thresh       = 8.52171362932e-22 |      |             |
| sink_accretion_radius     = 3.3753125e+18  |      |             |
| sink_softening_radius     = 3.3753125e+18  |      |             |
| sink_softening_type_gas   = "spline"       |      |             |
| sink_softening_type_sinks = "spline"       |      |             |
| sink_integrator           = "leapfrog"     |      |             |
| sink_subdt_factor         = 0.01           |      |             |
| sink_dt_factor            = 0.5            |      |             |
| sink_merging              = .true.          |      |             |
| sink_maxSinks             = 1000            |      |             |
| sink_convergingFlowCheck  = .true.          |      |             |
| sink_potentialMinCheck    = .true.          |      |             |
| sink_jeansCheck           = .true.          |      |             |
| sink_negativeEtotCheck    = .true.          |      |             |
| sink_GasAccretionChecks   = .true.          |      |             |

## Gravity
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| useGravity      = T               |      |             |
| updateGravity   = .true.           |      |             |
| grav_boundary_type      = "isolated" |      |             |
| grav_boundary_type_x    = "isolated" |      |             |
| grav_boundary_type_y    = "isolated" |      |             |
| grav_boundary_type_z    = "isolated" |      |             |

## Multigrid solver
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| mpole_lmax = 10             |      |             |
| mg_maxResidualNorm = 0.01    |      |             |
| mg_printNorm       = .false. |      |             |

## heating and cooling
| Variable (=default)           | Type | Description |
|:-----------------------------|:-----|:------------|
| useHeat        = T           |      |             |
| he_abundM = 0.0              |      |             |
| he_metal  = 4.               |      |             |
| subfactor       = 0.3        |      |             |
| he_int_method   = "Implicit" |      |             |
| he_pe_recipe    = "WD01"     |      |             |
| he_pe_norm      = 1.3e-24    |      |             |
| theatmin        = 1e-99      |      |             |
| theatmax        = 2.0E4      |      |             |
| tradmin         = 10.0       |      |             |
| tradmax         = 1.E15      |      |             |
| absTmin         = 10.        |      |             |
| absTmax         = 1e9        |      |             |
| dradmin         = 1e-99      |      |             |
| dradmax         = 1e10       |      |             |
| stratifyHeat    = .false.    |      |             |
| h_uv            = 9.25703274e20 |      |             |
| Gzero           = 1.69e0     |      |             |
| use_cr_heating  = .true.     |      |             |
| crIonRate       = 2.0e-17    |      |             |
| crIonExp        = 0.021      |      |             |
| crIonNH         = 1e20       |      |             |
| crIonEnergy     = 20.0e0     |      |             |
| useDustCool     = .true.     |      |             |
| T_cool_min      = 10.0       |      |             |
| nd_cool_min     = 1e1        |      |             |
| nd_cool_max     = 1e10       |      |             |
| dust_sputter_temp = 3e5      |      |             |

## radiative transfer options
| Variable (=default)            | Type | Description |
|:-------------------------------|:-----|:------------|
| useRadTrans    = T             |      |             |
| rt_rayTrace    = T             |      |             |
| useRadTransfer = T             |      |             |
| rt_maxHchange  = 0.1           |      |             |
| rt_ion_threshold = 0.1         |      |             |
| rt_ion_min    = 1e-8           |      |             |
| rt_neutral_min = 1e-10         |      |             |
| rt_dt_temp     = 1.0e5         |      |             |
| ph_sampling     = 2.0          |      |             |
| ph_initHPlevel  = 2            |      |             |
| ph_inBlockSplit = T            |      |             |
| ph_rotRays      = T            |      |             |
| ph_maxNRays     = 1000000      |      |             |
| ph_raysToBundle = 5            |      |             |
| ph_CommCheckInterval = 20      |      |             |
| ph_radPressure  = T            |      |             |
| cfl_radPressure = 0.3          |      |             |
| early_term_FUV  = T            |      |             |
| sigDust         = 1e-21        |      |             |
| dust_gas_ratio  = 0.01         |      |             |
| ph_EUVonDust    = T            |      |             |
| rt_useRadTransDt = .true.      |      |             |
| rt_useNumstepsRadTransDtOnStart = .false. | |         |
| rt_numStepsRadTransDt = 10     |      |             |

## wind options
| Variable (=default)            | Type | Description |
|:-------------------------------|:-----|:------------|
| ref_radius = -1.0              |      |             |
| min_radius = 0.0               |      |             |
| cons_quant = "momentum"        |      |             |
| min_wind_mass = 1.3923e34      |      |             |
| mass_load  = .true.            |      |             |
| var_radius = .false.           |      |             |
| wind_target_temp = 5e6         |      |             |
| perturb_velocity = .false.     |      |             |
| perturb_std_dev  = 0.05        |      |             |
| use_wind_compute_dt = .false.  |      |             |