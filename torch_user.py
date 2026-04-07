#!/usr/bin/env python
"""
User file for torch star formation code.
You must define the methods:
    user_initial_conditions(state, hydro)
    user_parameters()
    flash_par_parameters()
User parameters should have AMUSE units attached, where appropriate.
Design inspired by TRISTAN-MP, Athena++ architecture.
"""
from __future__ import division, print_function
from amuse.datamodel import Particles
from amuse.units import units
from torch_param import FlashPar
from torch_mainloop import run_torch

def flash_par_parameters():
    """
    User-configurable FLASH parameters.
    Define only the parameters you want to change from the template.
    All parameters here will be written to flash.par.
    """
    fp = {}
    
    # ==================================================================
    # Simulation/SimulationMain, monitors/Logfile/LogfileMain, IO/IOMain
    # ==================================================================
    fp['basenm'] = '"turbsph_"'
    fp['log_file'] = '"turbsph.log"'
    fp['run_comment'] = '"Turbulent sphere"'
    fp['output_directory'] = '"./data"'
    
    # ==============================
    # Simulation/SimulationMain/Cube
    # ==============================
    fp['sim_cubeFile'] = '"./cube128"'
    fp['sim_init_Hp'] = '1e-4'
    fp['sim_tdust'] = '18.0'
    fp['sim_A_n'] = '1.2784'
    fp['sim_A_i'] = '0.6698'
    
    # Initial magnetic field
    fp['bx0'] = '0.0e0'
    fp['by0'] = '0.0e0'
    fp['bz0'] = '3.0e-6'
    
    # ===================
    # physics/Eos/EosMain
    # ===================
    fp['gamma'] = '1.66666666666666667'
    
    # =============
    # Grid/GridMain
    # =============
    # Computational domain (for test problem in cube128: -7,+7 pc)
    fp['xmax'] = '2.1602e+19'
    fp['xmin'] = '-2.1602e+19'
    fp['ymax'] = '2.1602e+19'
    fp['ymin'] = '-2.1602e+19'
    fp['zmax'] = '2.1602e+19'
    fp['zmin'] = '-2.1602e+19'
    
    # Boundary types
    fp['xl_boundary_type'] = '"outflow"'
    fp['xr_boundary_type'] = '"outflow"'
    fp['yl_boundary_type'] = '"outflow"'
    fp['yr_boundary_type'] = '"outflow"'
    fp['zl_boundary_type'] = '"outflow"'
    fp['zr_boundary_type'] = '"outflow"'
    
    # Grid/GridMain/paramesh/paramesh4
    fp['nblockx'] = '1'
    fp['nblocky'] = '1'
    fp['nblockz'] = '1'
    fp['lrefine_max'] = '2'
    fp['lrefine_min'] = '2'
    
    fp['refine_var_1'] = '"pres"'
    fp['refine_filter_1'] = '1e-2'
    fp['refine_cutoff_1'] = '0.98'
    fp['derefine_cutoff_1'] = '0.6'
    
    # Derefinement outside rectangular region
    fp['use_deref'] = '.false.'
    fp['deref_lref'] = '2'
    fp['deref_xl'] = '-1.543e+19'
    fp['deref_xr'] = '1.543e+19'
    fp['deref_yl'] = '-1.543e+19'
    fp['deref_yr'] = '1.543e+19'
    fp['deref_zl'] = '-1.543e+19'
    fp['deref_zr'] = '1.543e+19'
    
    fp['gr_sanitizeDataMode'] = '3'
    fp['gr_sanitizeVerbosity'] = '4'
    fp['gr_pmrpConserve'] = '.true.'
    
    # ==================================================
    # Simulation time, output file intervals, restarting
    # ==================================================
    fp['restart'] = '.false.'
    fp['nend'] = '999999999'
    fp['tmax'] = '6.3113852e13'  # 2 Myr
    fp['dtinit'] = '3.15e8'
    fp['dtmin'] = '3.15e7'
    fp['dtmax'] = '3.15e12'
    fp['dr_shortenLastStepBeforeTMax'] = '.false.'
    fp['tstep_change_factor'] = '2.0'
    
    fp['wall_clock_checkpoint'] = '0.432e+05'
    fp['wall_clock_time_limit'] = '1.0e99'
    fp['checkpointFileNumber'] = '0'
    fp['checkpointFileIntervalStep'] = '0'
    fp['checkpointFileIntervalTime'] = '3.1556926e13'  # 1 Myr
    fp['plotFileNumber'] = '0'
    fp['plotFileIntervalStep'] = '0'
    fp['plotFileIntervalTime'] = '3.1556926e12'  # 0.1 Myr
    fp['particleFileNumber'] = '0'
    fp['particleFileIntervalStep'] = '0'
    fp['particleFileIntervalTime'] = '3.1556926e12'  # 0.1 Myr
    
    # Time step limiting for positive definiteness
    fp['dr_usePosdefComputeDt'] = '.true.'
    fp['dr_numPosdefVars'] = '4'
    fp['dr_posdefDtFactor'] = '0.9'
    fp['dr_posdefVar_1'] = '"ener"'
    fp['dr_posdefVar_2'] = '"eint"'
    fp['dr_posdefVar_3'] = '"dens"'
    fp['dr_posdefVar_4'] = '"pres"'
    
    # Variables in plotfile
    fp['plot_var_1'] = '"dens"'
    fp['plot_var_2'] = '"pres"'
    fp['plot_var_3'] = '"temp"'
    fp['plot_var_4'] = '"velx"'
    fp['plot_var_5'] = '"vely"'
    fp['plot_var_6'] = '"velz"'
    fp['plot_var_7'] = '"ihp"'
    fp['plot_var_8'] = '"iha"'
    fp['plot_var_9'] = '"eint"'
    fp['plot_var_10'] = '"ener"'
    fp['plot_var_11'] = '"tdus"'
    fp['plot_var_12'] = '"magx"'
    fp['plot_var_13'] = '"magy"'
    fp['plot_var_14'] = '"magz"'
    fp['plot_var_15'] = '"uvfl"'
    fp['plot_var_16'] = '"fufl"'
    fp['plot_var_17'] = '"auvf"'
    fp['plot_var_18'] = '"afuf"'
    fp['plot_var_19'] = '"magp"'
    fp['plot_var_20'] = '"gpot"'
    fp['plot_var_21'] = '"bgpt"'
    
    # ============
    # Hydro solver
    # ============
    fp['useHydro'] = 'T'
    fp['cfl'] = '0.8'
    fp['eintSwitch'] = '1.e-4'
    fp['killdivb'] = 'T'
    fp['UnitSystem'] = '"CGS"'
    
    # Field variable limits
    fp['smallt'] = '1.0d-99'
    fp['smallp'] = '1.3807e-19'
    fp['smallu'] = '1.0d-99'
    fp['smallx'] = '1.0d-99'
    fp['smalle'] = '1.0d-99'
    fp['smlrho'] = '1.67e-28'
    
    # Unsplit staggered mesh MHD solver
    fp['hydrocomputedtoption'] = '-1'
    fp['order'] = '3'
    fp['slopeLimiter'] = '"minmod"'
    fp['LimitedSlopeBeta'] = '1.'
    fp['charLimiting'] = '.true.'
    
    fp['use_avisc'] = '.true.'
    fp['cvisc'] = '0.1'
    fp['use_flattening'] = '.false.'
    fp['use_steepening'] = '.false.'
    fp['use_upwindTVD'] = '.true.'
    fp['use_gravhalfupdate'] = '.true.'
    
    fp['use_hybridOrder'] = '.true.'
    fp['hy_fallbackLowerCFL'] = '.true.'
    fp['hy_cflFallbackFactor'] = '0.9'
    
    # Magnetic and electric fields
    fp['E_modification'] = '.true.'
    fp['E_upwind'] = '.true.'
    fp['energyFix'] = '.true.'
    fp['ForceHydroLimit'] = '.false.'
    fp['prolMethod'] = '"injection_prol"'
    
    # Riemann solvers
    fp['RiemannSolver'] = '"HLLD"'
    fp['entropy'] = '.false.'
    fp['EOSforRiemann'] = '.false.'
    
    # Strong shock handling
    fp['shockDetect'] = '.true.'
    fp['shockLowerCFL'] = '.true.'
    
    # Super-time-stepping
    fp['useSTS'] = '.false.'
    fp['nstepTotalSTS'] = '5'
    fp['nuSTS'] = '0.2'
    
    # =======================
    # Particles/ParticlesMain
    # =======================
    fp['useParticles'] = 'T'
    fp['useSinkParticles'] = '.true.'
    fp['pt_maxPerProc'] = '10000'
    
    # Particles/ParticlesMain/active/SinkNoAdvance
    fp['jeans_ncells_ref'] = '12.0'
    fp['jeans_ncells_deref'] = '24.0'
    
    fp['sink_density_thresh'] = '8.52171362932e-22'
    fp['sink_accretion_radius'] = '3.3753125e+18'
    fp['sink_softening_radius'] = '3.3753125e+18'
    
    fp['sink_softening_type_gas'] = '"spline"'
    fp['sink_softening_type_sinks'] = '"spline"'
    fp['sink_integrator'] = '"leapfrog"'
    fp['sink_subdt_factor'] = '0.01'
    fp['sink_dt_factor'] = '0.5'
    fp['sink_merging'] = '.true.'
    fp['sink_maxSinks'] = '1000'
    fp['sink_convergingFlowCheck'] = '.true.'
    fp['sink_potentialMinCheck'] = '.true.'
    fp['sink_jeansCheck'] = '.true.'
    fp['sink_negativeEtotCheck'] = '.true.'
    fp['sink_GasAccretionChecks'] = '.true.'
    
    # =============================================
    # Gravity
    # =============================================
    fp['useGravity'] = 'T'
    fp['updateGravity'] = '.true.'
    fp['grav_boundary_type'] = '"isolated"'
    
    # Multigrid solver
    fp['mpole_lmax'] = '10'
    fp['mg_maxResidualNorm'] = '0.01'
    fp['mg_printNorm'] = '.false.'
    
    # ===================================================
    # Heating and cooling
    # ===================================================
    fp['useHeat'] = 'T'
    
    fp['he_abundM'] = '0.0'
    fp['he_metal'] = '4.'
    
    fp['subfactor'] = '0.3'
    fp['he_int_method'] = '"Implicit"'
    
    # Photoelectric heating
    fp['he_pe_recipe'] = '"WD01"'
    fp['he_pe_norm'] = '1.3e-24'
    
    # Heating and cooling thresholds
    fp['theatmin'] = '1e-99'
    fp['theatmax'] = '2.0E4'
    fp['tradmin'] = '10.0'
    fp['tradmax'] = '1.E15'
    fp['absTmin'] = '10.'
    fp['absTmax'] = '1e9'
    fp['dradmin'] = '1e-99'
    fp['dradmax'] = '1e10'
    
    fp['stratifyHeat'] = '.false.'
    fp['h_uv'] = '9.25703274e20'
    fp['Gzero'] = '1.69e0'
    
    # Cosmic ray heating
    fp['use_cr_heating'] = '.true.'
    fp['crIonRate'] = '2.0e-17'
    fp['crIonExp'] = '0.021'
    fp['crIonNH'] = '1e20'
    fp['crIonEnergy'] = '20.0e0'
    
    # Molecular + dust cooling
    fp['useDustCool'] = '.true.'
    fp['T_cool_min'] = '10.0'
    fp['nd_cool_min'] = '1e1'
    fp['nd_cool_max'] = '1e10'
    fp['dust_sputter_temp'] = '3e5'
    
    # ================================================
    # Radiative transfer
    # ================================================
    fp['useRadTrans'] = 'T'
    fp['rt_rayTrace'] = 'T'
    
    fp['rt_maxHchange'] = '0.1'
    fp['rt_ion_threshold'] = '0.1'
    fp['rt_ion_min'] = '1e-8'
    fp['rt_neutral_min'] = '1e-10'
    fp['rt_dt_temp'] = '1.0e5'
    
    fp['ph_sampling'] = '2.0'
    fp['ph_initHPlevel'] = '2'
    fp['ph_inBlockSplit'] = 'T'
    fp['ph_rotRays'] = 'T'
    fp['ph_maxNRays'] = '1000000'
    fp['ph_raysToBundle'] = '5'
    fp['ph_CommCheckInterval'] = '20'
    fp['ph_radPressure'] = 'T'
    fp['cfl_radPressure'] = '0.3'
    
    fp['early_term_FUV'] = 'T'
    fp['sigDust'] = '1e-21'
    fp['dust_gas_ratio'] = '0.01'
    fp['ph_EUVonDust'] = 'T'
    
    fp['rt_useRadTransDt'] = '.true.'
    fp['rt_useNumstepsRadTransDtOnStart'] = '.false.'
    fp['rt_numStepsRadTransDt'] = '10'
    
    # ======================================================================================
    # Wind options
    # ======================================================================================
    fp['ref_radius'] = '-1.0'
    fp['min_radius'] = '0.0'
    fp['cons_quant'] = '"momentum"'
    fp['min_wind_mass'] = '1.3923e34'
    fp['mass_load'] = '.true.'
    fp['var_radius'] = '.false.'
    fp['wind_target_temp'] = '5e6'
    fp['perturb_velocity'] = '.false.'
    fp['perturb_std_dev'] = '0.05'
    fp['use_wind_compute_dt'] = '.false.'
    
    return fp


def get_ntasks_from_run_script(name="run.sh"):
    """formally -n is --ntasks, de facto same as nprocs"""
    n = None
    nodes = None
    cores = None
    import os
    # Check for slurm ntasks
    n = int(os.getenv("SLURM_NTASKS"))
    if n is None:
        with open(name) as f:
            for line in f:
                w = line.split()
                if len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('--ntasks-per-node'):
                    assert cores is None  # throw error if #SBATCH -n occurs >1x
                    cores = int(''.join(char for char in w[1] if char.isdigit()))
                elif len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('-N'):
                    assert nodes is None  # throw error if #SBATCH -N or --nodes occurs >1x
                    nodes = int(''.join(char for char in w[1] if char.isdigit()))
                elif len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('--nodes'):
                    assert nodes is None  # throw error if #SBATCH -N or --nodes occurs >1x
                    nodes = int(''.join(char for char in w[1] if char.isdigit()))
                elif len(w) >= 2 and w[0] == '#SBATCH' and w[1].startswith('-n'):
                    assert n is None
                    n = int(''.join(char for char in w[1] if char.isdigit()))
        if n is None:
            n = nodes*cores
    assert n is not None
    return n


def user_initial_conditions(state, hydro):
    """
    User-provided method to set initial conditions for the simulation.
    Usually, this means adding star particles to the hydro code.
    We add stars to hydro only, not other Particles() structures such as
    state.stars or grav.particles.  The method torch.evolve(...) copies
    particles from hydro to other workers before it starts the evolution loop.
    """
    # ------------------------------------------------------------------------
    # Single star  test: plop a single star at the center of the box
    # ------------------------------------------------------------------------
    #    flashp = FlashPar("flash.par")
    #
    #    star          = Particles(1)
    #    star.mass     = 1 | units.MSun
    #    star.position = [0, 0, 0] | units.cm
    #    star.velocity = [0, 0, 0] | units.cm/units.s
    #
    #    star_tag = hydro.add_particles(star.x, star.y, star.z)
    #    hydro.set_particle_mass(star_tag, star.mass)
    #    hydro.set_particle_velocity(star_tag, star.vx, star.vy, star.vz)
    #    hydro.set_particle_oldmass(star_tag, star.mass) # Save initial stellar mass for SE code.
    # NOTE: This only works with a ZAMS star -- other parameters are needed to start 
    # with an evolved star.
    # ------------------------------------------------------------------------
    pass


def user_parameters():
    """
    User configurable parameters.  All parameters are currently required.
    """
    # First, generate flash.par from user-defined parameters
    flash_params = flash_par_parameters()
    write_flash_par(flash_params)
    
    p = {}
    flashp = FlashPar("flash.par")
    
    # <VorAMR>
    try:
        p['with_voramr'] = flashp['use_voramr']
    except KeyError:
        p['with_voramr'] = False
    
    if p['with_voramr']:    
        p['source_file'] = flashp['voramr_source']
        p['convert_file'] = True
        p['use_localRef'] = flashp['use_localRef']
        p['local_ref'] = [flashp['localRef_x'], flashp['localRef_y'], flashp['localRef_z'], flashp['localRef_r']]
        p['center_local_ref'] = flashp['center_localRef']
        p['input_file'] = flashp['voramr_input']
        p['pickle_kdtree'] = False
        p['pickle_file_name'] = "kdtree.pickle"
        p['numBlocks'] = 15000
        p['cellsPerBlock'] = 16
    
    # <bridge>
    p['npy_seed'] = 0
    p['restart_with_new_rng'] = False
    p['restart_with_user_ics'] = False
    p['evolve_async'] = True
    p['with_bridge'] = True
    p['with_multiples'] = True
    p['with_se'] = True
    p['remove_merged'] = True and p['with_se']
    
    # <timestepping>
    p['hy_dt_factor'] = 0.99999
    
    # <star/n-body gravity>
    p['with_ph4'] = True
    p['epsilon'] = 15.0 | units.RSun
    
    # <star/n-body gravity & binaries>
    p['with_petar'] = True
    p['r_bin'] = 50 | units.au
    p['r_out'] = 0.003 | units.pc
    p['dt_soft_max'] = 0.125 | units.kyr
    
    # <stellar evolution>
    p['with_lyc'] = True
    p['with_pe_heat'] = True
    p['sigd'] = float(flash_params['sigDust'])
    p['with_sn'] = True
    p['with_winds'] = True
    p['massloss_method'] = 'puls'
    p['min_feedback_mass'] = 7.0 | units.MSun
    p['remove_merged'] = True
    
    # <star particle creation>
    p['binaries'] = False
    p['mult_frac'] = 'field'
    p['pdist'] = 'inner'
    p['qdist'] = 'field'
    p['edist'] = 'field'
    p['min_imf_mass'] = 0.08 | units.MSun
    p['max_imf_mass'] = 100.0 | units.MSun
    p['sample_imf_mass'] = 10000.0 | units.MSun
    p['sample_imf_bins'] = 100
    p['sink_rad'] = float(flash_params['sink_accretion_radius']) | units.cm
    p['sum_small'] = False
    p['m_small'] = 1.0 | units.MSun
    
    # <amuse file overwrite>
    p['overwrite'] = True
    
    # <job>
    ntasks = get_ntasks_from_run_script("run.sh")
    p['num_grav_workers'] = 1
    p['num_hy_workers'] = ntasks - p['num_grav_workers'] - 1
    
    if p['with_petar']:
        p['with_ph4'] = False
        p['with_multiples'] = False
        p['epsilon'] = 0 | units.RSun
    
    if p['with_se']:
        p['num_hy_workers'] -= 1
    
    if p['with_multiples']:
        p['num_hy_workers'] -= 2
    
    return p


def write_flash_par(flash_params, output_file="flash.par"):
    """
    Generate flash.par file from user-defined parameters.
    
    Args:
        flash_params: Dictionary of parameter name -> value pairs
        output_file: Path to write the generated flash.par (default: "flash.par")
    """
    import os
    
    # Process each line
    output_lines = []
    for param, value in flash_params.items():
        line = f"{param} = {value}\n"
        output_lines.append(line)
    
    # Write the output file
    with open(output_file, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Generated {output_file} with user-defined parameters")

# ============================================================================
if __name__ == '__main__':
    run_torch(
        user_initial_conditions,
        user_parameters,
    )
