#!/usr/bin/env python
"""
Stellar evolution module for Torch code

Includes subroutines to

* Compute information about star luminosity (energy, photon count, etc) at
  requested wavelength and temperature

* Compute stellar wind mass loss rate, velocity
"""

from __future__ import division, print_function

import numpy as np
from scipy.integrate import quad

from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.ext.orbital_elements import *

from ionizingflux import ionizing_photon_flux
from torch_stdout import tprint

h = 6.6261e-27 # Planck's constant
c = 2.9979e10  # Speed of light
k = 1.3807e-16 # Boltzmann constant

sigSB = 5.6704e-5 | (units.g/((units.s)**3 * (units.K)**4)) # Stefan-Boltzmann constant, g s^-3 K^-4, CCC 26/04/2024
sig0 = 6.304e-18 # Photoionization cross section at threshold for hydrogen
E_ev = 1.60222497096e-12 # energy of 1 eV in erg
E_lyc = 13.6*E_ev  # 13.6 eV
# Cross section for dust per hydrogen atom.
# Value = tau / N_H where tau = gamma * Av (Draine and Bertoli 96)
# Av = N_H,tot / (1.87e21 cm^2) (Bohlin et al 78)
# gamma = 2.5 (Bergin et al 2004)

def binary_evolution(time, dt, se_restart_time, state, hydro, se,
                     with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                     massloss_method=None, min_feedback_mass=None, CE_method='wind', CE_alpha=1):
    """
    NOTE: time = target time to evolve TO, including the dt already.
    Chosen to follow AMUSE worker convention.
    """
    assert massloss_method is not None
    assert min_feedback_mass is not None
    
    # Set radius to physical radius for restart with user ICs
    # This assumes the stars are ZAMS, which may be incorrect 
    _attributes = state.stars.get_attribute_names_defined_in_store()
    if 'radius' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # It must be somewhat realistic in case there is a contact system
        # Empirical relation from https://articles.adsabs.harvard.edu/pdf/1991Ap%26SS.181..313D
        # Use linear MRR for upper mass range
        state.stars.radius = (1.01 * (state.stars.mass / (1 | units.MSun)) ** 0.57) | units.RSun
    if 'luminosity' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # Empirical relation from https://articles.adsabs.harvard.edu/pdf/1991Ap%26SS.181..313D
        # Use linear MLR for upper mass range
        state.stars.luminosity = (1.15 * (state.stars.mass / (1 | units.MSun)) ** 3.36) | units.LSun
    if 'temperature' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # Use BB luminosity and radius, luminosity
        state.stars.temperature = (state.stars.luminosity / (4 * np.pi * sigSB))**(1./4) * state.stars.radius**(-1./2)
    if "wind_mass_loss_rate" not in _attributes:
        # Save wind mass loss rate for comparison - CCC 09/04/2025
        # Set to 0 initially
        state.stars.wind_mass_loss_rate = np.zeros(len(state.stars)) | units.MSun / units.yr 

    
    # Update ALL the star properties in bulk for consistency.
    # Copy the old mass for the wind velocities (calculated outside of amuse) - CCC 04/11/2023
    old_mass = np.copy(state.stars.mass)
    
    dm_dt   = np.zeros(len(state.stars)) | units.g / units.s
    vterm   = np.zeros(len(state.stars)) | units.cm / units.s
    nion    = np.zeros(len(state.stars)) | units.s**-1
    eion    = np.zeros(len(state.stars)) | units.erg
    sigh    = np.zeros(len(state.stars)) | units.cm**2
    npe     = np.zeros(len(state.stars)) | units.s**-1
    epe     = np.zeros(len(state.stars)) | units.erg
    sigpe   = np.zeros(len(state.stars)) | units.cm**2

    # follow FLASH idiom; return dt after SN deposit
    se_dt = 1e99 | units.s

    # Time since star formation or restart, to use for stellar evolution
    se_time = time - max(se_restart_time, min(hydro.get_particle_creation_time(state.stars.tag)))
    
    # Pass information to SE
    # Check for new systems, important for binaries - CCC 02/08/2024
    state.binaries.synchronize_to(se.binaries)
    # Now pass attributes to binaries - CCC 03/08/2024
    for _attribute in state.binaries.get_attribute_names_defined_in_store():
        if _attribute in se.binaries.get_attribute_names_defined_in_store():
            state.binaries_to_se.copy_attributes([_attribute])
    # Evolve model
    se.evolve_model(se_time) #Time attached to se.particles.age, which is the "simulation" time
    
    # Pass information back to stars after end of SE loop - CCC 02/08/2024
    for _attribute in se.particles.get_attribute_names_defined_in_store():
        if _attribute in state.stars.get_attribute_names_defined_in_store():
            state.se_to_stars.copy_attributes([_attribute])

    # Reset the stars' age after the SE step
    state.stars.age = time - hydro.get_particle_creation_time(state.stars.tag)

    # Turn off wind if mass increased - CCC 11/09/2024 
    # Keep a list of stars for which feedback was calculated - CCC 25/10/2024
    evolved_stars = Particles()
    # Also keep a list of merged stars - CCC 08/01/2025
    merged_stars = Particles()
    # List primaries and companions - CCC 25/10/2024
    primaries = Particles()
    companions = Particles()
    if len(state.binaries) > 0:
        for binary in state.binaries:
            primaries.add_particle(binary.child1)
            companions.add_particle(binary.child2)

        # Update positions and velocities from binary interaction - CCC 03/04/2025
        orbital_elements = get_orbital_elements_from_binaries(primaries, companions)
    
    for i, s in enumerate(state.stars):
        
        in_binary = False # Set to true if in binary - CCC 25/10/2024
        
        if s in evolved_stars: # CCC 25/10/24
            continue
            
        if (s in primaries): # CCC 25/10/24
            in_binary = True
            # Find the binary
            j = np.where(primaries.tag == s.tag)[0]
            b = state.binaries[j]
            # and the other star
            k = np.where(state.stars.tag == companions[j].tag)[0][0]
            t = state.stars[k]
        
        if (s in companions):
            continue # Force choice of primary and companion for change in orbit - CCC 06/04/2025

        if went_supernova(s.stellar_type):
            continue

        if s.initial_mass >= min_feedback_mass:

            if with_sn and went_supernova_from_kick(s, state.stars): # Do not check if in binary for SN - CCC 29/10/2024

                inj_mass = old_mass[i] - s.mass  # minus stellar remnant's mass
                
                if inj_mass > 15.0|units.MSun:
                    # expected upper limit for SeBa tracks; see
                    # https://groups.google.com/forum/#!topic/torch-users/rWJd6l_mRBg/discussion
                    tprint("... flooring SN inj_mass {} MSun to 15 MSun".format(inj_mass.value_in(units.MSun)))
                    inj_mass = 15.0|units.MSun

                # Inject energy and mass onto grid                         
                # In SeBa, stars with CO core mass above 15 Msun are direct collapse, so don't inject SN          
                if s.COcore_mass <= 15 | units.MSun:
                    _tmp = hydro.energy_injection(1e51|units.erg, -1.0, inj_mass.in_(units.g), s.x, s.y, s.z)
                    se_dt = min(se_dt, _tmp)
                    tprint("... SN x={}, y={}, z={}, inj_mass={}, tag={}".format(s.x, s.y, s.z, inj_mass.value_in(units.MSun), s.tag))

                # implicitly zeros out feedback properties by not setting

                # Update velocity from kick velocity - CCC 06/04/2025
                s.vx += s.natal_kick_x
                s.vy += s.natal_kick_y
                s.vz += s.natal_kick_z
                
                # Set star to evolved after SN - CCC 29/10/2024
                evolved_stars.add_particle(s)
                
            # Check if companion went supernova - CCC 15/04/2025
            if in_binary and with_sn and went_supernova_from_kick(t, state.stars):
                
                inj_mass = old_mass[k] - t.mass  # minus stellar remnant's mass
                
                if inj_mass > 15.0|units.MSun:
                    # expected upper limit for SeBa tracks; see
                    # https://groups.google.com/forum/#!topic/torch-users/rWJd6l_mRBg/discussion
                    tprint("... flooring SN inj_mass {} MSun to 15 MSun".format(inj_mass.value_in(units.MSun)))
                    inj_mass = 15.0|units.MSun

                # Inject energy and mass onto grid                                                                                                                                                             
                # In SeBa, stars with CO core mass above 15 Msun are direct collapse, so don't inject SN                                                                                                       
                if t.COcore_mass <= 15 | units.MSun:
                    _tmp = hydro.energy_injection(1e51|units.erg, -1.0, inj_mass.in_(units.g), t.x, t.y, t.z)
                    se_dt = min(se_dt, _tmp)
                    tprint("... SN x={}, y={}, z={}, inj_mass={}, tag={}".format(t.x, t.y, t.z, inj_mass.value_in(units.MSun), t.tag))
                
                # implicitly zeros out feedback properties by not setting

                # Update velocity from kick velocity - CCC 06/04/2025
                t.vx += t.natal_kick_x
                t.vy += t.natal_kick_y
                t.vz += t.natal_kick_z
                
                # Set star to evolved after SN - CCC 29/10/2024
                evolved_stars.add_particle(t)
                
            else: # Check if in binary to determine the feedback mechanisms - CCC 29/10/2024
                
                if in_binary:

                    be = se.binaries[j] # Evolved binary

                    # Energy after minus energy before --> should be > 0 for CE - CCC 08/04/2025
                    _dE = (units.constants.G / 2) * (((s.mass * t.mass / be.semi_major_axis) - old_mass[i]*old_mass[k] / b.semi_major_axis))[0] # Try here - CCC 08/12/2024
                    inj_mass = old_mass[i] + old_mass[k] - s.mass - t.mass
                    # If _dE > 0, mass transfer or CE; if _dE < 0, wind mass loss

                    # Different ways to detect the interaction - CCC 06/04/2025
                    if (_dE > (0 | units.erg)) or (be.binary_type > 2) or (be.semi_major_axis < b.semi_major_axis) \
                    or accreted_mass(s.mass, old_mass[i]) or accreted_mass(t.mass, old_mass[k]): 
                        
                        tprint('... change in orbit from BE')

                        # Compare to wind mass loss rate - CCC 06/04/2025
                        dm_dt_wind_1, vterm_wind_1 = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                                        massloss_method=massloss_method)
                        dm_dt_wind_2, vterm_wind_2 = compute_dmdt_vterm(old_mass[k], t.temperature, t.radius, t.mass, t.luminosity, dt,
                                                                        massloss_method=massloss_method)

                        # If lost mass in excess of wind mass loss rate - CCC 06/04/2025
                        # Wind mass loss rate from SeBa is negative by definition
                        if dm_dt_wind_1 > -1*s.wind_mass_loss_rate or dm_dt_wind_2 > -1*t.wind_mass_loss_rate: 
                            tprint('... mass transfer')
                            tprint('Injected mass:', "{0:.2f}".format(inj_mass.value_in(units.MSun)), 'MSun')

                            # https://www.aanda.org/articles/aa/full_html/2021/04/aa40442-21/aa40442-21.html
                            # CCC 12/09/2024, 20/06/2023
                            E_bind = CE_alpha * _dE
                            #tprint('Ejecta energy:', "{0:.1e}".format(E_bind.value_in(units.erg)), 'erg')
                        
                        
                            if CE_method=='wind':
                            
                                if (old_mass[i] - s.mass) > (old_mass[k] - t.mass):
                                    # If donor is star s
                                    dm_dt[i] = inj_mass/dt
                                    vterm[i] = vterm_wind_1 # Cap the ejecta velocity at the wind velocity
                                    tprint('Ejecta velocity from wind', "{0:.1e}".format(vterm[i].value_in(units.km/units.s)), 'km/s')
                                else:
                                    dm_dt[k] = inj_mass/dt
                                    vterm[k] = vterm_wind_2 # Cap the ejecta velocity at the wind velocity
                                    tprint('Ejecta velocity from wind', "{0:.1e}".format(vterm[k].value_in(units.km/units.s)), 'km/s')
                            
                            # If there is no energy from the envelope ejection, do wind mass loss
                            # for the donor star instead
                            elif CE_method=='alpha':
                            
                                if (old_mass[i] - s.mass) > (old_mass[k] - t.mass):
                                    # If donor is star s
                                    dm_dt[i] = inj_mass/dt
                                    if E_bind > 0 | units.erg:
                                        vterm[i] = compute_vterm_binary(inj_mass, E_bind)
                                        tprint('Ejecta velocity from change in energy', "{0:.1e}".format(vterm[i].value_in(units.km/units.s)), 'km/s')
                                    else:
                                        vterm[i] = vterm_wind_1 # Cap the ejecta velocity at the wind velocity
                                        tprint('Ejecta velocity from wind', "{0:.1e}".format(vterm[i].value_in(units.km/units.s)), 'km/s')
                                else:
                                    dm_dt[k] = inj_mass/dt
                                    if E_bind > 0 | units.erg:
                                        vterm[k] = compute_vterm_binary(inj_mass, E_bind)
                                        tprint('Ejecta velocity from change in energy', "{0:.1e}".format(vterm[k].value_in(units.km/units.s)), 'km/s')
                                    else:
                                        vterm[k] = vterm_wind_2 # Cap the ejecta velocity at the wind velocity
                                        tprint('Ejecta velocity from wind', "{0:.1e}".format(vterm[k].value_in(units.km/units.s)), 'km/s')
                        
                            elif CE_method=='SN':
                                _tmp = hydro.energy_injection(E_bind, -1.0, inj_mass.in_(units.g), s.x, s.y, s.z)
                
                        else:
                            dm_dt[i] = dm_dt_wind_1
                            vterm[i] = vterm_wind_1
                            dm_dt[k] = dm_dt_wind_2
                            vterm[k] = vterm_wind_2

                        if with_lyc:
                            _tmp = compute_eion_nion_sigh(s.mass, s.temperature, s.radius)
                            eion[i] = _tmp[0]
                            nion[i] = _tmp[1]
                            sigh[i] = _tmp[2]
                            _tmp = compute_eion_nion_sigh(t.mass, t.temperature, t.radius)
                            eion[k] = _tmp[0]
                            nion[k] = _tmp[1]
                            sigh[k] = _tmp[2]
                        if with_pe_heat:
                            _tmp = compute_epe_npe(s.temperature, s.radius)
                            epe[i] = _tmp[0]
                            npe[i] = _tmp[1]
                            sigpe[i] = state.user['sigd'] | units.cm**2
                            _tmp = compute_epe_npe(t.temperature, t.radius)
                            epe[k] = _tmp[0]
                            npe[k] = _tmp[1]
                            sigpe[k] = state.user['sigd'] | units.cm**2
                        # Do not set wind properties for stars that have lost mass due to CE

                        # Update position and velocity - CCC 04/04/2025
                        stars_for_COM = Particles()
                        stars_for_COM.add_particle(primaries[j])
                        stars_for_COM.add_particle(companions[j])
                        COM = stars_for_COM.center_of_mass()
                        COV = stars_for_COM.center_of_mass_velocity()

                        rel_pos, rel_vel = rel_posvel_arrays_from_orbital_elements(s.mass, t.mass,
                                                                                   be.semi_major_axis, be.eccentricity,
                                                                                   orbital_elements[4][j], orbital_elements[5][j],
                                                                                   orbital_elements[6][j], orbital_elements[7][j])
                        m1_f = 1./(1 + s.mass/t.mass)                                                                                                                                     
                        m2_f = 1./(1 + t.mass/s.mass)
                        s.position = COM - m1_f*rel_pos[0]
                        t.position = COM + m2_f*rel_pos[0]
                        s.velocity = COV - m1_f*rel_vel[0]
                        t.velocity = COV + m2_f*rel_vel[0]
                        
                        # Set star to evolved after uMT/CE - CCC 22/11/2024
                        evolved_stars.add_particle(s)
                        evolved_stars.add_particle(t)
                        
                    else: # Normal feedback if the binary does not interact
                            
                        if with_winds:
                            _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                      massloss_method=massloss_method)
                            dm_dt[i] = _tmp[0]
                            vterm[i] = _tmp[1]
                            _tmp = compute_dmdt_vterm(old_mass[k], t.temperature, t.radius, t.mass, t.luminosity, dt,
                                                      massloss_method=massloss_method)
                            dm_dt[k] = _tmp[0]
                            vterm[k] = _tmp[1]
                        if with_lyc:
                            _tmp = compute_eion_nion_sigh(s.mass, s.temperature, s.radius)
                            eion[i] = _tmp[0]
                            nion[i] = _tmp[1]
                            sigh[i] = _tmp[2]
                            _tmp = compute_eion_nion_sigh(t.mass, t.temperature, t.radius)
                            eion[k] = _tmp[0]
                            nion[k] = _tmp[1]
                            sigh[k] = _tmp[2]
                        if with_pe_heat:
                            _tmp = compute_epe_npe(s.temperature, s.radius)
                            epe[i] = _tmp[0]
                            npe[i] = _tmp[1]
                            sigpe[i] = state.user['sigd'] | units.cm**2
                            _tmp = compute_epe_npe(t.temperature, t.radius)
                            epe[k] = _tmp[0]
                            npe[k] = _tmp[1]
                            sigpe[k] = state.user['sigd'] | units.cm**2
                    
                        # Set star to evolved after feedback
                        evolved_stars.add_particle(s)
                        evolved_stars.add_particle(t)
                    
                else: # Normal feedback if the star is not in a binary
                    
                    if with_lyc:
                        _tmp = compute_eion_nion_sigh(s.mass, s.temperature, s.radius)
                        eion[i] = _tmp[0]
                        nion[i] = _tmp[1]
                        sigh[i] = _tmp[2]
                    if with_pe_heat:
                        _tmp = compute_epe_npe(s.temperature, s.radius)
                        epe[i] = _tmp[0]
                        npe[i] = _tmp[1]
                        sigpe[i] = state.user['sigd'] | units.cm**2
                    if with_winds:
                        _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                  massloss_method=massloss_method)
                        dm_dt[i] = _tmp[0]
                        vterm[i] = _tmp[1]
                    
                    # Set star to evolved after feedback - CCC 29/10/2024
                    evolved_stars.add_particle(s)           

        elif (s.initial_mass < min_feedback_mass) and in_binary and (t.mass == (0. | units.MSun)):

            # Update position and velocity - CCC 04/04/2025
            tprint("... low-mass stars merged from BE")
            stars_for_COM = Particles()
            stars_for_COM.add_particle(primaries[j])
            stars_for_COM.add_particle(companions[j])
            s.position = stars_for_COM.center_of_mass()
            s.velocity = stars_for_COM.center_of_mass_velocity()

            # Set stars to evolved after change in orbit
            evolved_stars.add_particle(s)
            evolved_stars.add_particle(t)

            # Save tag of star it merged with
            t.merged_with = s.tag
            # Save merged time
            t.merger_time = hydro.get_time()
            merged_stars.add_particle(t)
            
    
    # Binaries sync'ed to stars later, do not sync here - CCC 12/09/2024

    # This assumes steps are relatively small in the mass loss rate of stars,
    # so that gravity can use the mass after all the wind mass loss has
    # occcured. Otherwise we'd have to average mass loss and keep up with old
    # and new masses and it just gets ugly.

    hydro.set_particle_mass(state.stars.tag, state.stars.mass)

    # TODO not sure if as_quantity_in(...) calls are actually needed.
    # FLASH worker has its own unit converter.  -AT, 2019Oct14
    hydro.set_particle_nion(state.stars.tag, nion)
    hydro.set_particle_eion(state.stars.tag, eion.as_quantity_in(units.erg))
    hydro.set_particle_sigh(state.stars.tag, sigh)

    hydro.set_particle_npep(state.stars.tag, npe)
    hydro.set_particle_epep(state.stars.tag, epe.as_quantity_in(units.erg)) # Set average energy of PE photon
    hydro.set_particle_sigd(state.stars.tag, sigpe) # Set cross section of dust to PE photons.

    hydro.set_particle_wind_mass(state.stars.tag, dm_dt.as_quantity_in(units.g/units.s))
    hydro.set_particle_wind_vel(state.stars.tag, vterm.as_quantity_in(units.cm/units.s))
    
    # Set SeBa properties for checkpoint - CCC 26/04/2024
    hydro.set_particle_rel_mass(state.stars.tag, state.stars.relative_mass)
    hydro.set_particle_rel_age(state.stars.tag, state.stars.relative_age)
    hydro.set_particle_co_corem(state.stars.tag, state.stars.COcore_mass)
    hydro.set_particle_corem(state.stars.tag, state.stars.core_mass)
    hydro.set_particle_radius(state.stars.tag, state.stars.radius)
    hydro.set_particle_stype(state.stars.tag, state.stars.stellar_type.value_in(units.stellar_type))

    if len(merged_stars) > 0:
        # Write merged stars
        state.out_merged_stars(merged_stars, overwrite=True)
        # Remove merged star from hydro
        hydro.remove_particles(merged_stars.tag)
        # Remove from SE
        se.particles.remove_particles(merged_stars)
        # Synchronize to state
        se.particles.synchronize_to(state.stars)

    return se_dt


def stellar_evolution(time, dt, se_restart_time, state, hydro, se,
    with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                      massloss_method=None, min_feedback_mass=None):
    """
    NOTE: time = target time to evolve TO, including the dt already.
    Chosen to follow AMUSE worker convention.
    """
    assert massloss_method is not None
    assert min_feedback_mass is not None
    
    # Set radius to physical radius for restart with user ICs
    # This assumes the stars are ZAMS, which may be incorrect 
    _attributes = state.stars.get_attribute_names_defined_in_store()
    if 'radius' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # It must be somewhat realistic in case there is a contact system
        # Empirical relation from https://articles.adsabs.harvard.edu/pdf/1991Ap%26SS.181..313D
        # Use linear MRR for upper mass range
        state.stars.radius = (1.01 * (state.stars.mass / (1 | units.MSun)) ** 0.57) | units.RSun
    if 'luminosity' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # Empirical relation from https://articles.adsabs.harvard.edu/pdf/1991Ap%26SS.181..313D
        # Use linear MLR for upper mass range
        state.stars.luminosity = (1.15 * (state.stars.mass / (1 | units.MSun)) ** 3.36) | units.LSun
    if 'temperature' not in _attributes:
        # Initial guess for the radius if running with user ICs - CCC 12/05/2023
        # Use BB luminosity and radius, luminosity
        state.stars.temperature = (state.stars.luminosity / (4 * np.pi * sigSB))**(1./4) * state.stars.radius**(-1./2)

    # Update ALL the star properties in bulk for consistency.
    # Copy the old mass for the mass loss rates (calculated outside of amuse) - CCC 04/11/2023
    old_mass = np.copy(state.stars.mass)
    
    dm_dt   = np.zeros(len(state.stars)) | units.g / units.s
    vterm   = np.zeros(len(state.stars)) | units.cm / units.s
    nion    = np.zeros(len(state.stars)) | units.s**-1
    eion    = np.zeros(len(state.stars)) | units.erg
    sigh    = np.zeros(len(state.stars)) | units.cm**2
    npe     = np.zeros(len(state.stars)) | units.s**-1
    epe     = np.zeros(len(state.stars)) | units.erg
    sigpe   = np.zeros(len(state.stars)) | units.cm**2

    # follow FLASH idiom; return dt after SN deposit
    se_dt = 1e99 | units.s

    # Time since star formation or restart, to use for stellar evolution
    se_time = time - max(se_restart_time, min(hydro.get_particle_creation_time(state.stars.tag)))
    
    # make list of remnant stars so we don't explode them again
    remnants = state.stars.tag[went_supernova(state.stars.stellar_type)]

    # Use evolve_model to evolve all stars at the same time
    # This allows us to restart from evolved stars and use the same structure for
    # binary evolution - CCC 26/04/2024
    state.stars_to_se.copy()
    se.evolve_model(se_time)
    state.se_to_stars.copy()
    
    # Reset the stars' age after the SE step, as the SeBa age is reset to 0
    # at each restart - CCC 22/11/2024
    state.stars.age = time - hydro.get_particle_creation_time(state.stars.tag)

    # Loop only over active stars while retaining the correct indexing for total star array
    for i, s in enumerate(state.stars):

        if s.tag in remnants or s.initial_mass < min_feedback_mass:
            continue

        if with_sn and went_supernova(s.stellar_type):

            inj_mass = old_mass[i] - s.mass  # minus stellar remnant's mass
            if inj_mass > 15.0|units.MSun:
                # expected upper limit for SeBa tracks; see
                # https://groups.google.com/forum/#!topic/torch-users/rWJd6l_mRBg/discussion
                tprint("... setting maximum SN inj_mass {} MSun to 15 MSun".format(inj_mass.value_in(units.MSun)))
                inj_mass = 15.0|units.MSun

            # Inject energy and mass onto grid
            # In SeBa, stars with CO core mass above 15 Msun are direct collapse, so don't inject SN
            if s.COcore_mass <= 15 | units.MSun:
                _tmp = hydro.energy_injection(1e51|units.erg, -1.0, inj_mass.in_(units.g), s.x, s.y, s.z)
                se_dt = min(se_dt, _tmp)
                tprint("... SN x={}, y={}, z={}, inj_mass={}, tag={}".format(s.x, s.y, s.z, inj_mass.value_in(units.MSun), s.tag))

            # implicitly zeros out feedback properties by not setting

        else:

            if with_lyc:
                _tmp = compute_eion_nion_sigh(s.mass, s.temperature, s.radius)
                eion[i] = _tmp[0]
                nion[i] = _tmp[1]
                sigh[i] = _tmp[2]
            if with_pe_heat:
                _tmp = compute_epe_npe(s.temperature, s.radius)
                epe[i] = _tmp[0]
                npe[i] = _tmp[1]
                sigpe[i] = state.user['sigd'] | units.cm**2
            if with_winds:
                _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                          massloss_method=massloss_method)
                dm_dt[i] = _tmp[0]
                vterm[i] = _tmp[1]

        # Evolutionary things besides winds could have reduced the stars mass.
        if dm_dt[i]*dt > 0.0|units.MSun:
            s.mass = min(s.mass, old_mass[i] - dm_dt[i]*dt)

    # This assumes steps are relatively small in the mass loss rate of stars,
    # so that gravity can use the mass after all the wind mass loss has
    # occcured. Otherwise we'd have to average mass loss and keep up with old
    # and new masses and it just gets ugly.

    hydro.set_particle_mass(state.stars.tag, state.stars.mass)

    # TODO not sure if as_quantity_in(...) calls are actually needed.
    # FLASH worker has its own unit converter.  -AT, 2019Oct14
    hydro.set_particle_nion(state.stars.tag, nion)
    hydro.set_particle_eion(state.stars.tag, eion.as_quantity_in(units.erg))
    hydro.set_particle_sigh(state.stars.tag, sigh)

    hydro.set_particle_npep(state.stars.tag, npe)
    hydro.set_particle_epep(state.stars.tag, epe.as_quantity_in(units.erg)) # Set average energy of PE photon
    hydro.set_particle_sigd(state.stars.tag, sigpe) # Set cross section of dust to PE photons.

    hydro.set_particle_wind_mass(state.stars.tag, dm_dt.as_quantity_in(units.g/units.s))
    hydro.set_particle_wind_vel(state.stars.tag, vterm.as_quantity_in(units.cm/units.s))

    # Set SeBa properties for checkpoint - CCC 26/04/2024
    hydro.set_particle_rel_mass(state.stars.tag, state.stars.relative_mass)
    hydro.set_particle_rel_age(state.stars.tag, state.stars.relative_age)
    hydro.set_particle_co_corem(state.stars.tag, state.stars.COcore_mass)
    hydro.set_particle_corem(state.stars.tag, state.stars.core_mass)
    hydro.set_particle_radius(state.stars.tag, state.stars.radius)
    hydro.set_particle_stype(state.stars.tag, state.stars.stellar_type.value_in(units.stellar_type))
    
    return se_dt


def compute_dmdt_vterm(prev_mass, se_temp, se_radius, se_mass, se_lum, dt, massloss_method=None):
    """
    Note: prev_mass = mass before dt update, NOT the ZAMS mass
    """
    if massloss_method == 'seba':

        dm_dt = (prev_mass - se_mass)/dt
        # Since we are using less certain mass loss rates anyway, just use velocity from Leitherer et al. 1992.
        vterm = 10**(1.23 - 0.30 * np.log10(se_lum.value_in(units.LSun))
                + 0.55*np.log10(se_mass.value_in(units.MSun))
                + 0.64*np.log10(se_temp.value_in(units.K))) | units.km/units.s
        
    elif massloss_method == 'seba_puls':
        # Mass loss rates from SeBa with velocities from Kudritzki & Puls winds
        # Added by CCC, 27/11/2024
        dm_dt = (prev_mass - se_mass)/dt
        star_wind   = PulsStellarWind(se_temp, prev_mass, se_lum, se_radius)
        vterm = star_wind.vterm

    # Note that Leitherer and Puls calculations use the old mass

    elif massloss_method == 'leit':
        # Leitherer et. al. 1992.
        dm_dt = 10**(-24.06 + 2.45 * np.log10(se_lum.value_in(units.LSun))
                    -1.10*np.log10(prev_mass.value_in(units.MSun))
                    + 1.31*np.log10(se_temp.value_in(units.K))) | units.MSun/units.yr
        vterm = 10**(1.23 - 0.30 * np.log10(se_lum.value_in(units.LSun))
                    + 0.55*np.log10(prev_mass.value_in(units.MSun))
                    + 0.64*np.log10(se_temp.value_in(units.K))) | units.km/units.s

    elif massloss_method == 'puls':
        # Kudritzki and Puls winds, see Kudritzki & Puls 2000, Markova & Puls 2004, 2008 and Vink 2000
        star_wind   = PulsStellarWind(se_temp, prev_mass, se_lum, se_radius)
        dm_dt = star_wind.dm_dt
        vterm = star_wind.vterm

    else:
        raise Exception("Invalid stellar mass loss method")

    return dm_dt, vterm

def compute_vterm_binary(inj_mass, E_bind):
    """
    Wind method for CE ejection and unstable mass transfer.
    Uses the unbiding energy of the envelope to calculate the ejecta velocity.
    """
    vterm = (2*E_bind/inj_mass)**(1./2)
    return vterm


def compute_eion_nion_sigh(se_mass, se_temp, se_radius):
    """Calculate the average ionizing photon energy based on the blackbody curve."""

    flux = ionizing_photon_flux(se_mass, se_radius, se_temp)

    l_min = 1e-7  # min wavelength, something really small.
    l_max = h*c/E_lyc  # wavelength of 13.6 eV photons, 9.116e-6 cm

    # First integrate the power from the BB curve at this stars temp.
    [power, err] = quad(lum_wl_cs, l_min, l_max, args=(l_max, se_temp.value_in(units.K)))
    # Now integrate to find the number of photons.
    [per_ph, err] = quad(lum_wl_cs_per_ph, l_min, l_max, args=(l_max, se_temp.value_in(units.K)))
    
    avg_E = power/per_ph / E_ev
    # Calculate the average frequency of an ionizing photon for this star
    avg_nu = avg_E*E_ev/h
    # Cross section calculation
    # Make sure you convert energy back to ergs if you
    # use it to calculate the frequency!
    sig = sig0*(h*avg_nu/E_lyc)**(-3)

    eion = (avg_E | units.eV) - (13.6 |units.eV)
    # Calculate total number of photons from stellar surface with stellar radius.
    # Since flux is interpolated from OSTAR2002, no extra factor of pi needed.
    nion = (flux*4*np.pi*se_radius**2).as_quantity_in(units.s**-1)
    sigh = sig | units.cm**2

    return eion, nion, sigh


def compute_epe_npe(se_temp, se_radius):
    """Calculate photoelectric heating parameters"""

    l_min_dust = h*c/E_lyc # wavelength at 13.6 eV
    l_max_dust = h*c/(5.6*E_ev) # wavelength at 5.6 eV

    # First integrate the power from the BB curve at this stars temp.
    [power, err] = quad(lum_wl, l_min_dust, l_max_dust, args=(se_temp.value_in(units.K)))
    # Now integrate to find the number of photons.
    [per_ph, err] = quad(lum_wl_per_ph, l_min_dust, l_max_dust, args=(se_temp.value_in(units.K)))

    avg_E = power/per_ph / E_ev

    # actual average energy of the photons WITH the ionizing potential still in there!
    epe = avg_E | units.eV # should be around 8 eV
    # Calculate total number of photons from stellar surface with stellar radius.
    # Extra factor of pi from solid angle integration of blackbody curve. 
    npe = (np.pi*(per_ph | units.cm**-2*units.s**-1)*4*np.pi*se_radius**2).as_quantity_in(units.s**-1)

    return epe, npe


def went_supernova(stellar_type):
    """                                                                                               
    Determines whether a SeBa star has went supernova or not. Types 13-15 are neutron star, black hole,                                                                                                   
    and disintegrated. This function returns an array or scalar based on input type.                                                                                                      
    """
    types = stellar_type.value_in(units.stellar_type)
    return (types >= 13) & (types <= 15)

def went_supernova_from_kick(star, stars):
    """                                                                                               
    Determines whether a SeBa star has undergone supernova or not, based on whether or not the
    natal kick is set for the star. Used to avoid superfluous SNe from e.g. CEE with BE.
    """
    kick_set = False
    if 'natal_kick_x' in stars.get_attribute_names_defined_in_store():
        if star.natal_kick_x != (0 | units.km/units.s):
            kick_set = True
    return kick_set

# Used to check for common envelope or unstable mass transfer
# This type of mass loss is triggered if >1% of the stars' mass was lost
# in one step - CCC 12/09/2024
def lost_envelope(new_mass, old_mass):
    dm = (old_mass - new_mass)/old_mass
    return dm >= 0.01

# Look for systems with high mass loss rates (> 1e-3 MSun/yr)
# See https://www.aanda.org/articles/aa/pdf/2001/14/aah2347.pdf 
# and https://articles.adsabs.harvard.edu/pdf/1986A%26A...168..111V
# for comparison to O stars and WR stars - CCC 08/12/2024
def high_dm_dt(new_mass, old_mass, dt):
    dm_dt = (old_mass - new_mass)/dt
    return dm_dt > (1.e-2 | units.MSun/units.yr)

# Used to check for stable mass transfer
# If a star accreted at one timestep, do not set the wind properties - CCC 12/09/2024
def accreted_mass(new_mass, old_mass):
    dm = (old_mass - new_mass)/old_mass
    return dm < 0

def lum_wl_cs(l, l_max, T):
    """
    Determine the stellar luminosity at a particular wavelength, temperature and cross section.
    Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
    Note I left out sig0 here b/c we divide this by lum_wl_cs_per_ph that would
    also have sig0 in it.
    """
    # suppress numpy "RuntimeWarning: overflow encountered in exp"
    # which occurs for some, but not all, numpy builds
    # (see, e.g., https://github.com/numpy/numpy/issues/11117).
    # exp(...) overflow to +inf is harmless for l > 1e-50
    if h*c/(l*k*T) > 709.7:  # e^709.7 ~ 1.7e+308 is just below overflow
        return 0
    L = (2*h*c**2/l**5) * (l/l_max)**3 / (np.exp(h*c/(l*k*T)) - 1)
    return L


def lum_wl_cs_per_ph(l, l_max, T):
    """
    Determine the number count of photons at a particular wavelength, temperature and cross section.
    Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
    """
    # suppress numpy "RuntimeWarning: overflow encountered in exp"
    # which occurs for some, but not all, numpy builds
    # (see, e.g., https://github.com/numpy/numpy/issues/11117).
    # exp(...) overflow to +inf is harmless for l > 1e-50
    if h*c/(l*k*T) > 709.7:  # e^709.7 ~ 1.7e+308 is just below overflow
        return 0
    L = (2*h*c**2/l**5) * (l/l_max)**3 / (np.exp(h*c/(l*k*T)) - 1) / (h*c/l)
    return L


def lum_wl(l, T):
    """
    Determine the stellar luminosity at a particular wavelength and temp.
    Uses the standard blackbody curve.
    """
    # suppress numpy "RuntimeWarning: overflow encountered in exp"
    # which occurs for some, but not all, numpy builds
    # (see, e.g., https://github.com/numpy/numpy/issues/11117).
    # exp(...) overflow to +inf is harmless for l > 1e-50
    if h*c/(l*k*T) > 709.7:  # e^709.7 ~ 1.7e+308 is just below overflow
        return 0
    L = (2*h*c**2/l**5) / (np.exp(h*c/(l*k*T)) - 1)
    return L


def lum_wl_per_ph(l, T):
    """
    Determine the number count of photons at a particular wavelength and temp.
    Uses the standard blackbody curve.
    """
    # suppress numpy "RuntimeWarning: overflow encountered in exp"
    # which occurs for some, but not all, numpy builds
    # (see, e.g., https://github.com/numpy/numpy/issues/11117).
    # exp(...) overflow to +inf is harmless for l > 1e-50
    if h*c/(l*k*T) > 709.7:  # e^709.7 ~ 1.7e+308 is just below overflow
        return 0
    L = (2*h*c**2/l**5) / (np.exp(h*c/(l*k*T)) - 1) / (h*c/l)
    return L


class PulsStellarWind(object):
    """Implementation of stellar winds based on Kudritzki and Puls ARAA 2000 and Vink A&A 2000."""

    def __init__(self, teff, mass, lum, radius):

        self.mass     = mass
        self.lum      = lum
        self.teff     = teff
        self.radius   = radius
        self.thom_sig()
        self.thom_Gam()
        self.vesc()
        self.vterm()
        self.dm_dt()

        return

    def thom_sig(self):
        if self.teff < 3e4|units.K:
            self.thom_sig = 0.31 # | units.cm**2.0 / units.g
        elif 3e4|units.K <= self.teff < 3.5e4|units.K:
            self.thom_sig = 0.32 # | units.cm**2.0 / units.g
        else:
            self.thom_sig = 0.33 # | units.cm**2.0 / units.g
        return

    def thom_Gam(self):
        self.thom_Gam = 7.66e-5*self.thom_sig/self.mass.value_in(units.MSun)*self.lum.value_in(units.LSun)
        self.thom_Gam = min(self.thom_Gam, 0.8)
        return

    def vesc(self):
        self.vesc = np.sqrt(2.0*units.constants.G*self.mass*(1-self.thom_Gam)
                            /(self.radius)).as_quantity_in(units.km / units.s)
        return

    def vterm(self):
        if self.thom_Gam >= 1:
            self.vterm = 100 | units.km/units.s # CCC 08/12/2024, set velocity to 100 km/s for test
        else:
            if self.teff <= 1.0e4|units.K:
                self.vterm = self.vesc
            elif 1.0e4|units.K < self.teff < 2.1e4|units.K:
                self.vterm = 1.4*self.vesc
            else:
                self.vterm = 2.65*self.vesc
        return

    def dm_dt(self):
        # Above the bi-stability jump (larger than B1).
        if self.teff > 2.75e4|units.K:
            self.dm_dt = 10**(self.mass_loss1()) | units.MSun / units.yr
        # Below the bi-stability jump (smaller than B1).
        elif self.teff < 2.25e4|units.K:
            self.dm_dt = 10**(self.mass_loss2()) | units.MSun / units.yr
        # Linear interpolation between the two.
        else:
            xp = np.array([2.25e4, 2.75e4])
            fp = np.array([self.mass_loss2(2.25e4), self.mass_loss1(2.75e4)])
            self.dm_dt = 10**(np.interp(self.teff.value_in(units.K), xp, fp))  | units.MSun / units.yr
        return

    # Note we make the temp passable so we can interpolate if we need to
    # and we return a value here for the same reason.

    # Above the bi-stability jump (larger than B1).
    def mass_loss1(self, teff=None):

        if teff is None:
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.697 + 2.194*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.313*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.226*np.log10(self.vterm/self.vesc/2.0) \
                            + 0.933*np.log10(teff/4e4) \
                            - 10.92*np.log10(teff/4e4)**2.0
        return log_dm_dt

    # Below the bi-stability jump (smaller than B1).
    def mass_loss2(self, teff=None):

        if teff is None:
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.688 + 2.210*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.339*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.601*np.log10(self.vterm/self.vesc/2.0) \
                            + 1.07*np.log10(teff/2e4)
        return log_dm_dt
    
    
# Merges stars with delta_r < r_1 + r_2, collisions not handled in current version of petar in amuse
def remove_merged_stars(remove, overwrite, state, hydro, grav, se):
    if remove:
        tprint("... checking for merged stars")
        
        x_mask = np.argsort(state.stars.x.value_in(units.pc))
        x_dist = state.stars.x[x_mask][1:] - state.stars.x[x_mask][:-1]
        r_both = state.stars.radius[x_mask][1:] + state.stars.radius[x_mask][:-1]
        r_mask = np.where(x_dist <= r_both)
        # Check 3D distance for those
        r_dist = ((state.stars.x[x_mask][1:][r_mask] - state.stars.x[x_mask][:-1][r_mask])**2
                 + (state.stars.y[x_mask][1:][r_mask] - state.stars.y[x_mask][:-1][r_mask])**2
                 + (state.stars.z[x_mask][1:][r_mask] - state.stars.z[x_mask][:-1][r_mask])**2)**(1./2)
        idx_w = np.where(r_dist <= r_both[r_mask])[0]
        idx_1 = x_mask[1:][r_mask][idx_w]
        idx_2 = x_mask[:-1][r_mask][idx_w]

        # loop over pairs of stars with identical positions
        if len(idx_w) > 0: # Check if array is empty
            stars_rem = Particles()
            for i in range(len(idx_w)):
                star1_idx = idx_1[i]
                star2_idx = idx_2[i]
                se.particles[star1_idx].merge_with_other_star(se.particles[star2_idx])
                # Save tag of star it merged with
                state.stars[star2_idx].merged_with = state.stars[star1_idx].tag
                # Save merged time
                state.stars[star2_idx].merger_time = hydro.get_time()
                stars_rem.add_particle(state.stars[star2_idx])

            # hydro requires sorted tags for removal
            # only the stars particle set has a tag attribute.
            t = stars_rem.tag
            t = np.sort(np.array(t).flatten())
            tprint("Removing ", len(t), "merged star(s)")
            # Remove from hydro
            hydro.remove_particles(t)
            # Remove from SE
            se.particles.remove_particles(stars_rem)
            # Synchronize to state and copy mass
            se.particles.synchronize_to(state.stars)
            state.se_to_stars.copy_attributes(["mass"])
            # Remove and re-add to grav
            state.stars.synchronize_to(grav.particles)
            state.stars_to_grav.copy_attributes(["mass"])
            if len(grav.particles) != len(state.stars):
                # See this issue: https://github.com/amusecode/amuse/issues/518
                tprint('... forced to re-sync grav from stars')
                grav.particles = Particles()
                grav.particles.add_particles(state.stars)
            state.out_merged_stars(stars_rem, overwrite)
               
        else:
            pass

    
if __name__ == '__main__':
    pass
