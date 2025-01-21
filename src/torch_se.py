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

from amuse.units import constants, units #constants added for tests, CCC 28/04/2023

from ionizingflux import ionizing_photon_flux
from torch_stdout import tprint

# CCC 28/04/2023, temporary, use same convention as for single stars later on
from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particles

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
sigDust = 1e-21 | units.cm**2.0 # Cross section for dust from Draine 2011
# TODO should sigDust be a user-controlled parameter? -AT, 2019oct14


def binary_evolution(time, dt, state, hydro, se,
                     with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                     massloss_method=None, min_feedback_mass=None, CE_method='wind', CE_alpha=1):
    """
    NOTE: time = target time to evolve TO, including the dt already.
    Chosen to follow AMUSE worker convention.
    """
    assert massloss_method is not None
    assert min_feedback_mass is not None

    # We call SeBa on indiv stars, but get/set hydro star props in bulk.

    # Always recompute star's age from hydro time and particle creation time.
    # Don't attach star age to particle.  Why?  (1) Repeated increment of star
    # age at each bridge step would introduce error.  (2) Multiple ways to
    # query star age may not agree exactly.

    # Age not used explicitly for SE - CCC 02/08/2024
    state.stars.age = (time - dt) - hydro.get_particle_creation_time(state.stars.tag)
    
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

    # Pass information to SE
    # Check for new systems, important for binaries - CCC 02/08/2024
    state.binaries.synchronize_to(se.binaries)
    #tprint('SeBa binaries:', se.binaries.binary_type)
    # Now pass attributes to binaries - CCC 03/08/2024
    for _attribute in state.binaries.get_attribute_names_defined_in_store():
        if _attribute in se.binaries.get_attribute_names_defined_in_store():
            state.binaries_to_se.copy_attributes([_attribute])
    # Print timesteps
    #tprint('Timesteps:', se.particles.time_step, se.binaries.time_step)
    # Evolve model
    se.evolve_model(time) #Time attached to se.particles.age, which is the "simulation" time
    #tprint('SeBa binaries:', se.binaries.binary_type)
    
    # Pass information back to stars after end of SE loop - CCC 02/08/2024
    for _attribute in se.particles.get_attribute_names_defined_in_store():
        if _attribute in state.stars.get_attribute_names_defined_in_store():
            state.se_to_stars.copy_attributes([_attribute])

    # Turn off wind if mass increased - CCC 11/09/2024 
    # Keep a list of stars for which feedback was calculated - CCC 25/10/2024
    evolved_stars = Particles()
    # List primaries and companions - CCC 25/10/2024
    primaries = Particles()
    companions = Particles()
    for binary in state.binaries:
        primaries.add_particle(binary.child1)
        companions.add_particle(binary.child2)

    
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
        
        if (s in companions): # CCC 25/10/24
            in_binary = True
            # Find the binary
            j = np.where(companions.tag == s.tag)[0]
            b = state.binaries[j]
            # and the other star
            k = np.where(state.stars.tag == primaries[j].tag)[0][0]
            t = state.stars[k]

        if went_supernova(s.stellar_type):
            continue

        if s.initial_mass >= min_feedback_mass:

            if with_sn and went_supernova(s.stellar_type): # Do not check if in binary for SN - CCC 29/10/2024

                inj_mass = old_mass[i] - s.mass  # minus stellar remnant's mass
                
                if inj_mass > 15.0|units.MSun:
                    # expected upper limit for SeBa tracks; see
                    # https://groups.google.com/forum/#!topic/torch-users/rWJd6l_mRBg/discussion
                    tprint("... flooring SN inj_mass {} MSun to 15 MSun".format(inj_mass.value_in(units.MSun)))
                    inj_mass = 15.0|units.MSun

                # inject energy and mass onto grid
                _tmp = hydro.energy_injection(1e51|units.erg, -1.0, inj_mass.in_(units.g), s.x, s.y, s.z)
                se_dt = min(se_dt, _tmp)
                tprint("... SN x={}, y={}, z={}, inj_mass={}, tag={}".format(s.x, s.y, s.z, inj_mass.value_in(units.MSun), s.tag))

                # implicitly zeros out feedback properties by not setting
                
                # Set star to evolved after SN - CCC 29/10/2024
                evolved_stars.add_particle(s)
                
            else: # Check if in binary to determine the feedback mechanisms - CCC 29/10/2024
                
                if in_binary:
                    
                    _dE = (units.constants.G / 2) * (((s.mass * t.mass / se.binaries[j].semi_major_axis) - old_mass[i]*old_mass[k] / b.semi_major_axis))[0] # Try here - CCC 08/12/2024
                    tprint('Change in binding energy:', "{0:.1e}".format(_dE.value_in(units.erg)), 'erg')
                    inj_mass = old_mass[i] + old_mass[k] - s.mass - t.mass
                    # If _dE > 0, mass transfer or CE; if _dE < 0, wind mass loss
            
                    if _dE > (0 | units.erg): 
                
                        tprint('... Do mass transfer')
                    
                        tprint('Injected mass:', "{0:.2f}".format(inj_mass.value_in(units.MSun)), 'MSun')
                
                        # https://www.aanda.org/articles/aa/full_html/2021/04/aa40442-21/aa40442-21.html
                        # CCC 12/09/2024, 20/06/2023
                        E_bind = CE_alpha * _dE
                        tprint('Ejecta energy:', "{0:.1e}".format(E_bind.value_in(units.erg)), 'erg')
                        
                        # If there is no energy from the envelope ejection, do wind mass loss
                        # for the donor star instead
                        
                        if CE_method=='wind':
                            
                            _vterm = compute_vterm_binary(inj_mass, E_bind)
                            tprint('Ejecta velocity from change in energy', "{0:.1e}".format(_vterm.value_in(units.km/units.s)), 'km/s')
                            
                            if (old_mass[i] - s.mass) > (old_mass[k] - t.mass):
                                # If donor is star s
                                _, _vterm_wind = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity,
                                                                    dt, massloss_method=massloss_method)
                                dm_dt[i] = inj_mass/dt
                                vterm[i] = min([_vterm, _vterm_wind]) # Cap the ejecta velocity at the wind velocity
                                tprint('Ejecta velocity', "{0:.1e}".format(vterm[i].value_in(units.km/units.s)), 'km/s')
                            else:
                                # If donor is star t
                                _, _vterm_wind = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity,
                                                                    dt, massloss_method=massloss_method)
                                dm_dt[k] = inj_mass/dt
                                vterm[k] = min([_vterm, _vterm_wind]) # Cap the ejecta velocity at the wind velocity
                                tprint('Ejecta velocity', "{0:.1e}".format(vterm[k].value_in(units.km/units.s)), 'km/s')
                        
                        elif CE_method=='SN':
                
                            _tmp = hydro.energy_injection(E_bind, -1.0, inj_mass.in_(units.g), s.x, s.y, s.z)
                
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
                            sigpe[i] = sigDust  # TODO magic constant -AT 2019Oct14
                            _tmp = compute_epe_npe(t.temperature, t.radius)
                            epe[k] = _tmp[0]
                            npe[k] = _tmp[1]
                            sigpe[k] = sigDust  # TODO magic constant -AT 2019Oct14
                        # Do not set wind properties for stars that have lost mass due to CE
                        
                        # Set star to evolved after uMT/CE - CCC 22/11/2024
                        evolved_stars.add_particle(s)
                        evolved_stars.add_particle(t)
                        
                    else:
                    
                        if accreted_mass(t.mass, old_mass[k]): #If accreted but still lost energy over dt
                            
                            _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                  massloss_method=massloss_method)
                            dm_dt[i] = inj_mass / dt
                            vterm[i] = _tmp[1]
                            tprint('Star 1 wind dm/dt:', dm_dt[i])
                            tprint('Star 1 wind velocity:', vterm[i])
                            
                        elif accreted_mass(s.mass, old_mass[i]):
                            
                            _tmp = compute_dmdt_vterm(old_mass[k], t.temperature, t.radius, t.mass, t.luminosity, dt,
                                                  massloss_method=massloss_method)
                            dm_dt[k] = inj_mass / dt
                            vterm[k] = _tmp[1]
                            tprint('Star 2 wind dm/dt:', dm_dt[k])
                            tprint('Star 2 wind velocity:', vterm[k])
                            
                        else: # Normal feedback if the binary does not interact
                            
                            if with_winds:
                                _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                      massloss_method=massloss_method)
                                dm_dt[i] = _tmp[0]
                                vterm[i] = _tmp[1]
                                tprint('Star 1 wind dm/dt:', dm_dt[i])
                                tprint('Star 1 wind velocity:', vterm[i])
                                _tmp = compute_dmdt_vterm(old_mass[k], t.temperature, t.radius, t.mass, t.luminosity, dt,
                                                      massloss_method=massloss_method)
                                dm_dt[k] = _tmp[0]
                                vterm[k] = _tmp[1]
                                tprint('Star 2 wind dm/dt:', dm_dt[k])
                                tprint('Star 2 wind velocity:', vterm[k])
                        
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
                            sigpe[i] = sigDust  # TODO magic constant -AT 2019Oct14
                            _tmp = compute_epe_npe(t.temperature, t.radius)
                            epe[k] = _tmp[0]
                            npe[k] = _tmp[1]
                            sigpe[k] = sigDust  # TODO magic constant -AT 2019Oct14
                    
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
                        sigpe[i] = sigDust  # TODO magic constant -AT 2019Oct14
                    if with_winds:
                        _tmp = compute_dmdt_vterm(old_mass[i], s.temperature, s.radius, s.mass, s.luminosity, dt,
                                                  massloss_method=massloss_method)
                        dm_dt[i] = _tmp[0]
                        vterm[i] = _tmp[1]
                    
                    # Set star to evolved after feedback - CCC 29/10/2024
                    evolved_stars.add_particle(s)           
            
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

    return se_dt


def stellar_evolution(time, dt, state, hydro, se,
    with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
                      massloss_method=None, min_feedback_mass=None):
    """
    NOTE: time = target time to evolve TO, including the dt already.
    Chosen to follow AMUSE worker convention.
    """
    assert massloss_method is not None
    assert min_feedback_mass is not None
    
    # We call SeBa on all stars at once but loop through stars for feedback
    # TO DO: Edit to loop only through feedback stars (see bitbucket) - CCC 04/11/2023

    # Always recompute star's age from hydro time and particle creation time.
    # Don't attach star age to particle.  Why?  (1) Repeated increment of star
    # age at each bridge step would introduce error.  (2) Multiple ways to
    # query star age may not agree exactly.
    state.stars.age = (time - dt) - hydro.get_particle_creation_time(state.stars.tag)
    
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
    
    # Structure changed to use evolve_model to evolve all stars at the same time
    # This allows us to restart from evolved stars and use the same structure for
    # binary evolution - CCC 04/11/2023
    state.stars_to_se.copy() # CCC 25/07/2024
    se.evolve_model(time)
    state.se_to_stars.copy()

    for i, s in enumerate(state.stars):

        if went_supernova(s.stellar_type):
            continue

        if s.mass >= min_feedback_mass:

            if with_sn and went_supernova(s.stellar_type):

                inj_mass = old_mass[i] - s.mass  # minus stellar remnant's mass
                
                if inj_mass > 15.0|units.MSun:
                    # expected upper limit for SeBa tracks; see
                    # https://groups.google.com/forum/#!topic/torch-users/rWJd6l_mRBg/discussion
                    tprint("... flooring SN inj_mass {} MSun to 15 MSun".format(inj_mass.value_in(units.MSun)))
                    inj_mass = 15.0|units.MSun

                # inject energy and mass onto grid
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
                    sigpe[i] = sigDust  # TODO magic constant -AT 2019Oct14
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
        tprint('Wind properties:', star_wind.thom_Gam, star_wind.vesc, star_wind.dm_dt)
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
    return 13 <= stellar_type.value_in(units.stellar_type) <= 15

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
        tprint('Gamma =', self.thom_Gam)
        self.thom_Gam = min(self.thom_Gam, 0.5)
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

    
if __name__ == '__main__':
    pass
