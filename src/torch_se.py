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

sigSB = 5.6704e-5 | (units.g/((units.s)**3 * (units.K)**4)) # Stefan-Boltzmann constant, g s^-3 K^-4
sig0  = 6.304e-18 # Photoionization cross section at threshold for hydrogen
E_ev  = 1.60222497096e-12 # energy of 1 eV in erg
E_lyc = 13.6*E_ev  # 13.6 eV
# Cross section for dust per hydrogen atom.
# Value = tau / N_H where tau = gamma * Av (Draine and Bertoli 96)
# Av = N_H,tot / (1.87e21 cm^2) (Bohlin et al 78)
# gamma = 2.5 (Bergin et al 2004)
sigDust = 1e-21 | units.cm**2.0 # Cross section for dust from Draine 2011
# TODO should sigDust be a user-controlled parameter? -AT, 2019oct14


def binary_evolution(time, dt, state, hydro, worker,
    with_lyc=True, with_pe_heat=True, with_winds=True, with_sn=True,
    massloss_method=None, min_feedback_mass=None):
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
    state.stars.age  = time - hydro.get_particle_creation_time(state.stars.tag)
    
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
    
    worker.evolve_model(time)
    state.se_to_stars.copy()
    state.se_to_binaries.copy()

    for i, s in enumerate(state.stars):

        if went_supernova(s.stellar_type):
            continue

        if s.initial_mass >= min_feedback_mass:

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

    return se_dt


def stellar_evolution(time, dt, state, hydro, worker,
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
    state.stars.age = time - hydro.get_particle_creation_time(state.stars.tag)
    
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
    worker.evolve_model(time)
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
        return

    def vesc(self):
        self.vesc = np.sqrt(2.0*units.constants.G*self.mass*(1-self.thom_Gam)
                            /(self.radius)).as_quantity_in(units.km / units.s)
        return

    def vterm(self):
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
