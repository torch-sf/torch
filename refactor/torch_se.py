"""
Torch code to do stellar evolution, JUST on AMUSE side.

Includes subroutines and stuff to...

* Compute information about star luminosity (energy, photon count, etc) at
  requested wavelength and temperature

* implement stellar wind

Joshua Wall, Drexel University
"""

from __future__ import division

import numpy as np

#from amuse.lab import *
from amuse.units import units

def do_stellar_evolution(stars, hydro):

    # must add dt, else new stars will evolve for 0.0 sec and things will break
    star_age  = hydro.get_time() + dt - hydro.get_particle_creation_time(stars.tag)
    star_mass = hydro.get_particle_mass(stars.tag)
    star_type = stars.stellar_type

    # for winds and ionizing radiation
    dm_dt     = np.zeros(num_particles) | units.g / units.s
    vterm     = np.zeros(num_particles) | units.cm / units.s
    nphot     = np.zeros(num_particles) | units.s**-1.0
    eion      = np.zeros(num_particles) | units.erg
    sigh      = np.zeros(num_particles) | units.cm**2.0
    npe       = np.zeros(num_particles) | units.s**-1.0
    epe       = np.zeros(num_particles) | units.erg
    sigpe     = np.zeros(num_particles) | units.cm**2.0

    part_inds = []

    print "Doing stellar evolution."

    for part in range(num_particles):

        if star_age[part].value_in(units.yr) < 1.0:
            star_age[part] = 1.0 | units.yr

        # Do stellar evolution unless I already went SN, in which case skip me.
        if (13 <=  stars.stellar_type[part].value_in(units.stellar_type) <= 15):
            print "Skipping this star that already went SN, current stellar type =", stars.stellar_type[part]
            continue

        # SE code accepts initial mass, not the current mass
        st_time, st_mass, star_radius, star_lum, star_temp, st_evol_time, st_type = se.evolve_star(stars.initial_mass[part], star_age[part], 0.02)
        star_type[part] = st_type

        if (with_massloss and (massloss_method == 'seba' or st_mass.value_in(units.MSun) < min_mass.value_in(units.MSun))):
            dm_dt[part] = ((stars.mass[part]-st_mass)/st_time).in_(units.g / units.s)
            # If below the mass cutoff for feedback, no wind present.
            if (st_mass.value_in(units.MSun) < min_mass.value_in(units.MSun)):
                vterm[part] = 0.0 | units.km/units.s
            # Since we are using less certain mass loss rates anyway, just use velocity from Leitherer et. al. 1992.
            else:
                vterm[part] = 10**(1.23 - 0.30 * np.log10(star_lum.value_in(units.LSun))
                        + 0.55*np.log10(st_mass.value_in(units.MSun))
                        + 0.64*np.log10(star_temp.value_in(units.K))) | units.km/units.s
        # Shouldn't Lietherer and Puls calculations use the old mass (stars[part].mass)?
        elif (with_massloss and massloss_method == 'leit'):
            # From Leitherer et. al. 1992.
            dm_dt[part] = 10**(-24.06 + 2.45 * np.log10(star_lum.value_in(units.LSun))
                        -1.10*np.log10(stars[part].mass.value_in(units.MSun))
                        + 1.31*np.log10(star_temp.value_in(units.K))) | units.MSun/units.yr
            vterm[part] = 10**(1.23 - 0.30 * np.log10(star_lum.value_in(units.LSun))
                        + 0.55*np.log10(stars[part].mass.value_in(units.MSun))
                        + 0.64*np.log10(star_temp.value_in(units.K))) | units.km/units.s
        elif (with_massloss and (massloss_method == 'puls' or st_mass.value_in(units.MSun) >= min_mass.value_in(units.MSun))):
            # Kudritzki and Puls winds, see Kudritzki & Puls 2000, Markova & Puls 2004, 2008 and Vink 2000
            star_wind   = StellarWind(star_temp, stars[part].mass, star_lum, star_radius)
            dm_dt[part] = star_wind.dm_dt.as_quantity_in(units.g / units.s)
            vterm[part] = star_wind.vterm.as_quantity_in(units.cm / units.s)
        else:
            dm_dt[part] = 0.0 | units.g / units.s
            vterm[part] = 0.0 | units.cm / units.s

        # If with energy injection, check to see if anything went supernova. If so, inject 10^51 ergs of
        # energy into the grid
        if (with_sn and 13 <= st_type.value_in(units.stellar_type) <= 15):

            print "A star just went SN on you. Should be calling that SN code now!"
            print "Going supernova at", stars.x[part], stars.y[part], stars.z[part]

            inj_x = stars.x[part]
            inj_y = stars.y[part]
            inj_z = stars.z[part]
            tot_e = 1e51 | units.erg
            fracKin = -1.0
            # injected mass is current mass minus remnant mass.
            inj_mass = (star_mass[part] - st_mass).in_(units.g)
            # This should never be more than 10 solar masses though.
            if (inj_mass.value_in(units.MSun) > 10.0):
                print "[bridge:SN]: WARNING! SN MASS > 10 MSun!"
                inj_mass = 10.0 | units.MSun

            sn_dt = hydro.energy_injection(tot_e, fracKin, inj_mass, inj_x, inj_y, inj_z)
            dt =  min(dt.value_in(units.s), sn_dt.value_in(units.s)) | units.s
            print "Timestep after SN is =", dt
            print "Now evolving until t =", t + dt

            # Set proper mass for remnant (SeBa does fine for this).
            star_mass[part] = st_mass
            # Set proper remnant stellar type so that we don't get any feedback from remnants.
            stars.stellar_type[part] = st_type
            # Switch off all feedback for this star.
            nphot[part] = 0.0 | units.s**-1
            eion[part]  = 0.0 | units.erg
            if (use_radiation):
                hydro.set_particle_nion(stars.tag[part], nphot[part])
                hydro.set_particle_eion(stars.tag[part], eion[part])
            npe[part] = 0.0 | units.s**-1
            epe[part] = 0.0 | units.eV
            if (pe_heat):
                hydro.set_particle_npep(stars.tag[part], npe[part])
                hydro.set_particle_epep(stars.tag[part], epe[part])
            # Winds don't matter for a SN star and we don't want to lose any more mass.
            dm_dt[part] = 0.0 | units.g / units.s
            vterm[part] = 0.0 | units.cm / units.s
            if (with_winds):
                hydro.set_particle_wind_mass(stars.tag[part], dm_dt[part])
                hydro.set_particle_wind_vel(stars.tag[part], vterm[part])
            # Do nothing else with this star, just jump to the next one.
            continue

        if (star_mass[part].in_(units.MSun) >= min_mass
            and not (13 <= st_type.value_in(units.stellar_type) <= 15)):

            part_inds.append(part)

            if (use_radiation):

                flux = ion.ionizing_photon_flux(st_mass, star_radius, star_temp)

                # Calculate the average ionizing photon energy based on the blackbody curve.

                # First integrate the power from the BB curve at this stars temp.
                # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                [power, err] = quad(lum_wl_cs, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                # Now integrate to find the number of photons.
                [per_ph, err] = quad(lum_wl_cs_per_ph, l_min, l_max, args=(l_max, star_temp.value_in(units.K)))
                avg_E = power/per_ph / E_ev
                # Calculate the average frequency of an ionizing photon for this star
                avg_nu = avg_E*E_ev/h
                # Cross section calculation
                # Make sure you convert energy back to ergs if you
                # use it to calculate the frequency!
                sig = sig0*(avg_nu/nu_min)**(-3.0)
                eion[part] = (avg_E | units.eV) - (13.6 |units.eV) #2.0 | units.eV #6.0 | units.eV
                sigh[part] = sig | units.cm**2.0 #6.3e-18 | units.cm**2.0
                nphot[part] = (flux*4*np.pi*star_radius**2.0).as_quantity_in(units.s**-1.0) #5e48 | units.s**(-1.0)

                if (pe_heat):

                    # First integrate the power from the BB curve at this stars temp.
                    # l_min=1e-7 (small enough), min wavelength, l_max=9.116e-6 cm, wavelength of 13.6 eV photons.
                    l_min_dust = h*c / E_min # wavelength at 13.6 eV
                    l_max_dust = h*c / (5.6*E_ev) # wavelength at 5.6 eV
                    [power, err] = quad(lum_wl, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))
                    # Now integrate to find the number of photons.
                    [per_ph, err] = quad(lum_wl_per_ph, l_min_dust, l_max_dust, args=(l_max_dust, star_temp.value_in(units.K)))

                    avg_E = power/per_ph / E_ev
                    # Cross section calculation
                    # We assume constant cross section for dust per hydrogen atom.
                    # Value = tau / N_H where tau = gamma * Av (Draine and Bertoli 96)
                    # Av = N_H,tot / (1.87e21 cm^2) (Bohlin et al 78)
                    # gamma = 2.5 (Bergin et al 2004)
                    sigpe[part] = sigDust
                    # Eion is the actual average energy of the photons WITH the ionizing potential still in there!
                    epe[part] = avg_E | units.eV # should be around 8 eV
                    # Calculate total number of photons from stellar surface with stellar radius.
                    npe[part] = ((per_ph | units.cm**-2*units.s**-1)*4*np.pi*star_radius**2.0).as_quantity_in(units.s**-1.0) #5e48 | units.s**(-1.0)

        if ((dm_dt[part]*dt).value_in(units.MSun) > 0.0):
            if (st_type.value_in(units.stellar_type) == 1):
                star_mass[part] = ((stars.mass[part] - dm_dt[part]*dt).value_in(units.MSun)) | units.MSun
            # Note other evolutionary things besides winds could have reduced the stars mass.
            else:
                star_mass[part] = min(st_mass.value_in(units.MSun), (stars.mass[part] - dm_dt[part]*dt).value_in(units.MSun)) | units.MSun
        else:
            star_mass[part] = st_mass

    # push stars mass -> multiples leaf, grav root
    # imitating bridge kick update of mult_grav,grav - AT 2019 Jul 01
    if (with_multiples):
        # QUESTION do we need mult_grav.channel_from_code_to_memory.copy() ?
        # I'm not sure what that does.
        # -AT 2019 Jul 01
        for st, st_mass in zip(stars, star_mass):
            for root, tree in mult_grav.root_to_tree.iteritems():
                leaves = tree.get_leafs_subset()
                if st in leaves:
                    st.as_particle_in_set(leaves).mass = st_mass
        # If the leaf mass changed, so does the com particle's properties.
        update_roots_from_leaves(mult_grav, grav) # If works, move outside se loop!

    # Are there any massive stars?

    if part_inds:  # empty list evaluates to False

        if use_radiation:

            print  stars.tag[part_inds]
            print "Stellar Mass and N photons=", star_mass[part_inds], nphot[part_inds]
            print "Eion (eV), SigH=", eion[part_inds], sigh[part_inds]
            hydro.set_particle_nion(stars.tag[part_inds], nphot[part_inds])
            hydro.set_particle_eion(stars.tag[part_inds], eion[part_inds].as_quantity_in(units.erg))
            hydro.set_particle_sigh(stars.tag[part_inds], sigh[part_inds])


        if pe_heat:
            print "Npe photons=", npe[part_inds]
            print "Eion PE (eV), SigD=", epe[part_inds], sigpe[part_inds]
            hydro.set_particle_npep(stars.tag[part_inds], npe[part_inds])
            # Set average energy of PE photon
            hydro.set_particle_epep(stars.tag[part_inds], epe[part_inds].as_quantity_in(units.erg))
            # Set cross section of dust to PE photons.
            hydro.set_particle_sigd(stars.tag[part_inds], sigpe[part_inds])


        if with_winds:
            hydro.set_particle_wind_mass(stars.tag[part_inds], dm_dt[part_inds])
            hydro.set_particle_wind_vel(stars.tag[part_inds], vterm[part_inds])

    # Remove any mass loss due to winds and update to this
    # mass. Note this assumes steps are relatively small
    # in the mass loss rate of stars, so that gravity
    # can use the mass after all the wind mass loss
    # has occcured. Otherwise we'd have to average
    # mass loss and keep up with old and new masses and
    # it just gets ugly.
    stars.mass = star_mass # - dm_dt*dt
    stars.age  = stars.age + dt
    stars.stellar_type = star_type


def lum_wl_cs(l, l_max, T):
    """
    Determine the stellar luminosity at a particular wavelength, temperature and cross section.
    Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
    Note I left out sig0 here b/c we divide this by lum_wl_cs_per_ph that would
    also have sig0 in it.
    """
    h = 6.6261e-27 # Planck's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzmann constant
    L = (2.0*h*(c**2.0)/(l**5.0)) * (l/l_max)**3.0 / (np.exp(h*c/(l*k*T)) - 1.0)
    return L


def lum_wl_cs_per_ph(l, l_max, T):
    """
    Determine the number count of photons at a particular wavelength, temperature and cross section.
    Uses the standard blackbody curve and incorporates the cross section as a function of wavelength.
    """
    h = 6.6261e-27 # Planck's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzmann constant
    L = (2.0*h*(c**2.0)/(l**5.0)) * (l/l_max)**3.0 / (np.exp(h*c/(l*k*T)) - 1.0) / (h*c/l)
    return L


def lum_wl(l, l_max, T):
    """
    Determine the stellar luminosity at a particular wavelength and temp.
    Uses the standard blackbody curve.
    """
    h = 6.6261e-27 # Planck's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzmann constant
    L = (2.0*h*(c**2.0)/(l**5.0)) / (np.exp(h*c/(l*k*T)) - 1.0)
    return L


def lum_wl_per_ph(l, l_max, T):
    """
    Determine the number count of photons at a particular wavelength and temp.
    Uses the standard blackbody curve.
    """
    h = 6.6261e-27 # Planck's constant
    c = 2.9979e10  # Speed of light
    k = 1.3807e-16 # Boltzmann constant
    L = (2.0*h*(c**2.0)/(l**5.0)) / (np.exp(h*c/(l*k*T)) - 1.0) / (h*c/l)
    return L


class StellarWind(object):
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
        if (self.teff.value_in(units.K) < 3e4):
            self.thom_sig = 0.31 # | units.cm**2.0 / units.g
        elif (3e4 <= self.teff.value_in(units.K) < 3.5e4):
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
        if (self.teff.value_in(units.K) <= 1.0e4):
            self.vterm = (self.vesc)
        elif (1.0e4 < self.teff.value_in(units.K) < 2.1e4):
            self.vterm = (1.4*self.vesc)
        else:
            self.vterm = (2.65*self.vesc)
        return

    def dm_dt(self):
        # Above the bi-stability jump (larger than B1).
        if (self.teff.value_in(units.K) > 2.75e4):
            self.dm_dt = 10**(self.mass_loss1()) | units.MSun / units.yr
        # Below the bi-stability jump (smaller than B1).
        elif (self.teff.value_in(units.K) < 2.25e4):
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

        if (teff is None):
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.697 + 2.194*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.313*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.226*np.log10(self.vterm/self.vesc/2.0) \
                            + 0.933*np.log10(teff/4e4) \
                            - 10.92*np.log10(teff/4e4)**2.0
        return log_dm_dt

    # Below the bi-stability jump (smaller than B1).
    def mass_loss2(self, teff=None):

        if (teff is None):
            teff = self.teff.value_in(units.K)

        log_dm_dt  = -6.688 + 2.210*np.log10(self.lum.value_in(units.LSun)/1e5) \
                            - 1.339*np.log10(self.mass.value_in(units.MSun)/30.0) \
                            - 1.601*np.log10(self.vterm/self.vesc/2.0) \
                            + 1.07*np.log10(teff/2e4)
        return log_dm_dt


if __name__ == '__main__':
    pass
