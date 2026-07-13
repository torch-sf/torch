### Metallicity Review ###
M-MML starting 23 Jun 2026, then edited by CC-C and M-MML on 29 Jun 2026. _Note that this was created based off of main, and the notes here will be moved to branch created from develop._

# DriverMain/Drive_sourceTerms.F90 #
This calls `Heat`, `Heatexchange`, `Cool`, `Ionize` in this order; `Heatexchange` is vanilla FLASH; `Heat` calls the `heating_and_cooling` subroutine so `Cool` is not used.

## RadHeat ##
active version I believe to be `src/flash/source/physics/sourceTerms/Heat/HeatMain/HeatCool/phenHeat/mol_and_dust/solver/RadHeat.F90`

The header file is `src/flash/source/physics/sourceTerms/Heat/HeatMain/HeatCool/phenHeat/Heat_Interface.f90`

### Cooling ###
Cooling uses dust cooling from Goldsmith (2001; https://scixplorer.org/abs/2001ApJ...557..736G/abstract).  This is implicitly solar.  To replace that, we can look at Dopcke et al. (2011;  https://scixplorer.org/abs/2011ApJ...729L...3D/abstract).  I've touched base with Ralf K, who suggested emailing Simon G., which I've done (24 June)

# heatCool.F90 #
The `cooling` function is called by `dei_dt`, which is called in the subroutine `heating_and_cooling`. The subroutine `cooling` begins at l. 1116. The 0.01 ionization fraction is hardcoded here, on l. 1143; *this needs to be updated*. This returns `emin_out`, which is cooling rate in erg cm-3 s-1, which then gets converted to erg g-1 s-1.

If `useDustCool` is true, then `molecular_cooling` and `dust_cooling` are also called.

_Atomic cooling_ : If temperature and densities within range set in flash.par, call `Radloss` subroutine, which calls `cooling`, which is the atomic cooling. This routine will have to be entirely replaced. 

On l.1525, calling *atomic cooling* from Dalgarno and McCray. Different functions are available for different ionization fractions; this will be replaced by Simon Glover (see notes below). This (combined with the below) will return an output variable used in `dei_dt`.

At low densities, this will be replaced by Simon Glover's updates; at high densities, *we need a set of metallicity-dependent cooling curves -- figure out which are the usual ones*.

_Molecular cooling_: On `heatCool.F90`, l.1211, the subroutine `molecular_cooling` is defined. This calls `cool_dat`, which gets the cooling table from `cool.dat`. This is based on Neufeld et al. 1995. 

This is valid over T=10-2500 K, H_2 densities 1e3-1e10 cm^-3, and assumes equilibrium. *This uses a CR ionization rate of 1e-17 per H2 molecular; in flash.par, we use 2e-17 -- is this by H atom?* This uses mu_mol and gets the cooling rate for the corresponding temperature and density from `cool.dat`.

Only one set of abundances is used, based on Galactic abundances.

At our densities, cooling is dominated by CO. CO is not self-shielding, only shielded by dust -- decreases more quickly than H2. *We should look up what did Enzo use before they introduced Grackle*. 

_Dust cooling_: This is deifined on l. 1305 in `heatCool.F90`. It uses cooling rates from Hollenbach and McKee 1989, as reported in Glover and Clark 2011. This includes several constants from Hollenbach & McKee 1989, as reported in Glover & Clark 2012b. 

# get_cooling_data.F90 #
This calls `cool.dat` and creates a variable called `cool_dat`, which containes temperature, densities, and cooling powers. This is called in `heatCool.F90`.

### Heating ###
Note: _There is a commented out function to fill the guard cells._
This uses constants.h, which is located at `/src/flash/source/Simulation/SimulationMain/StratBox/constants.h`. _Why is this in the StratBox directory?_ In vanilla FLASH, this is in the top-level `/Simulation` directory. There are no differences between the files.

*Heat_data* is a data file; this will likely be part of replacing the data files from Robi.

# Heat.F90 #
This call `RadHeat` from `RadHeat.F90` if (1) `IHP_SPEC` is not used or (2) `rt_heatInRad` or `rt_useRadTrans` is not used. `IHP_SPEC` is not defined; therefore, the `RadHeat` routine is used, although both `rt_heatInRad` and `rt_useRadTrans` are true.

# RadHeat.F90 #
On l. 278, the mean molecular weight mu_mol is hardcoded to 24/11, which is 10% He by number, which corresponds to 40% by weight. This number should be adjusted with metallicity.

# heatCool.F90 #
Note: _There is a commented out `use heatCool` call on l.48. Why?_
*Note also that we are only following along the `implicit` method, as this is the default choice; for other methods, we do not follow the line of references.*

`mu_mol` is used. It is defined in `RadHeat.F90` then passed around in `OdeData` (see `cool_vars`). This appears to be the only definition of `mu_mol` in the FLASH source files, *hardcoded elsewhere?* `mu_mol` is also used to calculate the dust temperature.

*On l.896, dust_gas_ratio = 0.01 is hardcoded.* This is used in f_ext, on l. 1718, which is a local approximation for self-shielding. *This is hardcoded everywhere else*, as the variable is not used anywhere else. f_ext is used in the local FUV flux calculation.

*Note also that the dust-to-gas ratio is used in VETTAM and likely defined separately there*. 

The CR ionization rate is hardcoded to the Milky Way value. See `he_crIonRate`, `he_crIonEnergy`,`he_crIonNH`, `he_crIonExp`; example of use is l. 928. It looks like heating is only CR heating and not background UV heating, which is defined in `cool_vars` as `Gflux`, which is taken from `PEFL`. 

Check W&D heating constants -- are those universal? They are on l. 985.

On l. 1022, the PE heating routine needs to checked. The reference is Bakes/Tielens 1994, then Wolfire 2003. Use the Wolfire paper to check the "magic numbers", and figure out which ones are metallicity dependent (see equation 20). Check for updates by Wolfire or others, *this needs more research as we currently do not have a metallicity-dependent heating rate*. 

On l. 1051, the he_pe_recipe also has several numbers hardcoded, which come from Weingarter & Draine 2001. *Also check those*. The rate should depend on G, T and n_e; the function included in the code is based on eq. 44. This is appropriate between 10-1e4 K and 1e2 K^1/2 cm^3 < G sqrt(T)/n < 1e6 K^1/2 cm^3; this is a safe temperature range because it corresponds to very low densities. The values provided in the table are stated as function of R_v, which is a ratio of visual extinction to reddening. The code currently uses R_v = 3.1, b=6e-5 and B0 (BB for T=3e4 K, cut off at 13.6 eV) from table 2 -- R_v = 3.1 is only appropriate for diffuse cloud. *Further research needed, we might need to replace this.*

On l.1060, a pre-factor of 1e-26 is hardcoded; check if this varies with Z. This pre-factor is also present on l. 1063. This is used in conjunction with he_pe_form (see list of exposed user parameters) and set from Hill+2012.

On l.703, in `get_dust_temperature`, phen_heat and dust_heat represent the same parameter; tdust and tgas cannot drip below he_absTmin, which is set in flash.par. The underlying assumption is that the dust is optically thin to itself; acceptable assumption in our regime from discussion with Simon Glover.
l. 801: The dust cooling rate is set to lambda_dust = 6.8 * tdust**6
There is a note there stating that his will overestimate cooling in wind bubbles and should be switched off above 1e5 K.
The same function also includes collisional cooling of the dust by the gas, for temperature
dust_t = dust_heat (flux - PE heating) + collisional cooling - dust radiative cooling

On l.1600, the piecewise power-law for radiative cooling is defined in the `Radloss` subroutine. It includes bremßtrahlung (with T^1/2) and other processes -- *look into this*. This is likely a piecewise power-law fit to Delgarno and McCray. 


### Wind routine ###

# inject_direct.F90 #
Look into `Particles/ParticlesMain/active/Sink/Couple_AMUSE/wind/inject_direct.F90`.

l. 126 defines a variable `gamma_` -- this is from the runtime parameters. 
There is a constant in the mass-loading routine on l.209, which is the post-shock temperature. This implicitly uses mu = 14/23 = 0.61, which is for fully ionized hydrogen and helium. This may not be appropriate, as several stars are not hot enough to ionize helium. 

The Weaver solution (Weaver+1977, eq. 12) implicitly sets gamma = 5/3, which is based on the adiabatic wind theory by Holzer and Axford 1970. This may have to be changed when we consider molecular hydrogen. See routine on l. 282. Note that this is only called if `variable_radius = .true.`, and the default is `false`. 

On l. 1029, in the sound speed calculation, the constant in the denominator appears to be a hard-coded value of mu. *Check this*. This is within a block with `if use_wind_compute_dt`, and this parameter is set to false by default.

*Check Eos_wrapped routine, as this may have references to metallicity*

### Equation of state ###
This is located in `/physics/Eos/EosMain/Eos_wrapped.F90`, which is vanilla FLASH. We use this wtih `MODE_DENS_EI`, which uses density and internal energy as inputs. *What does FLASH need to know about the metallicity for the equation of state?*

### Sink formation ###
This is located in `/src/flash/source/Particles/ParticlesMain/active/Sink/Couple_AMUSE/Couple_AMUSE_Sinks_and_Stars/Particles_sinkCreateAccrete.F90`.
Set gamma from example on inject_direct.F90 l. 126 AND make sure that the factor of gamma is applied everywhere in the routine. There is no gamma variable in this file. 
The value of gamma is hardcoded on l. 1045.

# Exposed user parameters that will need to be varied #
* Exposed CR parameters
* Gzero
* Scale height h_uv
* he_pe_form for photoelectric heating
* `dust_sputter_temp` may be relevant, but not changed for now

# Parameters that will need to be exposed #
* Dust-to-gas ratio
* Mean molecular weight -- current set in separate place?
* Gamma (also make sure it's used consistently)

# To review in flash.par #
* On l. 703, tolerance and smallt values are hardcoded -- are those also set in flash.par?
* Equation of state -- figure out the comment about double-counting mu

# Parameters to change in flash.par #
Several parameters will need to change in flash.par. It would be easiest to include a script changing them in a consistent manner when people wish to run at lower/higher metallicities.


# Notes from meeting with SCOG #

dust cooling is important for density > 1e5, which we don't really

low T fine structure lines have low critical densities so need density dependent cooling.

CO cooling is primary molecular cooling.

Below 1e4 K Hill et al can't be scaled with metallicity (OK above) because fine structure lines dominate and atomic to molecular transition happens

Check that Neufeld accounts for density

Glover & Clark 12 Fig 4 shows transition from atomic to molecular & fine structure cooling as a function of metallicity.

Treatment of this is available.

Some kind of chemical treatment needed?  Forming molecules in a free-fall time is good at solar, but H2 formation is 100x longer (scales with dust) so equilibrium may never be reached.

fine structure cooling may be sufficient, since we do have ionization.

However adding H2 tracer would require tracking Ly-W radiation as well as formation.

Neufeld+ 96 does act at solar Z (20K with fine structure to 10K with CO).  

Fig 8 of Glover & Clark 12 (the first one) shows transition across different cooling mechanisms.

So fine structure (C II) line is pretty good. Need to update atomic cooling to include Glover function dependent on temperature & density: solve equilibrium level populations.  Include stimulated emission and absorption so CMB at high z can be included.  (Would need to add CMB term).

For MW, effective T of ISR at CMB wavelength is still dominated by CMB. It's 158 micron so hard to get very bright emission from other mechanisms.

[RK: Primordial star clusters with puffy stars to watch collisional growth (which would require updated stellar metallicity-independent (because accretion dominated) mass-radius relation: Nandal et al, Osukawa?) But then need H2 chemistry without a doubt.  So a separate project.]

Dust cooling *can* be scaled by metallicity for our densities and temperatures.  Scaling dust to gas linearly or steeper dependence (transition is at SMC metallicity so should include).  Dust to gas should have a functional dependence therefore (don't worry about variability in dwarfs for now)

**Treatment of atomic cooling from SCOG -- a week or two given teaching, so mid-July**

Molecular cooling is marginally important at solar, less at lower metallicity so it can be treated approximately.  Transition is set by extinction (Av = 2) which means dependence isn't clear.  Density dependent threshold is crude approximation.

**mid-August or second half of September (after 11 Sept) for a C-C visit to HD** (Ralf will be in AUS end of Sept or early October; AG meeting in Garching is another target.)

metallicity propagating into binary evolution is interest of CC-C.

# Meeting notes (Claude) #
Likely metal dominated at our metallicities; don't need to forbid about dust coupling for our density range. 

Need to take both temperature and density into account for cooling. Molecular cooling should already account for that density-dependence. 

Atomic cooling is currently not metallicity-dependent; this will not work at low metallicity. We will need low temperature atomic cooling. See Glover & Clark 2012b, fig. 4 (MNRAS 421, 116). At low metallicity, the atomic to molecular transition is at higher densities; the low density gas needs a temperature and density dependence. We currently have no chemistry, *need a chemistry treatment* to account for this properly. Assuming chemistry tracks density is ok at solar abundance but at low Z, you cannot assume this. The chemistry timescales are longer at low Z, longer than e.g. cloud free-fall time. You therefore need some knowledge of history. 

We may be able to start from the tracer fields (e.g. by using one to trace HII). We would also need to account for the Lyman-Werner band. 

Can just use the ionization fraction, and assume that all the chemistry happens at higher densities than what we reach. Everything would then depend on local temperature and density. 

What about molecular cooling? Do we turn it off/taper it off? From Simon, you can likely turn it off; at solar Z, CO cooling brings you down to T~10 K rather than T~20 K. This becomes less important as you move to below solar.

See Fig. 8 of the same paper for contributions as a function of density. See what goes on around 1e4 cm-3 in our simulations; this is where CO would be important.

Way forward is to include the fine structure cooling (which would be an update of the atomic cooling); this would include stimulated emission and absoption terms. This would give us the high-z CMB contribution almost for free.

Ralf suggests looking at puffy stars and repeated collisions in primordial star clusters.

What to do:
* Update dust cooling: just make it metallicity-dependent. At our densities, the temperature is set by the radiation field (rather than collisions with the gas), so we only need to adjust how much we have. How should we scale the dust-to-gas ratio with metallcity? Observations suggest steeper than linear below SMC metallicity; also more scatter in low-Z galaxies. We could have a dust-to-gas ratio which is separate from the metallicity so that they could be varied separately.
* From Simon Glover, we will get a treatment of the atomic cooling, which will be packaged up and documented a bit (1-2 weeks).
* Molecular cooling can be removed/ignored at sub-solar metallicities. Could we scale Neufeld or taper it off? Simon: the densities don't shift directly in a metallicity-dependent way, but it depends on extinction (e.g. A_v = 2) -- for same gas morphology, the behaviour does scale. _Follow up on this once we understand better_, but taking a density-dependent threshold might be a good, although crude, start -- better than not doing anything.



## meeting on July 10
Shyam has developed simple molecular network tied to VETTAM that might serve to replace the options developed by Simon for fine structure atomic cooling.  This should go public within a week or two after which we can evaluate.


### VETTAM ###


Start in the top level directory, `RadTrans/RadTransMain/VETTAM` with the `Config` file.

# Semenov opacities #
Semenov opacities are available, but are not the default. The underlying assumption in using fixed opacitiies is that the opacity depends only on the gas-to-dust ratio. Note that this ratio scales superlinearly with metallicity at low metallicity; there is a kink in the distribution around 0.1 Zsun. See https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..31R/abstract (Fig. 4)

# Lyman-Werner band #
Adding the Lyman-Werner band would be necessary if we use Shyam Menon's simple chemistry model. 

# rt_ionHeatCool.F90 #
This doing radiation heating of the gas.
Hardcoded parameters: Tneutral = 10.0 K and Tion = 1.e4 K. Those are metallicity dependent and should be set as runtime parameters. 
This uses either ENER_VAR or TGAS_VAR (l. 53) --> Tion and Tneutral are used only when TGAS_VAR is used. Add if def for the defitions above.

# rt_ionise.F90 #
*Those notes were from the /OnlyH directory, which is not used.*
_Look for rt_ionisehydro.F90, which supplies the nenergy and momentum terms on gas due to photoionization_. 
We are using `IHA_SPEC` and `IHP_SPEC`; we might implicitly be using some pieces of KROME. The comments at the top of the file claim that IHA_SPEC points to the KROME network. On l. 117 to l. 120, mu is set from the mu of ionized and neutral hydrogen. Those values, hA and hpA, are set *somewhere unknown -- figure it out*. 

The `Config` file in the `/VETTAM/Photoionization` directory contains the default runtime parameters associated with VETTAM. 

*Now looking at the thermochemistry directory*
This uses a recombination coefficient which depends on temperature.
The subroutine SetIonRates calculates ionization rates based on UEUV_VAR (which is in Flash.h) and hnu.

The energy per ionization of H is set in `Photoionization/Config`. It is set from Kim+2023 as an average value for a star cluster. hnu is also set here. *This must be investigated in greater detail, and how it connects to the eion, epep values we set in Torch. This is critical, as if those values overwrite the eion, epep values, we are losing the spectral information we put in from the stellar evolution.* If we understand correctly what is happening there, the energy per photon does not vary with the stellar source, although the number of photons in each band does. Adding extra bands could alleviate those concerns, as the cross-section varies by band.

# rt_ioniseData.F90 #
Nothing relevant in this file.

# Photoionization/rt_ionisemodule.F90 #
This calculates the recombination coefficients. Those depend on temperature but not strongly on metallicity; the only effect would be the extra electrons from the ionized metals.

# Photoionization/rt_ionMomentum.F90 #
This calculates energies and opacities. *Check opacities*.

# VETTAM/IonHeatCool.F90#
On l. 236, `mu_mol = 1.3` is hardcoded.

# RadTrans_computeDt.F90 #
On l. 66, `mu = 0.61` is hardcoded; it is used on l. 215, 247.
On l. 229, a CFL of 0.3 is hardcoded in the dt_wind routine. Also see l. 271, there is a pre-factor of 0.3 hardcoded in the dt_min_local calculation. Note that on l. 267, cfl_radPressure is used, which is a runtime parameter. If we try to reduce the CFL is hot zone or raise it to speed up the code, this change will not propagate to the routine. On l. 267, the calculation for the dt from momentum injection includes dx**4 -- where is the power of 4 from?

The parameter `rt_gamma1` is set from `gamma - 1`, where `gamma` is a runtime parameter.

# RadTRans_data.F90 #
Nothing relevant in this file.

# RadTrans_finalize.F90 #
Nothing relevant in this file.

# RadTrans_init.F90 #
This gets the runtime parameters.

# RadTrans_interface.F90 #
Nothing relevant in this file.

# RadTrans.F90 #
Nothing relevant in this file.

# rt_data.F90 #
This sets several constants and conversion factors. 
The hydrogen recombination coefficients are set here (on l.71-72) but are not used anywhere as the values are overwritten by temperature-dependent values. Clean up.
On l. 96, `sigDust = 1e-21` is hardcoded; there might be a dependence on the dust-to-gas ratio here (see Draine 2011). On l. 97, the dust-to-gas ratio is hardcoded to 0.01. Note that both those parameters are already available in flash.par and called in  `rt_init.F90`.

# rt_dustTemperature.F90 #
 On l.494, `eos_gamma` from `Eos_data` and l. 496 `eos_gammam1` from `eos_idealGammaData` -- what are those, and how/hy are they different from the gamma - 1 above?

 On l. 585, we call `Eos_getAbarZbar` -- is this an average atomic weight? On l. 584, we use `GAMC_VAR` -- is this a field variable for gamma? Can we use it?

 _Note_: Look into the opacity calculation routine. How do we get `TAUP` and `TAUR`?

 # rt_dustTerms.F90 #
 Nothing relevant in this file.

 # rt_fillMatrix.F90 #
 Nothing relevant in this file.

 # rt_init.F90 #
 Constants and runtime parameters. On l.67, `rt_abar = 1.0 + rt_abundM*rt_metal` is set; does `Eos_getAbarZbar` reference this value?
 The carbon abundance is hardcoded to  `abu_c = 7.1e-7` on l.77 but *is never used in this directory*. It is used in calc_ionization.F90.

 # rt_petsc.F90 #
 Nothing relevant in this file.

 # rt_sedEddTensor.F90 #
Nothing relevant in this file.

# rt_setOpacity.F90 #
This calls `sim_A_n` and `sim_A_i`, which are set in flash.par. Those are atomic weights; the default values are 14/11 (10% He by number, neutral) and 14/21 (ionized hydrogen, neutral helium, assuming the same abundances). _Note that this is inconsistent with the values set for `rt_abundM` and `rt_metal` in the default flash.par_. Note also that the default value (25% by weight) holds across cosmic time; there may be a mild scaling with metallicity.

`dusttoGasRatio` is called here. It is called as a runtime parameter on l.82 in `RadTrans_init.F90`. Note that this is *a different dust-to-gas ratio from the one used above.* This ratio is set to 1 (i.e solar, since it is normalized to solar) by default, while the default value for the other one is 0.01 (also solar). 

# rt_sinkHydro.F90 #
Nothing relevant in this file.

# rt_sinkInject.F90 #
This calls luminosities per band. We are in the case where NION, NPEP, EION and EPEP are set, which means that we use the values set from SeBa.

# calc_ionization.F90 #
This does not appear to depend on metallicity. We may want to eventually double-check the implicit equation.

### Equation of state (in vanilla FLASH) ###

# /physics/Eos/EosMain/Gamma/Eos_idealGammaData.F90 #
`gamma` is a runtime parameter; anything coming from EOS will use that gamma.

# /physics/Eos/EosMain/Gamma/Eos_getAbarZbar.F90 #
Looking into Simulation_initSpecies.F90, we see that neutral and ionized species are set them. The machinery exists to have separate values of gamma for neutrals and ions but we have not used it so far. Note that we are not setting `sim_gamma_n` or `sim_gamma_i`, as those are the values used if we use a variable gamma.

Things to look at:
* Multispecies.h
* calc_ionization
* EOS






