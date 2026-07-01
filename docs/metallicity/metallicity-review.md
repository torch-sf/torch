### Metallicity Review ###
M-MML starting 23 Jun 2026, then edited by CC-C and M-MML on 29 Jun 2026. _Note that this was created based off of main, and the notes here will be moved to branch created from develop._

## RadHeat ##
active version I believe to be `src/flash/source/physics/sourceTerms/Heat/HeatMain/HeatCool/phenHeat/mol_and_dust/solver/RadHeat.F90`

The header file is `src/flash/source/physics/sourceTerms/Heat/HeatMain/HeatCool/phenHeat/Heat_Interface.f90`

### Cooling ###
Cooling uses dust cooling from Goldsmith (2001; https://scixplorer.org/abs/2001ApJ...557..736G/abstract).  This is implicitly solar.  To replace that, we can look at Dopcke et al. (2011;  https://scixplorer.org/abs/2011ApJ...729L...3D/abstract).  I've touched base with Ralf K, who suggested emailing Simon G., which I've done (24 June)


# Heating #
Note: _There is a commented out function to fill the guard cells._
This uses constants.h, which is located at `/src/flash/source/Simulation/SimulationMain/StratBox/constants.h`. _Why is this in the StratBox directory?_ In vanilla FLASH, this is in the top-level `/Simulation` directory. There are no differences between the files.

heat_data is a data file; this will likely be part of replacing the data files from Robi.

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

# Exposed user parameters that will need to be varied #
* Exposed CR parameters
* Gzero
* Scale height h_uv
* he_pe_form for photoelectric heating
* `dust_sputter_temp` may be relevant, but not changed for now

# Parameters that will need to be exposed #
* Dust-to-gas ratio
* Mean molecular weight -- current set in separate place?

# To review in flash.par #
On l. 703, tolerance and smallt values are hardcoded -- are those also set in flash.par?


# Notes from meeting with SCOG #

dust cooling is important for density > 1e5, which we don't really

low T fine structure lines have low critical densities so need density dependent cooling.

CO cooling is primary molecular cooling.

Below 1e4 K Hill et al can't be scaled with metallicity (OK above) because fine structure lines dominate and atomic to molecular transition happens

Check that Neufeld accounts for density

Glover & Clark 12 Fig 4 shows transition from atomic to molecular & (fine structure cooling as a function of metallicity.

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

Atomic cooling is currently not metallicity-dependent; this will not work at low metallicity. We will need low temperature atomic cooling. See Glover & Clark 2012b, fig. 4. At low metallicity, the atomic to molecular transition is at higher densities; the low density gas needs a temperature and density dependence. We currently have no chemistry, *need a chemistry treatment* to account for this properly. Assuming chemistry tracks density is ok at solar abundance but at low Z, you cannot assume this. The chemistry timescales are longer at low Z, longer than e.g. cloud free-fall time. You therefore need some knowledge of history. 

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







