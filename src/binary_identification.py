import numpy as np
from amuse.datamodel import Particles
from amuse.ext.orbital_elements import get_orbital_elements_from_binaries
from amuse.units import units
from sklearn.neighbors import KDTree

def get_binaries_from_stars(stars, num_stars_for_tree = 10, a_max = 1e4):

    _pmass    = np.copy(stars.mass.value_in(units.MSun))
    _cmass    = np.zeros(len(stars.tag))
    _p_id     = np.array(range(len(stars.tag)))
    _c_id     = np.zeros(len(stars.tag))
    _energies = np.zeros(len(stars.tag))
    _a        = np.zeros(len(stars.tag))
    _e        = np.zeros(len(stars.tag))

    def _relative_v2(star, stars):
        v2 = (star.vx-stars.vx)**2 + (star.vy-stars.vy)**2 + (star.vz-stars.vz)**2
        return v2

    def _relative_r(star, stars):
        r2 = (star.x-stars.x)**2 + (star.y-stars.y)**2 + (star.z-stars.z)**2
        r = np.sqrt(r2.value_in((units.m)**2)) | units.m
        return r

    def _binding_energy(star, stars):
        Mmult = (star.mass*stars.mass)
        Madd  = (star.mass+stars.mass)
        v2 = _relative_v2(star, stars)
        r  = _relative_r(star, stars)
        E = Mmult * v2 / (2 * Madd) - ((units.constants.G * Mmult) / r)
        return E

    def _get_companion(star, stars):
        E = _binding_energy(star, stars).flatten()
        c = np.argmin(E.value_in(units.erg))
        return c, E[c].value_in(units.erg), stars[c].mass.value_in(units.MSun)

    if len(stars) <= num_stars_for_tree:
        # Do N^2 search if fewer stars than number used for tree

        for s in _p_id:
            s_ids = np.where(_p_id != s)[0]
            _c_id_tmp, _energies[s], _cmass[s] = _get_companion(stars[s], stars[s_ids])
            if _energies[s] == -1*np.inf:
                print('Stars at same location!', _cmass[s], _pmass[s])
            _c_id[s] = s_ids[_c_id_tmp]
        _c_id = _c_id.astype('int')
        _a, _e = get_orbital_elements_from_binaries(stars[_p_id], stars[_c_id])[2:4]
        _a = _a.value_in(units.au)

    else:
        # Query the tree
        stars_for_tree = np.transpose(np.vstack((stars.x.value_in(units.pc), stars.y.value_in(units.pc), stars.z.value_in(units.pc))))
        tree = KDTree(stars_for_tree)

        nearest_dist, _nearest_ind = tree.query(stars_for_tree, k=num_stars_for_tree + 1)  # k=2 nearest neighbors where k1 = identity
        nearest_ind = _nearest_ind[:, 1:]

        for s in _p_id:
            s_ids = nearest_ind[s]
            _c_id_tmp, _energies[s], _cmass[s] = _get_companion(stars[s], stars[s_ids])
            if _energies[s] == -1*np.inf:
                print('Stars at same location!', _cmass[s], _pmass[s])
            _c_id[s] = s_ids[_c_id_tmp]
        _c_id = _c_id.astype('int')
        _a, _e = get_orbital_elements_from_binaries(stars[_p_id], stars[_c_id])[2:4]
        _a = _a.value_in(units.au)

    # Keep only bound pairs
    select_by_boundedness = np.intersect1d(np.where(_energies < 0)[0], np.where(_a <= a_max)[0])
    _pmass = _pmass[select_by_boundedness]
    _cmass = _cmass[select_by_boundedness]
    _p_id = _p_id[select_by_boundedness]
    _c_id = _c_id[select_by_boundedness]
    _energies = _energies[select_by_boundedness]
    _a = _a[select_by_boundedness]
    _e = _e[select_by_boundedness]
    # At this point, this seems OK but duplicated

    # Here select to ensure the primary is the most massive
    select_by_primary = np.where(_pmass > _cmass)[0]
    _pmass = _pmass[select_by_primary]
    _cmass = _cmass[select_by_primary]
    _p_id = _p_id[select_by_primary]
    _c_id = _c_id[select_by_primary]
    _energies = _energies[select_by_primary]
    _a = _a[select_by_primary]
    _e = _e[select_by_primary]

    # Sort by energies
    sort_by_energy = np.argsort(_energies)
    _pmass = _pmass[sort_by_energy]
    _cmass = _cmass[sort_by_energy]
    _p_id = _p_id[sort_by_energy]
    _c_id = _c_id[sort_by_energy]
    _energies = _energies[sort_by_energy]
    _a = _a[sort_by_energy]
    _e = _e[sort_by_energy]

    # Remove duplicate companions & primaries
    select_by_companion = []
    seen_c = set()
    seen_p = set()
    for i, (c_id, p_id) in enumerate(zip(_c_id, _p_id)):
        if p_id not in seen_c and c_id not in seen_p and c_id not in seen_c:
            # Accept system if primary is not previous companion and companion is not
            # in previously seen pair (including rejected pairs).
            select_by_companion.append(i)
        seen_c.add(c_id)
        seen_p.add(p_id)
    select_by_companion = np.array(select_by_companion, dtype=int)

    _pmass = _pmass[select_by_companion]
    _cmass = _cmass[select_by_companion]
    _p_id = _p_id[select_by_companion]
    _c_id = _c_id[select_by_companion]
    a = _a[select_by_companion] | units.au
    e = _e[select_by_companion]

    # Prepare particle set for saving
    binaries = Particles(semi_major_axis=a, eccentricity=e)
    binaries.child1 = list(stars[_p_id])
    binaries.child2 = list(stars[_c_id])
    
    return binaries
