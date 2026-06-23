from scipy.interpolate import RegularGridInterpolator

import numpy as np
import os
import re

class YieldSource:
    """
    Class for interpolating yields from a source.
    Creates a RegularGridInterpolator object for given elements, parameters
    and yields and provides functionallity to interpolate yields for either
    specific elements or total mass loss (sum of all elements).

    This class is intended to be used when creating structure for specific
    yields table.
    """

    def __init__(self, elements, params, yields):
        """
        Initialize RegularGridInterpolator object for element yields and total
        mass loss.

        Positional Arguments:
            elements - List of elements in provided yield table
            params - Parameter space to interpolate in
            yields - Array of yields, matching list of elements and parameters
        """

        self.params = params
        self.yields = {element: RegularGridInterpolator(self.params, yields[i],
                                                        fill_value=None, bounds_error=False)
                       for i, element in enumerate(elements)}
        self.mloss = RegularGridInterpolator(self.params, np.sum(yields[np.array(elements, dtype=str)!='Z'], axis=0),
                                             fill_value=None, bounds_error=False)

    def get_yld(self, elements, params, interpolate='nearest', extrapolate=False):
        """
        Return total yields for list of elements provided at parameters.

        Positional Arguments:
            elements - List of elements
            params - Points in parameters space to interpolate.
        Keyword Arguments:
            interpolate - Interpolation method, see scipy.interpolate.RegularGridInterpolator
            extrapolate - If true, allow extrapolation outside table parameter space. Use caution.
        """

        elements = np.atleast_1d(elements)

        if len(params) != len(self.params):
            raise ValueError(
                "Supplied parameters do not match yield set parameters.")

        points = self.convert2array(params)

        if not extrapolate:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip])
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip])
        else:
            warnings.warn("Extrapolating yields might lead to problematic behaviour (e.g., negative yields). Ensure that yields behave as expected.")

        if len(elements) == 1:
            try:
                return self.yields[elements[0]](points, method=interpolate)
            except KeyError:
                return np.nan
        else:
            if all(element in list(self.yields.keys()) for element in elements):
                return np.array([self.yields[element](points, method=interpolate) for element in elements])
            else:
                yld = []
                shape = self.yields['H'](points, method=interpolate).shape
                for element in elements:
                    try:
                        yld.append(self.yields[element](points, method=interpolate))
                    except KeyError:
                        yld.append(np.ones(shape) * np.nan)
                return yld

    def get_mloss(self, params, interpolate='nearest', extrapolate=False):
        """
        Return sum of all yields (mass loss) for at parameters.

        Positional Arguments:
            params - Points in parameters space to interpolate.
        Keyword Arguments:
            interpolate - Interpolation method, see scipy.interpolate.RegularGridInterpolator
            extrapolate - If true, allow extrapolation outside table parameter space. Use caution.
        """

        if len(params) != len(self.params):
            raise ValueError(
                "Supplied parameters do not match yield set parameters.")

        points = self.convert2array(params)

        if not extrapolate:
            for ip in range(len(params)):
                points[:, ip][(points[:, ip] < np.min(self.params[ip]))] = np.min(
                    self.params[ip])
                points[:, ip][(points[:, ip] > np.max(self.params[ip]))] = np.max(
                    self.params[ip])

        return self.mloss(points, method=interpolate)

    @staticmethod
    def convert2array(params):
        """
        Helper function to restructure parameters to data points for interpolation.
        """
        max_length = max([len(param) if isinstance(
            param, (list, np.ndarray)) else 1 for param in params])

        # Ensure all arguments are lists or arrays of the same length
        args = []
        for param in params:
            if isinstance(param, (list, np.ndarray)):
                if len(param) != max_length:
                    raise ValueError(
                        "All list of parameters must have the same length.")
                args.append(np.array(param))
            else:
                # Fill with the same value
                args.append(np.full(max_length, param))

        # Combine arguments into points for interpolation
        return np.stack(args, axis=-1)


class LimongiChieffi2018:
    """
    Structures and provides functions to work with yields from Limongi & Chieffi (2018).
    This assumes fall-back and mixing (Umeda & Nomoto, 2002) and direct collapse for stars
    with mass > 25Msun.

    If using this class, please cite Limongi & Chieffi, 2018, ApJS, 237, 13L

    To use this class, you must download the yield tables available at
        https://orfeo.iaps.inaf.it/

    Class properties:
        models            - List of models available
        mass              - Array of tabulated stellar masses
        metal             - Array of tabulated stellar metallicities
        rot               - Array of tabulated stellar rotation velocities
        ccsn_mmax         - Maximum mass for core collapse (otherwise direct collapse)
        filedir           - Directory of tabulated data
        yield_tablefile   - Filename for total yields
        wind_tablefile    - Filename for wind yields
        elements          - List of available elements
        atomic_num        - List of atomic numbers corresponding to element list
        wind              - Object holding interpolation points and data for winds
        ccsn              - Object holding interpolation points and data for core-collapse SNe

    Usage:
        Loads data from yield tables and provides function to simplify
        interpolation for list of elements at given points in mass, metallicity
        and rotation via functions:

            ccsn_yields(elements, mass, metal, rot)
            wind_yields(elements, mass, metal, rot)
            total_yields(elements, mass, metal, rot)

        Similarly, total mass loss rates are provided with functions:
            total_mloss(mass, metal, rot)
            ccsn_mloss(mass, metal, rot)
            wind_mloss(mass, metal, rot)
    """

    def __init__(self, model='R', path='$TORCH_DIR/data/yield_tables/'):
        """ Initialize class, loading tables and setting up objects 
            for interpolation.
        """

        self.models = ['F', 'I', 'M', 'R']
        self.mass = np.array(
            [13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0]) # [Msol]
        self.metal = np.array([1.345e-2, 3.236e-3, 3.2363e-4, 3.236e-5])
        self.rot = np.array([0, 150, 300]) # [km/s]
        self.ccsn_mmax = 25  # [Msol]

        if model not in self.models:
            raise ValueError("Model does not exist.")
        
        self.filedir = os.path.expandvars(path)

        self.yield_tablefile = self.filedir + \
            f'/LC2018/tab_{model}/tab_yieldstot_ele_exp.dec'
        self.wind_tablefile = self.filedir+f'LC2018/tab_{model}/tabwind.dec'

        self.elements, self.atomic_num = self.get_element_list()
        
        # Load tables
        wind_yld = self.load_wind_yields()
        ccsn_yld = self.load_ccsn_yields()
        
        # Add metallicity variable which is sum of all elements except first two (H, He)
        self.elements.insert(0, 'Z')
        self.atomic_num = np.insert(self.atomic_num, 0, 0)
        wind_yld = np.concatenate((wind_yld[2:].sum(axis=0, keepdims=True), wind_yld),axis=0)
        ccsn_yld = np.concatenate((ccsn_yld[2:].sum(axis=0, keepdims=True), ccsn_yld),axis=0)

        # Wind yield interpolation object
        self.wind = YieldSource(
                self.elements, [self.rot, self.metal, self.mass], wind_yld)

        # Core-collapse SNe yield interpolation object
        self.ccsn = YieldSource(
                self.elements, [self.rot, self.metal, self.mass], ccsn_yld)

    def ccsn_yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returns yields [Msol] for elements from core-collapse supernovae given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.ccsn.get_yld(elements, [rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def wind_yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returns yields [Msol] for elements from pre-supernovae wind given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.wind.get_yld(elements, [rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def total_yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning yields [Msol] for each element in elements, given
            stellar parameters mass [Msol], metal [mass fraction], and rot [km/s]."""
        args = (elements, mass, metal, rot, interpolate, extrapolate)
        return self.wind_yields(*args) + self.ccsn_yields(*args)

    def ccsn_mloss(self, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements from core-collapse supernovae given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.ccsn.get_mloss([rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def wind_mloss(self, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements from pre-supernovae wind given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.wind.get_mloss([rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def total_mloss(self, mass, metal, rot, source='all', interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements, given
            stellar parameters mass [Msol], metal [mass fraction], and rot [km/s]."""

        return self.wind_mloss(mass, metal, rot,
                               source=source,
                               interpolate=interpolate,
                               extrapolate=extrapolate) \
            + self.ccsn_mloss(mass, metal, rot,
                              source=source,
                              interpolate=interpolate,
                              extrapolate=extrapolate)

    def get_element_list(self):
        """ Helper function for reading element list during initialization.
        """

        with open(self.wind_tablefile, 'r') as file:
            lines = file.readlines()

        elements = []
        atomic_num = []
        for line in lines[1:]:
            if line.split()[0] == 'ele':
                return elements, np.array(atomic_num, dtype=float)
            element = ''.join(re.findall(r'[a-zA-Z]', line.split()[0]))
            if element not in elements:
                elements.append(element)
                atomic_num.append(int(line.split()[1]))

    def load_wind_yields(self):
        """ Loader function for wind yields. Used during initialization.
        """

        wind_yld = np.zeros([len(self.elements), self.rot.size,
                             self.metal.size, self.mass.size])

        with open(self.wind_tablefile, 'r') as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self.get_metal_index_from_model(model[3])
                ind_rot = self.get_rot_index_from_model(model[4:])

                data = np.genfromtxt(self.wind_tablefile,
                                     usecols=[1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                     skip_header=index+1,
                                     max_rows=142).T
                for i, (atom_nr, element) in enumerate(zip(self.atomic_num, self.elements)):
                    mask = (data[0] == atom_nr)
                    wind_yld[i, ind_rot, ind_metal] = np.sum(
                        data[1:, mask], axis=1)

        return wind_yld

    def load_ccsn_yields(self):
        """ Loader function for SNe yields. Used during initialization.
        """

        wind_yld = self.load_wind_yields()
        total_yld = np.zeros([len(self.elements), self.rot.size,
                              self.metal.size, self.mass.size])

        with open(self.yield_tablefile, 'r') as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self.get_metal_index_from_model(model[3])
                ind_rot = self.get_rot_index_from_model(model[4:])

                total_yld[:, ind_rot, ind_metal, :] = np.genfromtxt(self.yield_tablefile,
                                                                    usecols=[
                                                                        4, 5, 6, 7, 8, 9, 10, 11, 12],
                                                                    skip_header=index+1,
                                                                    max_rows=53)

        ccsn_yld = total_yld - wind_yld
        ccsn_yld[ccsn_yld < 0.0] = 0.0
        ccsn_yld[:, :, :, (self.mass > self.ccsn_mmax)] = 0.0

        return ccsn_yld

    @staticmethod
    def get_metal_index_from_model(model):
        """ Internal helper function.
        """
        if model == 'a':
            return 0
        elif model == 'b':
            return 1
        elif model == 'c':
            return 2
        elif model == 'd':
            return 3
        else:
            raise ValueError("Model does not exist.")

    @staticmethod
    def get_rot_index_from_model(model):
        """ Internal helper function.
        """
        if model == '000':
            return 0
        elif model == '150':
            return 1
        elif model == '300':
            return 2
        else:
            raise ValueError("Model does not exist.")
