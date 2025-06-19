# for the arepo port
from amuse.community import *
from amuse.community.interface.hydro import HydrodynamicsInterface
from amuse.support.options import OptionalAttributes, option
from amuse.units import generic_unit_system
from amuse.units import si
from amuse.community.interface.common import CommonCode
import numpy as np
import os

length = generic_unit_system.length
time = generic_unit_system.time
mass = generic_unit_system.mass
speed = generic_unit_system.speed
density = generic_unit_system.density
momentum =  generic_unit_system.momentum_density
energy =  generic_unit_system.energy_density
enerInt = generic_unit_system.length ** 2 / generic_unit_system.time ** 2
potential_energy =  generic_unit_system.energy
magnetic_field = generic_unit_system.mass / generic_unit_system.current / generic_unit_system.time ** 2
acc = generic_unit_system.acceleration
potential = generic_unit_system.potential

class FlashInterface(CodeInterface, HydrodynamicsInterface):

    use_modules = ['flash_run']

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="voramr_worker", **keyword_arguments)


    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def cleanup_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_grid_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx', 'rhovy', 'rhovz', 'rhoen']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rho','rhovx', 'rhovy', 'rhovz', 'rhoen']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function
# New setter, under active developement - SCL 11/08/22    
    @legacy_function
    def set_block_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('dataSize', dtype='i', direction=function.IN)
        for x in ['rho', 'vx', 'vy', 'vz', 'eint', 'gpot']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid', 'nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx', 'rhovy', 'rhovz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_momentum_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['rhovx', 'rhovy', 'rhovz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rhoen', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_energy_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rhoen', dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rho', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def set_grid_density():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('rho', dtype='d', direction=function.IN)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_number_of_grids():
        function = LegacyFunctionSpecification()
        function.addParameter('nproc', dtype='int32', direction=function.IN)
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_range():
        function = LegacyFunctionSpecification()
        for x in ['nx','ny','nz']:
            function.addParameter(x, dtype='int32', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='int32', direction=function.IN)
        function.addParameter('nproc', dtype='int32', direction=function.IN)
        function.result_type='int32'
        return function

    @legacy_function
    def get_position_of_index():
        """
        Retrieves the x, y and z position of the center of
        the cell with coordinates i, j, k in the grid specified
        by the index_of_grid
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('index_of_grid', dtype='i', direction=function.IN, default = 1)
        function.addParameter('nproc', dtype='i', direction=function.IN, default = 0)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_index_of_position():
        """
        Retrieves the i,j and k index of the grid cell containing the
        given x, y and z position, the index of the grid and the local
        processor number on which this grid resides.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['i','j','k']:
            function.addParameter(x, dtype='i', direction=function.OUT)
        function.addParameter('index_of_grid', dtype='i', direction=function.OUT)
        function.addParameter('nproc', dtype='i', direction=function.OUT)
        function.addParameter('n', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_leaf_indices():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('dummy', dtype='i', direction=function.IN)
        function.addParameter('ind', dtype='i', direction=function.OUT)
        function.addParameter('ret_cnt', dtype='i', direction=function.OUT)
        function.addParameter('num_of_blks', dtype='i', direction=function.OUT)
        function.addParameter('nparts',dtype='i', direction=function.LENGTH)
        function.result_type='i'
        return function

    @legacy_function
    def get_max_refinement():
        function = LegacyFunctionSpecification()
        function.addParameter('max_refine', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def grid_update_refinement():
        function = LegacyFunctionSpecification()
        function.addParameter('gridChanged', dtype='b', direction=function.OUT)
        function.result_type='i'
        return function

    # This needs to look like "get_grid_density" etc etc.
    @legacy_function
    def get_cell_volume():
        function = LegacyFunctionSpecification()
        for x in ['block','i','j','k']:
            function.addParameter(x, dtype='int32', direction=function.IN)
        function.addParameter('vol', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_number_of_procs():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_all_local_num_grids():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('num_grids_array', dtype='i', direction=function.INOUT)
        function.addParameter('nprocs', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    # WORK IN PROGRESS!!!!
    #@legacy_function
    #def get_data_all_local_blks():
        #function = LegacyFunctionSpecification()
        #function.must_handle_array = True
        #function.addParameter('data_array', dtype='d', direction=function.INOUT)
        #function.addParameter('numcells', dtype='i', direction=function.LENGTH)
        #function.result_type = 'i'
        #return function

    @legacy_function
    def get_1blk_cell_coords():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('axis', dtype='i', direction=function.IN)
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('limits', dtype='i', direction=function.IN)
        function.addParameter('coords', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

########################
# IO Stuff
########################

    @legacy_function
    def write_chpt():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def IO_out():
        function = LegacyFunctionSpecification()
        function.addParameter('output_type', dtype='s', direction=function.IN)
        function.addParameter('fileNumber', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def IO_num():
        function = LegacyFunctionSpecification()
        function.addParameter('output_type', dtype='s', direction=function.IN)
        function.addParameter('fileNumber', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_output_dir_wrapped():
        function = LegacyFunctionSpecification()
        function.addParameter('output_dir', dtype='s', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_runtime_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('rt_name', dtype='s', direction=function.IN)
        function.addParameter('rt_value', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def timer_summary():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

class Flash(CommonCode):



    def __init__(self, unit_converter = None, **options):
        self.unit_converter = unit_converter
        self.stopping_conditions = StoppingConditions(self)

        CommonCode.__init__(self,  FlashInterface(**options), **options)

#        self.set_parameters_filename(self.default_parameters_filename)

    def define_converter(self, object):
        if self.unit_converter is None:
            return

        object.set_converter(self.unit_converter.as_converter_from_si_to_generic())


    def get_index_range_inclusive(self, index_of_grid = 1, nproc=0):
        nx, ny, nz = self.get_grid_range(index_of_grid, nproc)
        return (1, nx, 1, ny, 1, nz)

    # I think all of these should start with a numpy array like:
    # three_vector = np.zeros((n,3))
    # even if n = 1, so that the returned shape is always the same
    # regardless of n and therefore the loop behavior is always the same
    # when looping over these. -JW

    def get_output_dir(self):
        output_dir = self.get_output_dir_wrapped()
        output_dir = output_dir.strip()

        return output_dir

    def define_methods(self, object):

        #length = units.cm
        #time = units.s
        #mass = units.g
        #speed = units.cm*units.s**-1
        #density = units.g*units.cm**-3
        #momentum =  density*speed
        #energy =  units.cm**2*units.g*units.s**-2
        #potential_energy =  energy
        #magnetic_field = units.g*0.1*units.C**-1*units.s**-1
        #acc = units.cm*units.s**-2

### These two are included in CommonCode

        #object.add_method(
            #'initialize_code',
            #(),
            #(object.ERROR_CODE)
        #)


        #object.add_method(
            #'cleanup_code',
            #(),
            #(object.ERROR_CODE)
        #)
        object.add_method(
            'evolve_model',
            (time,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_position_of_index',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length, length, length, object.ERROR_CODE,)
        )

        object.add_method(
            'get_index_of_position',
            (length, length, length),
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX,
             object.INDEX, object.ERROR_CODE,)
        )

        object.add_method(
            "get_max_refinement",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density, momentum, momentum, momentum, energy,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_state',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            density, momentum, momentum, momentum, energy),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'set_block_state',
            (object.INDEX, object.INDEX,
            density, speed, speed, speed, enerInt, potential),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            ( energy,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_energy_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
             energy),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (density,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            density),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (momentum, momentum, momentum,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_momentum_density',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            momentum, momentum, momentum),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_velocity',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (speed, speed, speed,
            object.ERROR_CODE,)
        )

        object.add_method(
            'set_grid_velocity',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX,
            speed, speed, speed),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_cell_volume',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length**3, object.ERROR_CODE,)
        )

        object.add_method(
            'get_number_of_procs',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_all_local_num_grids',
            (object.INDEX),
            (object.INDEX, object.ERROR_CODE)
        )

        object.add_method(
            'get_data_all_local_blks',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_1blk_cell_coords',
            (object.NO_UNIT,object.NO_UNIT,object.NO_UNIT),
            (length, object.ERROR_CODE)
        )

        object.add_method(
            'get_leaf_indices',
            (object.NO_UNIT),
            (object.NO_UNIT, object.INDEX, object.INDEX, object.ERROR_CODE)
        )

        object.add_method(
            "grid_update_refinement",
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )

    def specify_grid(self, definition, index_of_grid = 1, nproc=0):
        definition.set_grid_range('get_index_range_inclusive')

        definition.add_getter('get_position_of_index', names=('x','y','z'))

        definition.add_getter('get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_block_state', names=('rho', 'vx', 'vy', 'vz', 'energy', 'gpot')) # new setter - SCL 11/08/22

        definition.add_getter('get_grid_density', names=('rho',))
        definition.add_setter('set_grid_density', names=('rho',))

#       if self.mode == self.MODE_SCALAR:
#           definition.add_getter('get_grid_scalar', names=('scalar',))
#           definition.add_setter('set_grid_scalar', names=('scalar',))

        definition.add_getter('get_grid_momentum_density', names=('rhovx','rhovy','rhovz'))
        definition.add_setter('set_grid_momentum_density', names=('rhovx','rhovy','rhovz'))

        #definition.add_getter('get_grid_velocity', names=('vx','vy','vz'))
        #definition.add_setter('set_grid_velocity', names=('vx','vy','vz'))

        definition.add_getter('get_grid_energy_density', names=('energy',))
        definition.add_setter('set_grid_energy_density', names=('energy',))


#       definition.add_getter('get_grid_gravitational_potential', names=('gravitational_potential',))
#       definition.add_getter('get_grid_gravitational_acceleration', names=('gravitational_acceleration_x','gravitational_acceleration_y','gravitational_acceleration_z',))

#        definition.add_getter('get_grid_acceleration', names=('ax','ay','az'))
#        definition.add_setter('set_grid_acceleration', names=('ax','ay','az'))

        definition.define_extra_keywords({'index_of_grid':index_of_grid,'nproc':nproc})

    @property
    def grid(self):
        return self._create_new_grid(self.specify_grid, index_of_grid = 1, nproc=0)



    # Define an object that returns a list of all the blocks in the simulation.
    # This iterates over all processors and then loops over the blocks on the local processors.
    def itergrids(self):
        m = self.get_number_of_procs()

        for x in range(m): # Loop over processors.
            n = self.get_number_of_grids(x)
            #n = max(num_grids)
            #print "N =",n, "X =",x
            for y in range(1, n+1): # Loop over blocks.
                yield self._create_new_grid(self.specify_grid, index_of_grid = y, nproc=x)

    def define_state(self, object):
        CommonCode.define_state(self, object)
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)

        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')

        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')

        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')

        object.add_transition('EDIT', 'RUN', 'initialize_grid')
        object.add_method('RUN', 'evolve_model')
        object.add_method('RUN', 'get_hydro_state_at_point')

        for state in ['EDIT', 'RUN']:
            for methodname in [
#                    'set_particle_pointers',
                    'get_grid_state',
                    'set_grid_state',
                    'set_block_state', # new setter - SCL 11/08/22
#                    'get_potential_at_point',
#                    'get_potential',
#                    'get_gravity_at_point',
#                    'set_potential',
                    'get_grid_density',
                    'set_grid_density',
                    'set_grid_energy_density',
                    'get_grid_energy_density',
                    'get_grid_momentum_density',
                    'set_grid_momentum_density',
                    'get_grid_velocity',
                    'set_grid_velocity',
                    'get_position_of_index',
                    'get_index_of_position',
                    'get_max_refinement',
#                    'set_grid_scalar',
#                    'get_grid_scalar',
                    'get_number_of_grids',
                    'get_index_range_inclusive',
                    'get_cell_volume',
                    'get_number_of_procs',
                    'get_all_local_num_grids',
                    'get_data_all_local_blks',
                    'get_1blk_cell_coords',
                    'get_leaf_indices',
#                    'kick_grid',
#                    'kick_block',
#                    'get_gravity_gas_on_particles',
#                    'get_gravity_particles_on_gas',
                    "grid_update_refinement",
#                    'get_boundary_state',
#                    'set_boundary_state',
#                    'get_boundary_position_if_index',
#                    'get_boundary_index_range_inclusive',
#                    'get_timestep',
#                    'set_timestep',
#                    'get_end_time',
#                    'set_end_time',
#                    'get_time',
#                    'get_max_num_steps',
#                    'set_max_num_steps',
#                    'get_current_step',
#                    'get_restart',
#                    'energy_injection',
#                    'wind_injection',
#                    'get_number_of_particles',
#                    'get_accel_gas_on_particles',
#                    'get_particle_position',
#                    'set_particle_position',
#                    'get_particle_velocity',
#                    'set_particle_velocity',
#                    'get_sink_mean_cs',
#                    'get_sink_gas_mean_velocity',
#                    'get_sink_ang_mom',
#                    'get_sink_gas_var_velocity',
#                    'get_particle_acceleration',
#                    'make_sink',
#                    'add_particles',
#                    'remove_particles',
#                    'remove_all_particles',
#                    'set_starting_local_tag_numbers',
#                    'get_particle_creation_time',
#                    'set_particle_creation_time',
#                    'get_particle_mass',
#                    'set_particle_mass',
#                    'get_particle_oldmass',
#                    'set_particle_oldmass',
#                    'get_particle_nion',
#                    'set_particle_nion',
#                    'get_particle_eion',
#                    'set_particle_eion',
#                    'get_particle_sigh',
#                    'set_particle_sigh',
#                    'set_particle_npep',
#                    'set_particle_epep',
#                    'set_particle_sigd',
#                    'set_particle_gpot',
#                    'get_particle_gpot',
#                    'set_particle_wind_mass',
#                    'set_particle_wind_vel',
#                    'get_particle_tags',
#                    'get_particle_proc',
#                    'get_particle_block',
#                    'get_new_tags',
#                    'get_number_of_new_tags',
#                    'clear_new_tags',
#                    'particles_gather',
#                    'make_stars',
#                    'particles_sort',
                    #'make_particle_tree',
                    'write_chpt',
                    'IO_out',
                    'IO_num',
                    'get_output_dir',
                    'get_runtime_parameter',
                    'timer_summary'
                ]:
                object.add_method(state, methodname)

