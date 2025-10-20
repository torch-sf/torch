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
flux = generic_unit_system.energy / time / length ** 2

class FlashInterface(CodeInterface, HydrodynamicsInterface):

    use_modules = ['flash_run']

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="flash_worker", **keyword_arguments)


    @legacy_function
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_particle_pointers():
        function = LegacyFunctionSpecification()
        function.addParameter('part_type_in', dtype='s', direction=function.IN, default='mass')
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

    @legacy_function
    def set_block_state():
        # VorAMR addition - SCL
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('procID', dtype='i', direction=function.IN)
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
        function.must_handle_array = True
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
        function.must_handle_array = True
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
    def get_grid_flux_photoelectric():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('flux_pe', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_grid_flux_ionizing():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['i','j','k', 'index_of_grid','nproc']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('flux_ion', dtype='d', direction=function.OUT)
        function.addParameter('ngridpoints', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    @legacy_function
    def get_potential():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        for x in ['i','j','k', 'index_of_grid']:
            function.addParameter(x, dtype='i', direction=function.IN)
        function.addParameter('potential', dtype='d', direction=function.OUT)
        function.result_type='int32'
        return function

    @legacy_function
    def get_potential_at_point():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['eps','x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type='int32'
        return function

    #@legacy_function
    #def get_gravity_at_point():
    #    function = LegacyFunctionSpecification()
    #    function.must_handle_array = True
    #    for x in ['eps','x','y','z']:
    #        function.addParameter(x, dtype='d', direction=function.IN)
    #    for x in ['gax','gay','gaz']:
    #        function.addParameter(x, dtype='d', direction=function.OUT)
    #    function.addParameter('nparts',dtype='i',direction=function.LENGTH)
    #    function.result_type='int32'
    #    return function

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
    def set_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def get_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_end_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN, default=0.0)
        function.result_type='i'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def get_end_time():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_max_num_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def get_max_num_steps():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def get_current_step():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_begin_iter_step():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.IN, default=1)
        function.result_type='i'
        return function

    @legacy_function
    def get_begin_iter_step():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def initialize_restart():
        function = LegacyFunctionSpecification()
        function.result_type='i'
        return function

    @legacy_function
    def get_restart():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='b', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def set_restart():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='b', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def grid_update_refinement():
        function = LegacyFunctionSpecification()
        function.addParameter('gridChanged', dtype='b', direction=function.OUT)
        function.result_type='i'
        return function

    @legacy_function
    def get_hydro_state_at_point():
        function = LegacyFunctionSpecification()
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, default = 0)
        for x in ['rho','rhovx','rhovy','rhovz','rhoen']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
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
        function.addParameter('procID', dtype='i', direction=function.IN)
        function.addParameter('limits', dtype='i', direction=function.IN)
        function.addParameter('coords', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def kick_grid():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def kick_block():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['accel_x', 'accel_y', 'accel_z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('blockID', dtype='i', direction=function.IN)
        function.addParameter('block_arr', dtype='i', direction=function.IN)
        function.addParameter('limits', dtype='i', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type='i'
        return function

    @legacy_function
    def get_gravity_gas_on_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.addParameter('kick_number', dtype='i', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def get_gravity_particles_on_gas():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.addParameter('kick_number', dtype='i', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def energy_injection():
        function = LegacyFunctionSpecification()
        #function.must_handle_array = True
        for x in ['energy', 'fracKin', 'mass','xloc','yloc','zloc']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('dt', dtype='d', direction=function.OUT)
        #function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    #@legacy_function
    #def wind_injection():
        #function = LegacyFunctionSpecification()
        ##function.must_handle_array = True
        #for x in ['starMass', 'injectMass','xloc','yloc','zloc']:
            #function.addParameter(x, dtype='d', direction=function.IN)
        #function.addParameter('dt', dtype='d', direction=function.OUT)
        ##function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        #function.result_type = 'i'
        #return function

########################################
# Particle Stuff
########################################

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_accel_gas_on_particles():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('eps', dtype='d', direction=function.IN)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_num_part_prop():
        function = LegacyFunctionSpecification()
        function.addParameter('n', dtype='int32', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_position_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN,unit=NO_UNIT)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=length)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_position():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=length)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_velocity():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_velocity_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_sink_mean_vel_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_sink_var_vel_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=speed**2.0)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_acceleration_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['ax','ay','az']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=acc)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_sink_ang_mom_array():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.OUT, unit=length*mass*speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_ang_mom():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN, unit=NO_UNIT)
        for x in ['lx','ly','lz']:
            function.addParameter(x, dtype='d', direction=function.IN, unit=speed)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def make_sink():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('tags', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def add_particles():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        for x in ['x','y','z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.addParameter('tags', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def remove_particles():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type='i'
        return function

    @legacy_function
    def remove_all_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('clear_local_tag', dtype='b', default='False', direction=function.IN)
        function.result_type='i'
        return function

    @legacy_function
    def set_starting_local_tag_numbers():
        function = LegacyFunctionSpecification()
        function.result_type='i'
        return function

    @legacy_function
    def set_particle_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_oldmass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_nion():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('nion', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_eion():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('eion', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_rel_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('rel_mass', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_rel_age():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('rel_age', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_corem():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('corem', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_co_corem():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('co_corem', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_particle_stype():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('stype', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_particle_radius():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_sigh():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('sigh', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_npep():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('nion', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_epep():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('eion', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_sigd():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('sigh', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_wind_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('dmdt', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_wind_vel():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('velw', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_sink_mean_cs():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('cs', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_nion():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('nion', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_eion():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('eion', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_rel_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('rel_mass', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_rel_age():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('rel_age', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_corem():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('corem', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_co_corem():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('co_corem', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
    
    
    @legacy_function
    def get_particle_stype():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('stype', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_particle_radius():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('radius', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_sigh():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('sigh', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_mass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_oldmass():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_creation_time():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('creation_time', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_creation_time():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('creation_time', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_gpot():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.OUT)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_particle_gpot():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('gpot', dtype='d', direction=function.IN)
        function.addParameter('nparts',dtype='i',direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_tags():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('n', dtype='int32', direction=function.IN)
        function.addParameter('tags', dtype='d', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_proc():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('procs', dtype='int32', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_particle_block():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('tags', dtype='d', direction=function.IN)
        function.addParameter('blocks', dtype='int32', direction=function.OUT)
        function.addParameter('nparts', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_new_tags():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('new_tags_length', dtype = 'i', direction=function.IN)
        function.addParameter('tags', dtype = 'd', direction=function.OUT)
        function.addParameter('nparts', dtype = 'i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_number_of_new_tags():
        function = LegacyFunctionSpecification()
        function.addParameter('new_tag_num', dtype = 'int32', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def clear_new_tags():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def particles_gather():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def make_stars():
        function = LegacyFunctionSpecification()
        function.addParameter('dt', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def particles_sort():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

#    @legacy_function
#    def make_particle_tree():
#        function = LegacyFunctionSpecification()
#        function.result_type = 'i'
#        return function


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

###############################################
# Default implemenation made by build.py - Josh
###############################################

#class Flash(InCodeComponentImplementation):

#    def __init__(self, unit_converter = None, **options):
#        InCodeComponentImplementation.__init__(self,  FlashInterface(**options), **options)

#        object.set_converter(self.unit_converter.as_converter_from_si_to_generic())

#####################################################
# Attempt to copy amrvac's class implemenation, working so far! - Josh
#####################################################

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
    def get_particle_position(self, tags):

        [x, y, z] = self.get_particle_position_array(tags)

        pos_array = np.array([x.value_in(units.cm), y.value_in(units.cm),
                  z.value_in(units.cm)]).transpose() | units.cm

        if (hasattr(x,"__len__") == False):

            pos_array = pos_array.flatten()

        return pos_array

    def get_particle_velocity(self, tags):

        [x, y, z] = self.get_particle_velocity_array(tags)

        vel_array = np.array([x.value_in(units.cm / units.s), y.value_in(units.cm / units.s),
                  z.value_in(units.cm / units.s)]).transpose() | units.cm / units.s

        if (hasattr(x,"__len__") == False):

            vel_array = vel_array.flatten()

        return vel_array

    def get_particle_acceleration(self, tags):

        [x, y, z] = self.get_particle_acceleration_array(tags)

        acc_array = np.array([x.value_in(units.cm / (units.s**2.0)), y.value_in(units.cm / (units.s**2.0)),
                  z.value_in(units.cm / (units.s**2.0))]).transpose() | units.cm / (units.s**2.0)

        if (hasattr(x,"__len__") == False):

            acc_array = acc_array.flatten()

        return acc_array

    def get_sink_gas_mean_velocity(self, tags):

        [x, y, z] = self.get_sink_mean_vel_array(tags)

        vel_array = np.array([x.value_in(units.cm / units.s), y.value_in(units.cm / units.s),
                  z.value_in(units.cm / units.s)]).transpose() | units.cm / units.s

        if (hasattr(x,"__len__") == False):

            vel_array = vel_array.flatten()

        return vel_array

    def get_sink_gas_var_velocity(self, tags):

        [x, y, z] = self.get_sink_var_vel_array(tags)

        vel_array = np.array([x.value_in(units.cm**2.0 / units.s**2.0), y.value_in(units.cm**2.0 / units.s**2.0),
                  z.value_in(units.cm**2.0 / units.s**2.0)]).transpose() | units.cm**2.0 / units.s**2.0

        if (hasattr(x,"__len__") == False):

            vel_array = vel_array.flatten()

        return vel_array

    def get_sink_ang_mom(self, tags):

        [lx, ly, lz] = self.get_sink_ang_mom_array(tags)

        #print "In get_sink_ang_mom"
        #print lx, ly, lz

        ang_mom_array = np.array([lx.value_in(units.cm**2.0 * units.g / units.s), ly.value_in(units.cm**2.0 * units.g / units.s),
                  lz.value_in(units.cm**2.0 * units.g / units.s)]).transpose() | units.cm**2.0 * units.g / units.s

        #print ang_mom_array

        #print "array has len attr?", hasattr(lx,"__len__")

        if (hasattr(lx,"__len__") == False):

            ang_mom_array = ang_mom_array.flatten()

        #print ang_mom_array
        #print "Leaving get_sink_ang_mom"

        return ang_mom_array

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
            'set_particle_pointers',
            object.NO_UNIT,
            object.ERROR_CODE
        )

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
        # VorAMR addition - SCL
        object.add_method(
            'set_block_state',
            (object.INDEX, object.INDEX, object.INDEX,
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
            'get_potential',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (potential, object.ERROR_CODE,)
        )

        object.add_method(
            'get_potential_at_point',
            (length, length, length, length),
            (potential, object.ERROR_CODE,)
        )

        #object.add_method(
        #    'get_gravity_at_point',
        #    (length, length, length, length),
        #    (acc, acc, acc, object.ERROR_CODE,)
        #)

        object.add_method(
            'get_cell_volume',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (length**3, object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_flux_photoelectric',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (flux, object.ERROR_CODE,)
        )

        object.add_method(
            'get_grid_flux_ionizing',
            (object.INDEX, object.INDEX, object.INDEX, object.INDEX, object.INDEX),
            (flux, object.ERROR_CODE,)
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
            (object.NO_UNIT,object.NO_UNIT,object.NO_UNIT,object.NO_UNIT),
            (length, object.ERROR_CODE)
        )

        object.add_method(
            'get_leaf_indices',
            (object.NO_UNIT),
            (object.NO_UNIT, object.INDEX, object.INDEX, object.ERROR_CODE)
        )

        object.add_method(
            'kick_grid',
            (time),
            (object.ERROR_CODE,)
        )
        object.add_method(
            'kick_block',
            (acc, acc, acc, object.INDEX, object.INDEX, object.NO_UNIT, time),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_gravity_gas_on_particles',
            (time, object.NO_UNIT),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_gravity_particles_on_gas',
            (time, object.NO_UNIT),
            (object.ERROR_CODE)
        )

        object.add_method(
            "grid_update_refinement",
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            "get_timestep",
            (),
            (time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_timestep",
            (time),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_end_time",
            (),
            (time, object.ERROR_CODE,)
        )

        object.add_method(
            "set_end_time",
            (time, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_time',
            (),
            (time, object.ERROR_CODE,)
        )

        object.add_method(
            'get_max_num_steps',
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            'set_max_num_steps',
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_current_step',
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            'get_restart',
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )

        object.add_method(
            'set_restart',
            (object.NO_UNIT,),
            (object.ERROR_CODE,)
        )

        object.add_method(
            'get_hydro_state_at_point',
            (length, length, length,
                speed, speed, speed),
            (density, momentum, momentum,
                momentum, energy, object.ERROR_CODE)
        )

        object.add_method(
            'energy_injection',
            (generic_unit_system.energy, object.NO_UNIT,
             mass, length, length, length),
            (time, object.ERROR_CODE)
        )

        #object.add_method(
            #'wind_injection',
            #(mass, mass, length, length, length),
            #(time, object.ERROR_CODE)
        #)

###################################
######### Particles
###################################

        object.add_method(
            'get_number_of_particles',
            (),
            (object.NO_UNIT,object.ERROR_CODE,)
        )

        object.add_method(
            'get_accel_gas_on_particles',
            (length,length,length,length),
            (acc,acc,acc, object.ERROR_CODE,)
        )

        object.add_method(
            'get_num_part_prop',
            (),
            (object.NO_UNIT,object.ERROR_CODE)
        )

        ### I'm implementing this with my own defined function
        ### so that the structure of the array return looks
        ### right.

        #object.add_method(
            #'get_particle_position',
            #(object.NO_UNIT),
            #(length, length, length, object.ERROR_CODE)
        #)

        object.add_method(
            'set_particle_position',
            (object.NO_UNIT, length, length, length),
            (object.ERROR_CODE)
        )

        ### Same as above.

        #object.add_method(
            #'get_particle_velocity',
            #(object.NO_UNIT),
            #(speed, speed, speed, object.ERROR_CODE)
        #)

        object.add_method(
            'set_particle_velocity',
            (object.NO_UNIT, speed, speed, speed),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_ang_mom',
            (object.NO_UNIT,
             units.cm**2.0 * units.g / units.s,
             units.cm**2.0 * units.g / units.s,
             units.cm**2.0 * units.g / units.s),
            (object.ERROR_CODE)
        )

        object.add_method(
            'make_sink',
            (length, length, length),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'add_particles',
            (length, length, length),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'remove_particles',
            (object.NO_UNIT),
            (object.ERROR_CODE)
        )

        object.add_method(
            'remove_all_particles',
            (object.NO_UNIT),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_starting_local_tag_numbers',
            (),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_mass',
            (object.INDEX),
            (mass, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_oldmass',
            (object.INDEX),
            (mass, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_creation_time',
            (object.INDEX),
            (time, object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_creation_time',
            (object.INDEX, time),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_nion',
            (object.INDEX),
            ((time**-1.0), object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_eion',
            (object.INDEX),
            (mass*length**2.0*(time**-2.0), object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_rel_mass',
            (object.INDEX), 
            (mass, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_rel_age',
            (object.INDEX), 
            (time, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_corem',
            (object.INDEX), 
            (mass, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_co_corem',
            (object.INDEX), 
            (mass, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_stype',
            (object.INDEX), 
            (object.NO_UNIT, object.ERROR_CODE)
        )
        
        object.add_method(
            'get_particle_radius',
            (object.INDEX), 
            (length, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_sigh',
            (object.INDEX),
            (length*length, object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_mass',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_oldmass',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_nion',
            (object.INDEX, (time**-1.0)),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_eion',
            (object.INDEX, mass*length**2.0*(time**-2.0)),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_rel_mass',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_rel_age',
            (object.INDEX, time),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_corem',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_co_corem',
            (object.INDEX, mass),
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_stype',
            (object.INDEX, object.NO_UNIT), 
            (object.ERROR_CODE)
        )
        
        object.add_method(
            'set_particle_radius',
            (object.INDEX, length), 
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_sigh',
            (object.INDEX, length*length),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_npep',
            (object.INDEX, (time**-1.0)),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_epep',
            (object.INDEX, mass*length**2.0*(time**-2.0)),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_sigd',
            (object.INDEX, length*length),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_wind_mass',
            (object.INDEX, mass/time),
            (object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_wind_vel',
            (object.INDEX, length/time),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_sink_mean_cs',
            (object.INDEX),
            (length/time, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_gpot',
            (object.INDEX),
            (potential, object.ERROR_CODE)
        )

        object.add_method(
            'set_particle_gpot',
            (object.NO_UNIT, potential),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_tags',
            (object.INDEX),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_proc',
            (object.INDEX),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_particle_block',
            (object.INDEX),
            (object.INDEX, object.ERROR_CODE)
        )

        object.add_method(
            'write_chpt',
            (),
            (object.ERROR_CODE)
        )

        object.add_method(
            'IO_out',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'IO_num',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_output_dir',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_runtime_parameter',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'timer_summary',
            (),
            (object.ERROR_CODE)
        )

        object.add_method(
            'get_new_tags',
            (object.NO_UNIT),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'get_number_of_new_tags',
            (),
            (object.NO_UNIT, object.ERROR_CODE)
        )

        object.add_method(
            'clear_new_tags',
            (),
            (object.ERROR_CODE)
        )

        object.add_method(
            'particles_gather',
            (),
            (object.ERROR_CODE)
        )

        object.add_method(
            'make_stars',
            (time),
            (object.ERROR_CODE)
        )

        object.add_method(
            'particles_sort',
            (),
            (object.ERROR_CODE)
        )

#        object.add_method(
#            'make_particle_tree',
#            (),
#            (object.ERROR_CODE)
#        )

    def specify_grid(self, definition, index_of_grid = 1, nproc=0):
        definition.set_grid_range('get_index_range_inclusive')

        definition.add_getter('get_position_of_index', names=('x','y','z'))

        definition.add_getter('get_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_grid_state', names=('rho', 'rhovx','rhovy','rhovz','energy'))
        definition.add_setter('set_block_state', names=('rho', 'vx', 'vy', 'vz', 'energy', 'gpot')) # VorAMR addition - SCL
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

        definition.add_getter('get_grid_flux_photoelectric', names=('flux_photoelectric',))
        definition.add_getter('get_grid_flux_ionizing', names=('flux_ionizing',))


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


    #def define_particle_sets(self, object):
        #object.define_set('particles', 'index_of_the_particle')
        #object.set_new('particles', 'new_particle')
        #object.set_delete('particles', 'delete_particle')
        #object.add_setter('particles', 'set_state')
        #object.add_getter('particles', 'get_state')
        #object.add_setter('particles', 'set_mass')
        #object.add_getter('particles', 'get_mass', names = ('mass',))
        #object.add_setter('particles', 'set_position')
        #object.add_getter('particles', 'get_position')
        #object.add_setter('particles', 'set_velocity')
        #object.add_getter('particles', 'get_velocity')
        #object.add_setter('particles', 'set_radius')
        #object.add_getter('particles', 'get_radius')
        #object.add_query('particles', 'get_indices_of_colliding_particles', public_name = 'select_colliding_particles')

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
                    'set_particle_pointers',
                    'get_grid_state',
                    'set_grid_state',
                    'set_block_state', # VorAMR addition - SCL
                    'get_potential_at_point',
                    'get_potential',
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
                    'get_grid_flux_photoelectric',
                    'get_grid_flux_ionizing',
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
                    'kick_grid',
                    'kick_block',
                    'get_gravity_gas_on_particles',
                    'get_gravity_particles_on_gas',
                    "grid_update_refinement",
#                    'get_boundary_state',
#                    'set_boundary_state',
#                    'get_boundary_position_if_index',
#                    'get_boundary_index_range_inclusive',
                    'get_timestep',
                    'set_timestep',
                    'get_end_time',
                    'set_end_time',
                    'get_time',
                    'get_max_num_steps',
                    'set_max_num_steps',
                    'get_current_step',
                    'get_restart',
                    'energy_injection',
#                    'wind_injection',
                    'get_number_of_particles',
                    'get_accel_gas_on_particles',
                    'get_particle_position',
                    'set_particle_position',
                    'get_particle_velocity',
                    'set_particle_velocity',
                    'get_sink_mean_cs',
                    'get_sink_gas_mean_velocity',
                    'get_sink_ang_mom',
                    'get_sink_gas_var_velocity',
                    'get_particle_acceleration',
                    'make_sink',
                    'add_particles',
                    'remove_particles',
                    'remove_all_particles',
                    'set_starting_local_tag_numbers',
                    'get_particle_creation_time',
                    'set_particle_creation_time',
                    'get_particle_mass',
                    'set_particle_mass',
                    'get_particle_oldmass',
                    'set_particle_oldmass',
                    'get_particle_nion',
                    'set_particle_nion',
                    'get_particle_eion',
                    'set_particle_eion',
                    'get_particle_rel_mass',
                    'set_particle_rel_mass',
                    'get_particle_rel_age',
                    'set_particle_rel_age',
                    'get_particle_corem',
                    'set_particle_corem',
                    'get_particle_co_corem',
                    'set_particle_co_corem',
                    'get_particle_stype',
                    'set_particle_stype',
                    'get_particle_radius',
                    'set_particle_radius',
                    'get_particle_sigh',
                    'set_particle_sigh',
                    'set_particle_npep',
                    'set_particle_epep',
                    'set_particle_sigd',
                    'set_particle_gpot',
                    'get_particle_gpot',
                    'set_particle_wind_mass',
                    'set_particle_wind_vel',
                    'get_particle_tags',
                    'get_particle_proc',
                    'get_particle_block',
                    'get_new_tags',
                    'get_number_of_new_tags',
                    'clear_new_tags',
                    'particles_gather',
                    'make_stars',
                    'particles_sort',
                    #'make_particle_tree',
                    'write_chpt',
                    'IO_out',
                    'IO_num',
                    'get_output_dir',
                    'get_runtime_parameter',
                    'timer_summary'
                ]:
                object.add_method(state, methodname)

