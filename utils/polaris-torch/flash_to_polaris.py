# Converter for Torch-FLASH grid to polaris grid. 
# Written by Stefan Reissl, 2025
# Edited by Donglin Wu, July 2026
# DISCLAIMER: this script is a template converter script made for
# a specific set of Torch simulations. User discretion advised.


# IMPORTANT: We assume that the Torch-FLASH grid is an octree grid.

import numpy as np
import struct
import sys

from amuse.lab import SeBa, Particles
from amuse.units import units


import yt
from dust_helper import tau_sp, sublimation_radius


au_to_m = 149597870700.0
con_mh = 1.6735575e-24 #hydrogen mass


#octree grid header
grid_id = 20           #grid ID (20 = octree) 


CLR_LINE = "                                                                  \r"
cell_counter = 0
nr_of_cells = 0

class cell_oct:
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level = _level
    
        self.isleaf = 0
        self.data = []
        self.branches = []      

class cell_oct:
    def __init__(self, _x_min, _y_min, _z_min, _length, _level):
        self.x_min = _x_min
        self.y_min = _y_min
        self.z_min = _z_min
        
        self.length = _length
        self.level = _level
    
        self.isleaf = 0
        self.data = []
        self.branches = []      

class OcTree:
    def __init__(self, _x_min, _y_min, _z_min, _length):
        self.root = cell_oct(_x_min, _y_min, _z_min, _length, 0)

    def initCellBoundaries(self, cell,_level):
        x_min = cell.x_min
        y_min = cell.y_min
        z_min = cell.z_min
        l = 0.5 * cell.length

        level = _level

        cell.isleaf = 0
        cell.data = []
        cell.branches = [None, None, None, None, None, None, None, None]
        cell.branches[0] = cell_oct(x_min, y_min, z_min, l, level)
        cell.branches[1] = cell_oct(x_min + l, y_min, z_min, l, level)
        cell.branches[2] = cell_oct(x_min, y_min + l, z_min, l, level)
        cell.branches[3] = cell_oct(x_min + l, y_min + l, z_min, l, level)

        cell.branches[4] = cell_oct(x_min, y_min, z_min + l, l, level)
        cell.branches[5] = cell_oct(x_min + l, y_min, z_min + l, l, level)
        cell.branches[6] = cell_oct(x_min, y_min + l, z_min + l, l, level)
        cell.branches[7] = cell_oct(x_min + l, y_min + l, z_min + l, l, level)     
        
    def insertInTree(self, cell_pos, cell, _level, _limit):    
        x_pos = cell.x_min
        y_pos = cell.y_min
        z_pos = cell.z_min
        
        if cell_pos.level == cell.level:
            cell_pos.data=cell.data  
            cell_pos.isleaf=1
            
            #print("inserted")
        else:
            if cell_pos.level == _limit:
              
                if len(cell_pos.data)==0:
                    cell_pos.data=[0.0]
            
                d_level=-float(cell_pos.level - cell.level)
                fc=8.0**(-d_level)
              
                data_len = len(self.data_ids)
            
                for i in range(0,data_len):
                    cell_pos.data[i]+=fc*cell.data[i]
              
                cell_pos.isleaf=1
                print("inserted max!")
              
            else:
                
                #print("branch",len(cell_pos.branches))
                
                if len(cell_pos.branches)==0:
                    self.initCellBoundaries(cell_pos,_level+1)
                    
                #print("branch",len(cell_pos.branches))

                x_mid = cell_pos.x_min+0.5*cell_pos.length
                y_mid = cell_pos.y_min+0.5*cell_pos.length
                z_mid = cell_pos.z_min+0.5*cell_pos.length
              
                new_cell_pos = cell_pos
                
                found = False

                if(z_pos < z_mid): #z 0 1 2 3

                    if(y_pos < y_mid): #y 0 1

                        if(x_pos < x_mid): #x 0
                            new_cell_pos = cell_pos.branches[0]
                            found = True
                        else: #x 1
                            new_cell_pos = cell_pos.branches[1]
                            found = True

                    else: #y 2 3

                        if(x_pos < x_mid): #x 2
                            new_cell_pos = cell_pos.branches[2]
                            found = True
                        else: #x 3
                            new_cell_pos = cell_pos.branches[3]
                            found = True

                else: #z 4 5 6 7

                    if(y_pos < y_mid): #y 4 5

                        if(x_pos < x_mid): #x 4
                            new_cell_pos = cell_pos.branches[4]
                            found = True
                        else: #x 5
                            new_cell_pos = cell_pos.branches[5]
                            found = True

                    else: #y 6 7

                        if(x_pos < x_mid): #x 6
                            new_cell_pos = cell_pos.branches[6]
                            found = True
                        else: #x 7
                            new_cell_pos = cell_pos.branches[7]
                            found = True

                if found==False:
                    print("Not found")
                    exit()
                    
                    
                self.insertInTree(new_cell_pos, cell, _level+1,_limit)


    def writeOcTree(self, file, cell):
        global cell_counter
        global nr_of_cells
                       
        file.write(struct.pack("H", cell.isleaf))
        file.write(struct.pack("H", cell.level))   

        if cell.isleaf == 1:    
            data_len = len(cell.data)
            
            if cell_counter % 5000 == 0:
                sys.stdout.write('-> Checking octree integrity : %.3f '%(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1 
         
            for i in range(0, data_len):
                file.write(struct.pack("f", cell.data[i]))
        else:
            
            #print(cell_counter,cell.branches[i])
            
            for i in range(8):
                self.writeOcTree(file, cell.branches[i])
                
                
    def checkOcTree(self, cell):
        global cell_counter
        global nr_of_cells

        if cell.isleaf == 1:    
            length = len(cell.data)
            
            if length == 0:
                print("Wrong data lengths")
                return False
            
            
            if cell_counter % 100 == 0:
                sys.stdout.write('-> Checking octree integrity : %.3f '%(100.0 * cell_counter / nr_of_cells) + ' %     \r')
                sys.stdout.flush()
                
            cell_counter += 1    
            
        else:
            length = len(cell.branches)
            
            if length == 0:
                print("Wrong branch lengths")
                return False
            
            for i in range(8):
                self.checkOcTree(cell.branches[i])                
                
        return True    


# particle_file: path to the Torch particle file
# minimum_mass: minimum mass of stars to consider (in solar masses)
# stellar_metallicity: metallicity of the stars (dimensionless)
def load_stars_from_torch(particle_file, minimum_mass, stellar_metallicity):		
    ds = yt.load(particle_file)
    ad = ds.all_data()

	# Find stars (which has particle_csgm == 0.0) that are more massive than the minimum mass
    star_idx = np.logical_and(ad['all', 'particle_csgm'] == 0.0, (ad['all', 'particle_old_pmass']*yt.units.g).to("Msun").v >= minimum_mass)
    stars = Particles(len(ad['all','particle_mass'][star_idx]))
    stars.initial_mass = (ad['all', 'particle_old_pmass'][star_idx]*yt.units.g).to("Msun").v | units.MSun
    stars.stellar_type = np.ones(len(stars)) | units.stellar_type

    # stellar positions in meters
    stars_x = ad['all', 'particle_posx'][star_idx].to('m').v
    stars_y = ad['all', 'particle_posy'][star_idx].to('m').v
    stars_z = ad['all', 'particle_posz'][star_idx].to('m').v
    
    # Calculate the stellar age in Myr based on the current simulation time and the particle creation time
    t_evol  = ds.current_time.in_units('Myr').v - (ad['all', 'particle_creation_time'][star_idx]*yt.units.s).to('Myr').v

    se = SeBa()
    se.initialize_code()

    # Evolve the star with SeBa based on stellar metallicity and age
    _tmp = se.evolve_star(stars.initial_mass, t_evol | units.Myr, stellar_metallicity) 
    se_time, se_mass, se_radius, se_lum, se_temp, se_evol_time, se_type = _tmp

    # Stellar luminosity in solar luminosities
    stars_L = se_lum.value_in(units.LSun)

    return stars_x, stars_y, stars_z, stars_L



# Create a mask for the cells that are within the sublimation radius of any star
#
# ds_grid: path to the Torch grid data (plot data)
# particle_file, minimum_mass, stellar_metallicity: same as above
# T_sub: sublimation temperature in Kelvin
def mask_sublimation_regions(ds_grid, particle_file, minimum_mass, stellar_metallicity, T_sub):
    ad_grid = ds_grid.all_data()
    stars_x, stars_y, stars_z, stars_L = load_stars_from_torch(particle_file, minimum_mass, stellar_metallicity)

    mask_sublimation = np.zeros(len(ad_grid['index', 'x']), dtype=bool)

    for index_star in range(len(stars_x)):
        star_x = stars_x[index_star]
        star_y = stars_y[index_star]
        star_z = stars_z[index_star]
        star_L = stars_L[index_star]

        r_sub = sublimation_radius(star_L, T_sub) * au_to_m

        # Create a mask for the cells within the sublimation radius
        distance_squared = (ad_grid['index', 'x'].to('m').v - star_x)**2 + (ad_grid['index', 'y'].to('m').v - star_y)**2 + (ad_grid['index', 'z'].to('m').v - star_z)**2
        mask_sublimation_star = distance_squared < r_sub**2

        mask_sublimation = np.logical_or(mask_sublimation, mask_sublimation_star)
    
    return mask_sublimation





# The primary function to convert a Torch-FLASH grid to a POLARIS grid
#
# file_input: path to the Torch-FLASH grid data (plot data)
# file_output: path to the output POLARIS grid file
# mass_fraction: dust-to-gas mass fraction (to save the dust densities into the grid)
#                Using this value is preferred over using dust-to-gas ratio directly in POLARIS, 
#                as it allows for more flexibility in post-processing (e.g. dust sputtering).
# mu: mean molecular weight (default is 1.3)
# timescale_sputtering: threshold for the sputtering timescale in Gyr (default is 1e-2 Gyr). 
#                       If set to 'time_snapshot', it will use the current simulation time as the sputtering timescale.
# clear_sublimation_regions: if True, cells within the sublimation radius of any star will have their dust density set to zero.
# The following parameters are only relevant if clear_sublimation_regions is True:
#   file_particle_input: path to the Torch particle file 
#   minimum_mass_stars: minimum mass of stars to consider for sublimation (in solar masses)
#   stellar_metallicity: metallicity of the stars (dimensionless)
#   T_sub: sublimation temperature in Kelvin (default is 2000 K)

def convert(file_input, file_output, mass_fraction, mu=1.3, 
            timescale_sputtering=1e-2,     # timescale_sputtering in Gyr
            clear_sublimation_regions=False,
            file_particle_input=None,
            minimum_mass_stars=20.0,
            stellar_metallicity=0.014,
            T_sub=2000.0):
    
    global cell_counter, nr_of_cells

    if not mass_fraction is None:
        data_ids = [0, 29, 3]
    else:
        data_ids = [0, 3]


    ds = yt.load(file_input,unit_system="cgs")
    ad = ds.all_data()

    if timescale_sputtering == 'time_snapshot':
        timescale_sputtering = float(ds.current_time.in_units('Gyr'))
    if clear_sublimation_regions:
        if file_particle_input is None:
            raise ValueError("file_particle_input must be provided when clear_sublimation_regions is True.")
        
        mask_sublimation = mask_sublimation_regions(ds, file_particle_input, minimum_mass_stars, stellar_metallicity, T_sub)




    # List of center positions over all data cells
    # cm
    lvl_idx = np.where(ad['index', 'grid_level'].v > .0)
    
    lst_c_x = ad['index', 'x'][lvl_idx].v
    lst_c_y = ad['index', 'y'][lvl_idx].v
    lst_c_z = ad['index', 'z'][lvl_idx].v
    lst_d_z = ad['index', 'dz'][lvl_idx].v
    lst_ng = ad['flash', 'dens'][lvl_idx].v / con_mh / mu
    lst_rhog = ad['flash', 'dens'][lvl_idx].v
    lst_temp = ad['flash', 'temp'][lvl_idx].v

    # lst_dust_temp = ad['flash', 'tdus'][lvl_idx].v
    
    side_length = lst_c_z.max()-lst_c_z.min() + lst_d_z.min()
    
    lst_level= np.array(np.log(side_length / lst_d_z ) / np.log(2.0),np.int32)

    nr_of_cells=len(lst_ng)

    if clear_sublimation_regions:
        assert len(mask_sublimation) == nr_of_cells
    
    print("cx   : %.6e - %.6e cm"%(lst_c_x.min(),lst_c_x.max()))
    print("cy   : %.6e - %.6e cm"%(lst_c_y.min(),lst_c_y.max()))
    print("cz   : %.6e - %.6e cm"%(lst_c_z.min(),lst_c_z.max()))
    
    print("nd   : %.6e - %.6e cm^-3"%(lst_ng.min(),lst_ng.max()))
    # print("dust temp : %.6e - %.6e K"%(lst_dust_temp.min(),lst_dust_temp.max()))
    print("temp : %.6e - %.6e K"%(lst_temp.min(),lst_temp.max()))
    print("level: %.d - %.d "%(lst_level.min(),lst_level.max()))
    print("sl   : %.6e   cm"%side_length)
    print("cells: %d   "%nr_of_cells)
    
    tree = OcTree(-0.5*side_length, -0.5*side_length, -0.5*side_length, side_length)    
    
    
    cell_root = tree.root
    
    for i in range(nr_of_cells):
        
        if i%5000==0:
            sys.stdout.write('-> Inserting cells :  %.3f '%(100.0 * i / nr_of_cells) + ' %     \r')
            
        
        x=lst_c_x[i]
        y=lst_c_y[i]
        z=lst_c_z[i]
        nd=lst_ng[i]
        lvl=lst_level[i]
        
        # temp_dust = lst_dust_temp[i]

        
        cell = cell_oct(x, y, z, 0, lvl)
        if mass_fraction is not None:
            #The following divison by 1000 is necessary because POLARIS converts the length/volume part, but doesn't automatically convert g to kg. 
            rhog = lst_rhog[i]/1000  #in kg cm-3 
            
            # Set the dust density given a dust-to-gas ratio
            rhod = mass_fraction*rhog #in kg cm-3
            temp=lst_temp[i]
            
            # If the sputtering timescale is less than the threshold, set the dust density to effectively zero. 
            # In practice, 1e-30 times the gas density is used. This is done to avoid numerical issues in POLARIS.
            # tau_sp returns sputtering timescale in Gyr, takes density in g cm-3 and temperature in K
            if tau_sp(lst_rhog[i], temp) < timescale_sputtering: 
                rhod = 1e-30*rhog
                print(f'Dust in cell {i} sputtered.')

            # If the cell is within the sublimation radius of any star, set the dust density to effectively zero.
            elif clear_sublimation_regions and mask_sublimation[i]:
                rhod = 1e-30*rhog
                print(f'Dust in cell {i} sublimated.')

            cell.data = [nd, rhod, temp]
        else:
            cell.data = [nd, temp]
        
        tree.insertInTree(cell_root, cell,0,100)
       
    sys.stdout.write(CLR_LINE)    
    print("Constructing octree:    done   ")

    #check octree integrity
    cell_counter=0
    check = tree.checkOcTree(cell_root)
    
    sys.stdout.write(CLR_LINE)  
    if check == False:
        print("ERROR: Octree integrity is inconsistent!   \n\n")
        exit ()
    else:
        print("Octree structure   :    OK      ")
        
        
    #write octree file header
    data_len = len(data_ids)
    file = open(file_output, "wb")
        
    file.write(struct.pack("H", grid_id))
    file.write(struct.pack("H", data_len))

    for d_ids in data_ids:
        file.write(struct.pack("H", d_ids))

    file.write(struct.pack("d", side_length))
    
    #write octree
    cell_counter = 0.0
    tree.writeOcTree(file, tree.root)
    sys.stdout.write(CLR_LINE)

    print("Writing octree     :    done   \n")
    
    print("Octree successfully created")
        
