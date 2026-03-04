# Converter for Torch-FLASH grid to polaris grid. 
# Written by Stefan Reissl, 2025
# DISCLAIMER: this script is a template converter script made for 
# a specific set of Torch simulations. User discretion advised.

import h5py
import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
from collections import defaultdict
import struct
import sys
import bisect
import random
from functools import reduce
from matplotlib import pyplot as plt
import yt

con_mh = 1.6735575e-25 #hydrogen mass

#octree grid header
grid_id = 20           #grid ID (20 = octree) 
data_ids = [0]         #s. POLARIS manual page 11

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
              
                data_len = len(data_ids)
            
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


def convert(file_input, file_output):
    global cell_counter, nr_of_cells

    ds = yt.load(file_input,unit_system="cgs")
    ad = ds.all_data()

    # List of center positions over all data cells
    # cm
    lvl_idx = np.where(ad['index', 'grid_level'].v > .0)
    
    lst_c_x = ad['index', 'x'][lvl_idx].v
    lst_c_y = ad['index', 'y'][lvl_idx].v
    lst_c_z = ad['index', 'z'][lvl_idx].v
    lst_d_z = ad['index', 'dz'][lvl_idx].v
    lst_ng = ad['flash', 'dens'][lvl_idx].v / con_mh
    
    side_length = lst_c_z.max()-lst_c_z.min() + lst_d_z.min()
    
    lst_level= np.array(np.log(side_length / lst_d_z ) / np.log(2.0),np.int32)

    nr_of_cells=len(lst_ng)
    
    print("cx   : %.6e - %.6e cm"%(lst_c_x.min(),lst_c_x.max()))
    print("cy   : %.6e - %.6e cm"%(lst_c_y.min(),lst_c_y.max()))
    print("cz   : %.6e - %.6e cm"%(lst_c_z.min(),lst_c_z.max()))
    
    print("nd   : %.6e - %.6e cm^-3"%(lst_ng.min(),lst_ng.max()))
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
        
        
        cell = cell_oct(x, y, z, 0, lvl)
        cell.data = [nd]
        
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
        
