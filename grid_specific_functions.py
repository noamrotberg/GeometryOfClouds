##### import modules #####
# standard
import numpy as np
import math
import xarray as xr

# for connected components module
import sys
sys.path.append('/pf/b/b380906/python-3.6-anaconda3-bleedingedge')
import cc3d # connected components algorithm

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



##### functions #####

# path_input:
#   create globale variables to open the xArray dataset from gridfile and datafile
# INPUT:
#   path        string      location where gridfile and datafile are stored
#   gridfile    string      name of the file containing the grid
#   datafile    string      name of the file containing the field

def path_input(path, gridfile, datafile):
    
    global ds
    ds = xr.open_dataset(path + gridfile)
    
    global df 
    df = xr.open_dataset(path + datafile)

    
    
# get_neighbors_of_cell:
#   returns edge neighbors of a given triangle

def get_neighbors_of_cell (triangle):
    
    neighbors = ds.neighbor_cell_index[:,triangle]-1 # index shift: Python starts to count from 0, grid from 1
    neighbors = neighbors.values
    
    for triangle in range(len(neighbors)):
        if neighbors[triangle] == -1:   # triangles at the grid border claim to be adjacent to -1
            neighbors[triangle] = -9999 # set to -9999
            
    return neighbors



# get_vertices_of_cell:
#   returns the three vertices of a given triangle

def get_vertices_of_cell (triangle):
    vertices = ds.vertex_of_cell[:, triangle ]-1    # index shift: Python starts to count from 0, grid from 1
    return vertices.values



# get_edges_of_cell:
#   returns the three edges of a given triangle

def get_edges_of_cell (triangle):
    edges = ds.edge_of_cell[:, triangle ]-1 # index shift: Python starts to count from 0, grid from 1
    return edges.values



# get_vertices_of_edge:
#   returns the two vertices of a given edge

def get_vertices_of_edge (edge):
    vertices = ds.edge_vertices[:, edge ]-1 # index shift: Python starts to count from 0, grid from 1
    return vertices.values



# get_field:
#   returns whether a given triangle is covered by a cloud (in our case) or not
# INPUT:
#   triangle    int     the triangle index
#   height      int     only used for 3d grids, specifies the considered height-layer
#   rounding    int     threshold for cloud coverage to be considered as covered/ not covered
# OUTPUT:
#   cloud       0/1     0 if the triangle is not covered (by a cloud), 1 if it is 
def get_field (triangle, height, rounding = 0.6):
    
    if height == None:  # 2d case
        cloud = df.clct[0, triangle].values # getting cloud data of triangle 
    else: # 3d case
        cloud = df.clc[0, height, triangle].values # getting cloud data of (triangle, height)-tuple
        
    cloud = cloud/100   # rescale percent to decimal number
    
    if cloud >= rounding:
        cloud = 1
    elif cloud < rounding:
        cloud = 0
    
    return cloud



# get_vlon:
#   returns the longitude of a given vertex

def get_vlon(vertex):
    vlon = ds.vlon[vertex]
    return vlon.values



# get_vlat:
#   returns the latitude of a given vertex

def get_vlat(vertex):
    vlat = ds.vlat[vertex]
    return vlat.values



# get_height:
#   returns the number of height-layers for the grid

def get_height():
    height = df.height[:]
    return height
