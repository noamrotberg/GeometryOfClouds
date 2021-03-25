import sys
sys.path.append('/pf/b/b380906/python-3.6-anaconda3-bleedingedge')
import numpy as np
import math
import matplotlib.pyplot as plt
import xarray as xr
import cc3d
from mpl_toolkits.mplot3d import axes3d

def path_input(path, gridfile, datafile):
    gridfile = path+gridfile
    global ds
    global df 
    
    ds = xr.open_dataset(gridfile) 
    
    df = xr.open_dataset(path+datafile)

def get_neighbors_of_cell (triangle):
    neighbors = ds.neighbor_cell_index[:,triangle]-1
    neighbors = neighbors.values
    for triangle in range(len(neighbors)):
        if neighbors[triangle] == -1:
            neighbors[triangle] = -9999
            
            
    return neighbors

def get_vertices_of_cell (triangle):
    vertices = ds.vertex_of_cell[:, triangle ]-1
    return vertices.values


def get_edges_of_cell (triangle):
    edges = ds.edge_of_cell[:, triangle ]-1
    return edges.values


def get_vertices_of_edge (edge):
    vertices = ds.edge_vertices[:, edge ]-1
    return vertices.values


def get_field (triangle, height, rounding = 0.6):
    if height == None:
        cloud = df.clct[0, triangle].values
    else:
        cloud = df.clc[0, height, triangle].values
        
    cloud = cloud/100
    
    if cloud >= rounding:
        cloud = 1
    elif cloud < rounding:
        cloud = 0
    
    return cloud

def get_vlon(vertex):
    vlon = ds.vlon[vertex]
    return vlon.values

def get_vlat(vertex):
    vlat = ds.vlat[vertex]
    return vlat.values

def get_height():
    height = df.height[:]
    return height