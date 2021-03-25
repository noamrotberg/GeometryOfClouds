import numpy as np
import math
import matplotlib.pyplot as plt
import xarray as xr
import cc3d
from mpl_toolkits.mplot3d import axes3d
from grid_specific_functions import *


# constant to transform radian to degree
rad2deg=180.0/np.pi
    
    
def print_component_list (component_list):
    for component in range (len(component_list)):
        print (component + 1, ". component contains following triangles: ")
        print (component_list[component], '\n')


class EnhancedTri:
    def __init__(self, triangle, cube):
        self.triangle = triangle
        self.cube = cube
    def print(self):
        print("Triangle index:", self.triangle, " Cube coordinate:", self.cube)

        
# determine new_edge's colour from the colours of first_known_edge and second_known_edge
def colour_exactly_one_new_edge (first_known_edge, second_known_edge, new_edge, Edge_Colours):
    
    # determine the two known colours in current triangle
    first_colour = Edge_Colours[first_known_edge]
    second_colour = Edge_Colours[second_known_edge]
    
    # determine new colour, e.g. the one not taken from (0,1,2) by first_colour or second_colour
    for i in (0,1,2): 
        if not (first_colour == i or second_colour == i):
            Edge_Colours[new_edge] = i
            return
        
        
# determine colour of new triangle's edges with colours of old triangle's edges
def colour_new_edges (old_triangle_edges, new_triangle_edges, joint_edge, Edge_Colours):
    
    # determine the edges of new triangle that are not shared with old triangle
    new_edges = new_triangle_edges.difference( old_triangle_edges ) 
    first_edge = new_edges.pop()
    second_edge = new_edges.pop()
    
    # are the found edges uncoloured?
    first_is_coloured = first_edge in Edge_Colours
    second_is_coloured = second_edge in Edge_Colours
    
    # case: first_edge is coloured
    if first_is_coloured: 
        if second_is_coloured: 
            return # both known -> nothing to colour
        else: 
            # colour second_edge using first_edge's colour
            colour_exactly_one_new_edge( first_edge, joint_edge, second_edge, Edge_Colours )
    
    # case: second_edge is coloured
    if second_is_coloured: 
        # colour first_edge using second_edge's colour
        colour_exactly_one_new_edge( second_edge, joint_edge, first_edge, Edge_Colours )
    
    # case: both edges are uncoloured
    else: 
        # get old triangle's edges that are not joined with new triangle
        comparison_edges = old_triangle_edges.difference( new_triangle_edges )
        first_comparison_edge = comparison_edges.pop()
        second_comparison_edge = comparison_edges.pop()
        
        # get vertices of first_edge and first_comparison_edge
        first_vertices = get_vertices_of_edge( first_edge )
        first_comparison_vertices = get_vertices_of_edge( first_comparison_edge )
        
        # first_edge and first_comparison_edge have no shared vertices? -> are they parallel?
        parallel = set(first_vertices).isdisjoint( set(first_comparison_vertices) )
        
        if parallel: # parallel -> same colour
            Edge_Colours[ first_edge ] = Edge_Colours[ first_comparison_edge ]
            Edge_Colours[ second_edge ] = Edge_Colours[ second_comparison_edge ] 
        
        else: # not parallel -> first_edge and second_comparison_edge parallel -> interchange colours
            Edge_Colours[ first_edge ] = Edge_Colours[ second_comparison_edge ]
            Edge_Colours[ second_edge ] = Edge_Colours[ first_comparison_edge ]
            
            
# compute cube coordinate of new triangle with old_triangles one and the colour of their joint edge
def determine_cube_coordinates (old, old_triangle_edges, new_triangle_edges, joint_edge, Edge_Colours):
    
    # determine direction in which the two cube coordinates differ
    direction = Edge_Colours[joint_edge] 
    
    # determine direction-vector (1,0,0), (0,1,0) or (0,0,1)
    direction_vector = np.array([0,0,0])
    direction_vector[direction] = 1 
    
    # use invariant: sum of coordinates musst be 0 or 1
    coordinate_sum = sum(old.cube[i] for i in (0,1,2))
    if coordinate_sum == 0:
        # old triangle's coordinate sum is 0 thus new ones has to be 1 -> add
        new_CubeCoord = old.cube + direction_vector
    else: 
        # old triangle's coordinate sum is 1 thus new ones has to be 0 -> subtract
        new_CubeCoord = old.cube - direction_vector
    
    return new_CubeCoord


# compute cube coordinates of outmost triangles
def cubing_next_round (coordinates, visitedTri, outmost, Edge_Colours):
    
    # list for triangles that will be outmost in next iteration
    new_outmost = []
    
    for old in outmost: # consider all at the border of visited
        for neigh in get_neighbors_of_cell( old.triangle ): # consider all neighbours
            new = EnhancedTri(neigh, 0)
            if new.triangle != -9999: # bypass bugs/feature in icon-grid_nawdex_78w40e23n80n_R80000m.nc grid
                if new.triangle not in visitedTri: # new not visited yet

                    # add new to outmost for next round hereafter
                    new_outmost.append( new )

                    # obtain joint edge of new_triangle and old_triangle 
                    old_triangle_edges = set( get_edges_of_cell( old.triangle )) # this is a set for python reasons
                    new_triangle_edges = set( get_edges_of_cell( new.triangle ))
                    joint_edge = ( new_triangle_edges & old_triangle_edges ).pop()

                    # colour all edges of new
                    colour_new_edges (old_triangle_edges, new_triangle_edges, joint_edge, Edge_Colours)

                    # get cube coordinates for new
                    new.cube = determine_cube_coordinates (old, old_triangle_edges, new_triangle_edges, joint_edge, Edge_Colours)
                    visitedTri.append ( new.triangle )
                    coordinates.append(np.array([new.triangle, new.cube]))

            else:
                break
                
    # update outmost triangls
    outmost = new_outmost
    
    return outmost


# shift coordinates such that they are all positive
def shift_coordinates (coordinates, radius):
    shift = int(radius/2)
    
    # update all coordinates
    for entry in coordinates:
        cube_coordinates = entry[1]
        
        # consider all three directions
        for direction in (0,1,2):
            cube_coordinates[direction] += shift # shift coordinates
    
    return coordinates


# determines color for each component such that the different components are distict in plot
def set_color(component):
    
    # list of possible good differentiable colors. Could be changed or new colours added
    color = ['red', 'blue', 'green', 'cyan', 'magenta', 'springgreen', 'yellow', 'indigo', 'darkorange', 
             'olive', 'saddlebrown', 'lightcoral', 'sandybrown', 'fuchsia', 'lightblue', 'teal', 'orangered', 'deepskyblue', 
             'mediumorchid'] #,'chocolate']
    
    # triangle is filled with white if there's no cloud
    if component == 0:
        return 'white'
    
    # more components than colors -> use colors twice or more: modulo calculation
    counter = component%len(color)
    return color[counter]


# function 2
# create 3d array with cloud data
def make_field_array (cube_coordinates, radius, height = None):
     
    # determine size and initialize array
    array_size = radius + 1
    field_array = np.zeros((array_size, array_size, array_size), dtype = 'int')
    
    # write clout data into field_array
    for entry in cube_coordinates:
        triangle = entry[0]     # index of triangle
        cube = entry[1]         # associated cube coordinate
        field = get_field (triangle, height)
        field_array [cube[0]][cube[1]][cube[2]] = field
        
    
    return field_array


# function 3
# if there's a cloud, in which cloud component is it?
def make_connected_components (field_array, connectivity = 'vertex'):
    
    # translate connectivity into input value for 3d connected component labeling
    if connectivity == 'edge':
        connectivity_value = 6
    else:
        connectivity_value = 26     # default
        if connectivity != 'vertex':
            print ("Invalid input: only 'vertex' or 'edge' connectivity are allowed. \
                   Continuing by assuming vertex connectivity.")
            
    # calling external 3d connected component labeling
    component_array = cc3d.connected_components(field_array, connectivity = connectivity_value)
    
    return component_array


# function 4
# which triangle belongs to which component?
def make_list_of_components (cube_coordinates, component_array, print_list=True):
    
    # initialize list for components
    component_list = []
    
    # determine component
    for entry in cube_coordinates:
        triangle = entry[0]    # index of triangle
        cube = entry[1]        # associated cube coordinate
        component = component_array [cube[0]][cube[1]][cube[2]]
        
        # no cloud, no component
        if component != 0:
        
            # append new lists for new component
            while len(component_list) < component:
                component_list.append([])
        
            # add triangle to its component list
            component_list[component - 1].append(triangle)
    
    # show component_list if wanted
    if print_list==True:
        print_component_list (component_list)
    
    return component_list
