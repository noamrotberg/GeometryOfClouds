##### import modules #####

import numpy as np
import math
import matplotlib.pyplot as plt
import xarray as xr
import sys
sys.path.append('/pf/b/b380906/python-3.6-anaconda3-bleedingedge')
import cc3d
from mpl_toolkits.mplot3d import axes3d

from helpfunctions import *
from outsourcing import *

    
##### main functions #####

# function 1
# determine cube coordinates of triangles
def make_cube_coordinates (start_triangle, radius, print_progress=False): 
    
    # initialization
    init = EnhancedTri(start_triangle, np.array([0,0,0])) # assign coordinates (0,0,0) to start_triangle
    coordinates = [np.array([init.triangle, init.cube])]
    visitedTri = [init.triangle]
    outmost = [init]
    
    # colour the edges of start_triangle with 0,1,2
    Edge_Colours = { get_edges_of_cell(start_triangle)[x] : x for x in (0,1,2) }
    
    # iterating
    for n in range (radius):
        
        # update outmost triangles
        outmost = cubing_next_round (coordinates, visitedTri, outmost, Edge_Colours)
        
        # show progress if print_info is True
        if print_progress:
            print('Round ', n+1, 'has finished. Visited ', len(coordinates), \
                  'triangles, thereof ', len(outmost), 'new.')
    
    # return coordinates after shifting them to positivity
    return shift_coordinates(coordinates, radius)
    




##########################################################
def make_connected_2d_components (cube_coordinates, radius, height = None, connectivity = 'vertex', 
                                  make_plot = False, save_plot = False, save_name = 'connected_components.pdf', 
                                  plot_title = 'connected components on icon'): 
    
    field_array = make_field_array (cube_coordinates, radius, height)
    component_array = make_connected_components (field_array, connectivity)
    component_list = make_list_of_components (cube_coordinates, component_array)
    if make_plot == True:
        plotting(cube_coordinates, component_array, save_plot, save_name, plot_title)
    
    return component_list

# function 5
# compute connected components for data of several altitudes/heights
def make_connected_3d_components (cube_coordinates, radius, connectivity = 'vertex', make_plot = False, save_plot = False, 
                                  save_name = 'connected_components.pdf', plot_title = 'connected components on icon'):
    
    # initialize layers
    layers = len(get_height()) # number of layers
    
    # initialize component list
    component_list=[]
    
    # consider each height
    for height in range(layers):
        # compute connected components in each height
        field_array= make_field_array (cube_coordinates, radius, height)
        component_array = make_connected_components (field_array, connectivity = connectivity)
        height_component_list = make_list_of_components (cube_coordinates, component_array, print_list = False)
    
        # update horizontal components
        for h_component in height_component_list:
            
            # append h_component to component_list
            component_list.append ([])   
            for triangle in h_component:      
                component_list[len(component_list)-1].append ((triangle, height))
            
            # compare with layer height-1
            if height > 0:
                # initialize new_list as h_component
                new_list = component_list[len(component_list)-1]
                
                # merge new_list with triangles of height-1
                for triangle in h_component:
                    for component in component_list:
                        
                        # if there's a cloud in height-1
                        if (triangle, height-1) in component:
                            if (triangle, height-1) not in new_list:
                                
                                # update component_list and new_list
                                # VERBESSERUNGSWÜRDIG
                                component.extend(new_list)
                                component_list.remove(new_list)
                                new_list = component
                                break
                                
    print_component_list(component_list)
      
    # show plot of final components if wanted
    if make_plot == True:
        plotting_3d(component_list,  save_name, plot_title, save_plot)
        
    return component_list


# function 6.1
# highlighting components in grid
def plotting(cube_coordinates, component_array, 
             save_plot = False, save_name = 'connected_components.pdf', plot_title = 'connected components on icon'):
    
    # initialization
    vertex = []
    vlat = []
    vlon = []
    fig = plt.figure()
    ax = fig.gca()
    
    # determine longitute and latitute for each vertex of every triangle
    for entry in cube_coordinates:
        triangle = entry[0]
        cube = entry[1]
        component = component_array [cube[0]][cube[1]][cube[2]]
        vertex = np.array(get_vertices_of_cell(triangle))
        
        # arrays with longitude and latitude for one triangle
        # vertex 0 is used twice for a closed line around the triangle
        vlon = np.array([rad2deg*get_vlon(vertex[0]),rad2deg*get_vlon(vertex[1]),
                         rad2deg*get_vlon(vertex[2]),rad2deg*get_vlon(vertex[0])])
        vlat = np.array([rad2deg*get_vlat(vertex[0]),rad2deg*get_vlat(vertex[1]),
                         rad2deg*get_vlat(vertex[2]),rad2deg*get_vlat(vertex[0])])

        # determines color of the filled triangle such that every triangle belonging to the same component has the same color
        ax.fill(vlon,vlat,color=set_color(component))
        
        # plotting the lines between the vertices
        ax.plot(vlon,vlat,color='black', linestyle='-')
        
    # size of the showed picture    
    plt.rcParams["figure.figsize"] = [60,20]
    
    # title and axis labels
    ax.set_title(plot_title, fontsize=30)
    plt.xlabel('longitude',fontsize=20); plt.ylabel('latitude',fontsize=20)
    
    # save plot with name <save_name> if wanted 
    if save_plot:
        plt.savefig(save_name)
    
    plt.show()


# function 6.2
# highlighting components in grid
def plotting_3d (component_list,  save_name = 'connected_components_3d.pdf', plot_title = 'connected components 3d', save_plot = False):
    
    # initialization
    vertex = []
    vlat = []
    vlon = []
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.set_xlabel('longitude', fontsize=25)
    ax.set_ylabel('latitude',fontsize=25)
    ax.set_zlabel('height',fontsize=25)

    # for each componet plot all it's (triangle, height)-tuples
    comp = 0
    for component in component_list:
        comp += 1
        for entry in component:
            triangle = entry[0]
            height = np.full(4,[entry[1]])
            
            vertex = np.array(get_vertices_of_cell(triangle))

            # arrays with longitude and latitude for one triangle
            # vertex 0 is used twice for a closed line around the triangle
            vlon = np.array([rad2deg*get_vlon(vertex[0]),rad2deg*get_vlon(vertex[1]),
                             rad2deg*get_vlon(vertex[2]),rad2deg*get_vlon(vertex[0])])
            vlat = np.array([rad2deg*get_vlat(vertex[0]),rad2deg*get_vlat(vertex[1]),
                             rad2deg*get_vlat(vertex[2]),rad2deg*get_vlat(vertex[0])])

            # plotting the lines between the vertices
            # determines color of the triangle line such that every triangle belonging to the same component has the same color
            ax.plot3D(vlon,vlat,height, color = set_color(comp))

    # size of the showed picture    
    plt.rcParams["figure.figsize"] = [60,20]

    # title and axis labels
    ax.set_title(plot_title, fontsize=30)
    #plt.xlabel('longitude',fontsize=15); plt.ylabel('latitude',fontsize=15); plt.zlabel('height', frontsize=15)
 
    # save plot with name <save_name> if wanted
    if save_plot:
        plt.savefig(save_name)

    plt.show()