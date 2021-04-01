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

# 
from grid_specific_functions import * # 
from helpfunctions import * ## rename


    
##### main functions #####

# make_cube_coordinates:
#   Assign to each triangle a 3d-coordinate called cube coordinate.
# INPUT: 
#   start_triangle      int         the triangle index where the algorithm starts
#   radius              int         the radius where the iteration terminates
#   print_progress      Boolean     print progress of iteration (default False)
# OUTPUT: 
#   Returns a list of arrays of the following form:
#   np.array([<triangle index>, <cube array>])
#   where <cube array> is an array with three entries encoding the cube coordinate of the triangle with <trinagle index>. 

def make_cube_coordinates (start_triangle, 
                           radius, 
                           print_progress=False): 
    
    # initialization
    init = EnhancedTri(start_triangle, np.array([0,0,0])) # assign coordinates (0,0,0) to start_triangle
    cube_coordinates = [np.array([init.triangle, init.cube])]
    visited_triangles = [init.triangle]
    outmost = [init] # list for EnhancedTri that were added last round and thus lie "at the boundary of discovery"
    
    # colour the edges of start_triangle with 0,1,2
    edge_colours = { get_edges_of_cell(start_triangle)[x] : x for x in (0,1,2) }
    
    # iterating
    for n in range (radius):
        
        # update outmost triangles
        outmost = cubing_next_round (cube_coordinates, visited_triangles, outmost, edge_colours)
        
        # show progress if print_progress is True
        if print_progress:
            print('Round ', n+1, 'has finished. Visited ', len(cube_coordinates), \
                  'triangles, thereof ', len(outmost), 'new.')
    
    # return coordinates after shifting them to positivity
    return shift_coordinates(cube_coordinates, radius)
    

    
# make_connected_2d_components:
#   Directly from cube coordinates, compute components-list containing triangle indices.
# INPUT:
#   cube_coordinates    list of arrays      output of make_cube_coordinates 
#   radius              int                 same as in make_cube_coordiantes
#   height              int                 only used for 3d grids
#   connectivity        string              use 'vertex' (default) or 'edge' connectivity
#   make_plot           Boolean             plot the components (default False)
#   save_plot           Boolean             save plot (default False)
#   save_name           string              name for save-file
#   plot_title          string              title for plot
# OUTPUT: 
#   Returns a list of lists. 
#   An inner list correspond to a connected component and contains all belonging triangle-indices.

def make_connected_2d_components (cube_coordinates, 
                                  radius, 
                                  height = None, 
                                  connectivity = 'vertex', 
                                  make_plot = False, 
                                  save_plot = False, 
                                  save_name = 'connected_components.pdf', 
                                  plot_title = 'connected components on icon'): 
    
    # generates a 3-dimensional array containing the field:
    # if a triangle is covered, the field entry at the triangles cube coordinate is set to 1
    field_array = make_field_array (cube_coordinates, radius, height)
    
    # applies cc3d-module to compute connected components on field_array:
    # yields a new 3-dimensional array 
    # containing the connected component index of a triangle with cube coordinates (x,y,z) at position component_array[x][y][z]
    component_array = make_connected_components (field_array, connectivity)
    
    # a component list of triangle indices is generated for each connected component:
    # returns a list containing all component lists
    component_list = make_list_of_components (cube_coordinates, component_array)
   
    # show plot if wanted
    if make_plot == True:
        plot_connected_2d_components (cube_coordinates, component_array, save_plot, save_name, plot_title)
    
    return component_list



# make_connected_3d_components:
#   For each height layer: Compute local connected component with functions make_connected_2d_components on this layer. 
#   Append the found components to (global) component_list by component-wise insertion of tuples (triangle index, height). 
#   Merge components if one contains the tuple (triangle index, height) and the other (triangle index, height-1) and update the list in this way.
# INPUT:
#   cube_coordinates    list of arrays      output of make_cube_coordinates
#   radius              int                 same as in make_cube_coordinates
#   connectivity        string              use 'vertex' (default) or 'edge' connectivity
#   make_plot           Boolean             plot the components (default False)
#   save_plot           Boolean             save plot (default False)
#   save_name           string              name for save-file
#   plot_title          string              title for plot
# OUTPUT:
#   Returns a list of lists.  An inner list correspond to a connected component and contains all belonging (triangle index, height)-tuples

def make_connected_3d_components (cube_coordinates, 
                                  radius, 
                                  connectivity = 'vertex', 
                                  make_plot = False, 
                                  save_plot = False, 
                                  save_name = 'connected_components.pdf', 
                                  plot_title = 'connected components on icon'):
    
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
                                component.extend(new_list)
                                component_list.remove(new_list)
                                new_list = component
                                break
    # show comonent_list                            
    print_component_list (component_list)
      
    # show plot of final components if wanted
    if make_plot == True:
        plot_connected_3d_components (component_list,  save_name, plot_title, save_plot)
        
    return component_list



# plot_connected_2d_components:
#   Show plot of coloured components in 3d-grid
# INPUT: 
#   component_list    list of lists   output of make_connected_2d_components
#   save_name         string          name of save-file
#   plot_title        string          title of the plot
#   save_plot         Boolean         save plot (default False)

def plot_connected_2d_components (cube_coordinates, 
                                  component_array, 
                                  save_plot = False, 
                                  save_name = 'connected_components.pdf', 
                                  plot_title = 'connected components on icon'):
    
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


# plot_connected_3d_components
#   Show plot of coloured components in 3d-grid
# INPUT: 
#   component_list    list of lists   output of make_connected_3d_components
#   save_name         string          name of save-file
#   plot_title        string          title of the plot
#   save_plot         Boolean         save plot (default False)

def plot_connected_3d_components (component_list,  save_name = 'connected_components_3d.pdf', plot_title = 'connected components 3d', save_plot = False):
    
    # initialization
    vertex = []
    vlat = []
    vlon = []
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    
    # label axes
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

    # title
    ax.set_title(plot_title, fontsize=30)
     
    # save plot with name <save_name> if wanted
    if save_plot:
        plt.savefig(save_name)

    plt.show()
