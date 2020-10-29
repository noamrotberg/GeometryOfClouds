import numpy as np
import math
import matplotlib.pyplot as plt
import xarray as xr
import sys
sys.path.append('/pf/b/b380459/python-3.6-anaconda3-bleedingedge')
import cc3d
from mpl_toolkits.mplot3d import axes3d


rad2deg=180.0/np.pi


def path_input(path, file, datafile):
    file = path+file
    
    global ds
    global df 
    
    ds = xr.open_dataset(file) 
    
    df = xr.open_dataset(path+datafile)
        


def get_neighbors_of_cell (triangle):
    neighbors = ds.neighbor_cell_index[:,triangle]-1
    return neighbors.values

def get_vertices_of_cell (triangle):
    vertices = ds.vertex_of_cell[:, triangle ]-1
    return vertices.values

def get_edges_of_cell (triangle):
    edges = ds.edge_of_cell[:, triangle ]-1
    return edges.values

def get_vertices_of_edge (edge):
    vertices = ds.edge_vertices[:, edge ]-1
    return vertices.values

def print_3d_array (array):
    n = array.shape
    for i in range(n[0]):
        for j in range(n[1]):
            print(array[i][j])
        print()
        
def print_4d_array (array):
    n=array.shape
    länge = len(n)
    if länge == 3:
        print_3d_array(array)
    else:
        for i in range(n[0]):
            for j in range(n[1]):
                for k in range(n[2]):
                    print(array[i][j][k])
        print()
        
def get_cloud (triangle, height, rounding = 0.6):
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
        for neigh in get_neighbors_of_cell( old.triangle ): # considere all neighbours
            new = EnhancedTri(neigh, 0)
            if new.triangle != -1: # bypass bugs in icon-grid_nawdex_78w40e23n80n_R80000m.nc grid
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
                    new.cube = determine_cube_coordinates (old, old_triangle_edges, new_triangle_edges, joint_edge,
                                                           Edge_Colours)
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
             'olive', 'saddlebrown', 'lightcoral', 'sandybrown', 'fuchsia', 'lightblue', 'teal']
    
    # triangle is filled with white if there's no cloud
    if component == 0:
        return 'white'
    
    # more components than colors -> use colors twice or more: modulo calculation
    counter = component%len(color)
    return color[counter]
    


#########################################################################
# function 1
# Determine cube coordinates of triangles
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
        
        
    
    return shift_coordinates(coordinates, radius)
    #return coordinates


# function 2
# create 3d array with cloud data

def make_cloud_datafield (cube_coordinates, radius, height = None, print_datafield = False):
     
    # determine size and initialize array
    array_size = radius + 1
    cloud_datafield = np.zeros((array_size, array_size, array_size), dtype = 'int')
    
    # write clout data into cloud_datafield
    for entry in cube_coordinates:
        triangle = entry[0]     # index of triangle
        cube = entry[1]         # associated cube coordinate
        cloud = get_cloud (triangle, height)
        cloud_datafield [cube[0]][cube[1]][cube[2]] = cloud
        
    # show cloud_datafiled if wanted
    if print_datafield:
        print_3d_array(cloud_datafield)
    
    return cloud_datafield



# function 3
# if there's a cloud, in which cloud component is it?
def make_connected_components (cloud_datafield, connectivity = 'vertex', print_datafield = False):
    
    # translate connectivity into input value for 3d connected component labeling
    if connectivity == 'edge':
        connectivity_value = 6
    else:
        connectivity_value = 26     # default
        if connectivity != 'vertex':
            print ("Error in input: only 'vertex' or 'edge' connectivity are allowed. \
                   Assuming vertex connectivity.")
            
    
    # calling external 3d connected component labeling
    component_datafield = cc3d.connected_components(cloud_datafield, connectivity = connectivity_value)

    
    # show component_datafiled if wanted
    if print_datafield:
        print_xd_array(component_datafield)
    
    return component_datafield



# function 4
# which triangle belongs to which component?
def make_list_of_components (component_datafield, cube_coordinates, print_list = True):
    
    # initialize list for components
    component_list = []
    
    # determine component
    for entry in cube_coordinates:
        triangle = entry[0]     # index of triangle
        cube = entry[1]        # associated cube coordinate
        component = component_datafield [cube[0]][cube[1]][cube[2]]
        
        # no cloud, no component
        if component != 0:
        
            # append new lists for new component
            while len(component_list) < component:
                component_list.append([])
        
            # add triangle to its component list
            component_list[component - 1].append(triangle)
    
    # show component_list if wanted
    if print_list:
        print_component_list (component_list)
    
    return component_list
   



# function 5
# compute connected components for data of several altitudes/heights
def make_connected_3d_components (cube_coordinates, radius, 
                                  connectivity = "vertex", plotting_result = False, save_plot = False):
    
    # initialize 4d-matrix
    layers = len(df.height[:]) # number of layers
    array_size = radius + 1
    four_d_matrix = np.zeros((layers, array_size, array_size, array_size), dtype = 'int')
    
    # initialize component list
    component_list=[]
    
    # consider each height
    for height in range(layers):
        # compute connected components in each height
        cloud_datafield= make_cloud_datafield (cube_coordinates, radius, height, print_datafield = False)
        component_datafield = make_connected_components (cloud_datafield, connectivity = connectivity, print_datafield = False)
        height_component_list = make_list_of_components (component_datafield, cube_coordinates, print_list = False)
        
        
        # write cloud data into 4d-matrix
        four_d_matrix [height,:,:,:] = cloud_datafield[:,:,:]
    
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
                                
    print_component_list(component_list)
    
    
    # show plot of final components if wanted
    if plotting_result == True:
        plotting_3d(component_list, save_plot)
        
    return component_list


# function 6.1
# highlighting components in grid
def plotting(component_datafield, cube_coordinates, 
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
        component = component_datafield [cube[0]][cube[1]][cube[2]]
        vertex = np.array(get_vertices_of_cell(triangle))
        
        # arrays with longitude and latitude for one triangle
        # vertex 0 is used twice for a closed line around the triangle
        vlon = np.array([rad2deg*ds.vlon[vertex[0]],rad2deg*ds.vlon[vertex[1]],
                         rad2deg*ds.vlon[vertex[2]],rad2deg*ds.vlon[vertex[0]]])
        vlat = np.array([rad2deg*ds.vlat[vertex[0]],rad2deg*ds.vlat[vertex[1]],
                         rad2deg*ds.vlat[vertex[2]],rad2deg*ds.vlat[vertex[0]]])

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
def plotting_3d (component_list,  save_name = 'connected_components_3d.pdf', 
                 plot_title = 'connected components 3d', save_plot = False):
    # initialization
    vertex = []
    vlat = []
    vlon = []
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.set_xlabel('X Axes')
    ax.set_ylabel('Y Axes')
    ax.set_zlabel('Z Axes')

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
            vlon = np.array([rad2deg*ds.vlon[vertex[0]],rad2deg*ds.vlon[vertex[1]],
                             rad2deg*ds.vlon[vertex[2]],rad2deg*ds.vlon[vertex[0]]])
            vlat = np.array([rad2deg*ds.vlat[vertex[0]],rad2deg*ds.vlat[vertex[1]],
                             rad2deg*ds.vlat[vertex[2]],rad2deg*ds.vlat[vertex[0]]])

            
            # plotting the lines between the vertices
            # determines color of the triangle line such that every triangle belonging to the same component has the same color
            ax.plot3D(vlon,vlat,height, color = set_color(comp))



    # size of the showed picture    
    plt.rcParams["figure.figsize"] = [60,20]

    # title and axis labels
    ax.set_title(plot_title, fontsize=30)
    plt.xlabel('longitude',fontsize=20); plt.ylabel('latitude',fontsize=20)
 
    # save plot with name <save_name> if wanted
    if save_plot:
        plt.savefig(save_name)

    plt.show()