##### import modules #####
#
from grid_specific_functions import *


##### helpfunctions for make_cube_coordinates #####

# class to enhance a triangle with its cube coordinate:
#   for a given triangle create an object consisting of the triangle index and the triangle's cube coordinate
class EnhancedTri:
    def __init__(self, triangle, cube):
        self.triangle = triangle
        self.cube = cube
    def print(self):
        print("Triangle index:", self.triangle, " Cube coordinate:", self.cube)

        

# colouring edges of a triangle
#   situation:  the colours of two edges are already known
#               the third edge must be coloured
#   called from colour_new_edges in this special situation
# INPUT:
#   first_known_edge      int     index of a already coloured edge
#   second_known_edge     int     index of a already coloured edge
#   new_edge              int     index of the not yet coloured edge
#   edge_colours          dict    dictionary storing to each edge index the colour
# OUTPUT:
#   None, but updates edge_colours

def colour_exactly_one_new_edge (first_known_edge, second_known_edge, new_edge, edge_colours):
    
    # get the colours of the already coloured edges
    first_colour = edge_colours[first_known_edge]
    second_colour = edge_colours[second_known_edge]
    
    # determine new colour, i.e. the one not taken from (0,1,2) by first_colour and second_colour
    for i in (0,1,2): 
        if not (first_colour == i or second_colour == i):
            edge_colours[new_edge] = i  # colour found -> update edge_colours
            return
        
        
# colouring edges of a triangle
#   from old_triangle with already coloured edges determine the edge_colours of the adjacent new_triangle
# INPUT:
#   old_triangle_edges      set     contains the edge indices of old_triangle's edges
#   new_triangle_edges      set     contains the edge indices of new_triangle's edges
#   joint_edge              int     index of the shared edge from old_triangle and new_triangle
#   edge_colours            dict    dictionary storing to each edge index the colour
# OUTPUT:
#   None, but updates edge_colours

def colour_new_edges (old_triangle_edges, new_triangle_edges, joint_edge, edge_colours):
    
    # determine the edges of new_triangle that are not shared with old_triangle
    new_edges = new_triangle_edges.difference( old_triangle_edges ) 
    first_edge = new_edges.pop()
    second_edge = new_edges.pop()
    
    # are the found edges uncoloured?
    first_is_coloured = first_edge in edge_colours
    second_is_coloured = second_edge in edge_colours
    
    # case: first_edge is coloured
    if first_is_coloured: 
        if second_is_coloured: 
            return # both known -> nothing to colour
        else: 
            # colour second_edge using the colours of first_edge and joint_edge
            colour_exactly_one_new_edge( first_edge, joint_edge, second_edge, edge_colours )
    
    # case: second_edge is coloured
    if second_is_coloured: 
        # colour first_edge using the colours of second_edge and joint_edge
        colour_exactly_one_new_edge( second_edge, joint_edge, first_edge, edge_colours )
    
    # case: both edges are uncoloured
    else: 
        # get old_triangle's edges that are not joined with new_triangle
        comparison_edges = old_triangle_edges.difference( new_triangle_edges )
        first_comparison_edge = comparison_edges.pop()
        second_comparison_edge = comparison_edges.pop()
        
        # get vertices of first_edge and first_comparison_edge
        first_vertices = get_vertices_of_edge( first_edge )
        first_comparison_vertices = get_vertices_of_edge( first_comparison_edge )
        
        # if first_edge and first_comparison_edge have no shared vertices:
        #   the two edges must be parallel and thus, have the same colour
        #   else, first_edge must be parallel to second_comparison_edge and thus, have the same colour
        parallel = set(first_vertices).isdisjoint( set(first_comparison_vertices) )
        
        if parallel: # parallel: same colour
            edge_colours[ first_edge ] = edge_colours[ first_comparison_edge ]
            edge_colours[ second_edge ] = edge_colours[ second_comparison_edge ] 
        
        else: # not parallel: first_edge and second_comparison_edge are parallel and same-coloured
            edge_colours[ first_edge ] = edge_colours[ second_comparison_edge ]
            edge_colours[ second_edge ] = edge_colours[ first_comparison_edge ]
            
         
        
# compute cube coordinate of a new triangle:
#   from the old_triangle's cube_coordinate 
#       derive the cube_coordinate of the adjacent new_triangle 
#       by use of the colour of their joint_edge 
# INPUT:
#   old                     EnhancedTri     the old EnhancedTri with known cube_coordinates
#   joint_edge              int             index of the shared edge from old_triangle and new_triangle
#   edge_colours            dict            dictionary storing to each edge index the colour
# OUTPUT:
#   new_cube_coordinate, the cube_coordinate of new_triangle

def determine_cube_coordinates (old, joint_edge, edge_colours):
    
    # get colour of the joint_edge and call it direction
    direction = edge_colours[joint_edge]
    # the direction encodes the parallel class of the joint edge
    # this encodes also the direction in which one has to walk from old_triangle to new_triangle
    # this gives the coordinate in which the cube_coordinates of old_- and new_triangle differ
    
    # determine direction-vector (1,0,0), (0,1,0) or (0,0,1) by which the cube_coordinates will differ
    direction_vector = np.array([0,0,0])
    direction_vector[direction] = 1 
    
    # use invariant: sum of coordinates musst be 0 or 1
    # determine wether the direction vector has to be added or subtracted
    # by using the invariant that the sum of all valid cube_coordinates has to be 0 or 1
    old_coordinate_sum = sum(old.cube[i] for i in (0,1,2))
    if old_coordinate_sum == 0:
        # old triangle's coordinate sum is 0 thus new ones has to be 1 -> add
        new_cube_coordinate = old.cube + direction_vector
    else: 
        # old triangle's coordinate sum is 1 thus new ones has to be 0 -> subtract
        new_cube_coordinate = old.cube - direction_vector
    
    return new_cube_coordinate



# iterate over the outmost triangles and compute the cube_coordinates of their adjacent triangles
# INPUT:
#   cube_coordinates        array       EnhancedTri's of the triangles whose cube_coordinates are already computed
#   visited_triangles       list        indices of the triangles whose cube_coordinates are already computed
#   outmost                 list        EnhancedTri's of the triangles considered in the round before   
#   edge_colours            dict        dictionary storing to each edge index the colour
# OUTPUT:
#   updated outmost, EnhancedTri's of the new triangles considered this round

def cubing_next_round (cube_coordinates, visited_triangles, outmost, edge_colours):
    
    # list for triangles that will be outmost in next iteration
    new_outmost = []
    
    for old in outmost: # consider all EnhancedTri's at the border of visited
        for neigh in get_neighbors_of_cell( old.triangle ): # consider all neighbours
            new = EnhancedTri(neigh, 0)
            if new.triangle == -9999: # use feature: triangles at the grid border claim to be adjacent to -9999
                break
            else:
                if new.triangle not in visited_triangles: # EnhancedTri new has no cube_coordinate yet

                    # add new to new_outmost for next round hereafter and update visited_triangles
                    new_outmost.append( new )
                    visited_triangles.append ( new.triangle )

                    # preparation: obtain edges of new_triangle and old_triangle and derive their joint_edge
                    old_triangle_edges = set( get_edges_of_cell( old.triangle )) # this is a set for python reasons...
                    new_triangle_edges = set( get_edges_of_cell( new.triangle ))
                    joint_edge = ( new_triangle_edges & old_triangle_edges ).pop()

                    # colour the edges of new.triangle
                    colour_new_edges (old_triangle_edges, new_triangle_edges, joint_edge, edge_colours)

                    # get cube coordinates for EnhancedTri new and update the array of cube_coordinates
                    new.cube = determine_cube_coordinates (old, old_triangle_edges, new_triangle_edges, joint_edge, edge_colours)
                    cube_coordinates.append(np.array([new.triangle, new.cube]))
    
    return new_outmost



# shift coordinates such that they are all positive
# INPUT:
#   cube_coordinates    array   stores triangle - cube_coordinate pairs
#   radius              int     same as in make_cube_coordiantes
# OUTPUT:
#   shifted cube_coordinates

def shift_coordinates (cube_coordinates, radius):
    shift = int(radius/2)
    
    # update all coordinates
    for entry in cube_coordinates:
        cube = entry[1] # get cube_coordinate 
        
        # consider all three directions
        for direction in (0,1,2):
            cube[direction] += shift # shift coordinate
    
    return cube_coordinates



##### helpfunctions for make_connected_2/3d_components #####

# create 3d array containing the field data
# INPUT:
#   cube_coordinates    array   output of make_cube_coordinates   
#   radius              int     same as in make_cube_coordiantes
#   height              int     only used when considering only one height layer of a 3d grid
# OUTPUT:
#   field_array, 3d matrix containing the field data
def make_field_array (cube_coordinates, radius, height = None):
     
    # determine size and initialize array
    array_size = radius + 1
    field_array = np.zeros((array_size, array_size, array_size), dtype = 'int')
    
    # write the field data into field_array
    for entry in cube_coordinates:
        triangle = entry[0]     # index of triangle
        cube = entry[1]         # cube coordinate
        field = get_field (triangle, height) # get the field data 
        field_array [cube[0]][cube[1]][cube[2]] = field
        
    return field_array



# compute connected components
# INPUT:
#   field_array    array    output of make_field_array   
#   connectivity   string   use 'vertex' (default) or 'edge' connectivity
# OUTPUT:
#   component_array, 3d matrix containing the component indices
def make_connected_components (field_array, connectivity = 'vertex'):
    
    # translate connectivity - string into input value for 3d connected component labeling
    if connectivity == 'edge':
        connectivity_value = 6
    else:
        connectivity_value = 26     
        if connectivity != 'vertex':    # use default 'vertex' connectivity for invalid input
            print ("Invalid input: only 'vertex' or 'edge' connectivity are allowed. \
                   Continuing by assuming vertex connectivity.")
            
    # calling external 3d connected component labeling
    component_array = cc3d.connected_components(field_array, connectivity = connectivity_value)
    
    return component_array



# takes a list of components and prints for each component 
#   the component number followed by its triangle indices
def print_component_list (component_list):
    for component in range (len(component_list)):
        print (component + 1, ". component contains following triangles: ")
        print (component_list[component], '\n')
        
        

# which triangle belongs to which component?
# INPUT:
#   cube_coordinates    array       output of make_cube_coordinates   
#   component_array     string      output of make_connected_components
#   print_list          Boolean     print output (default True)
# OUTPUT:
#   Returns a list of lists. 
#   An inner list correspond to a connected component and contains all belonging triangle-indices.

def make_list_of_components (cube_coordinates, component_array, print_list=True):
    
    # initialize list for components
    component_list = []
    
    # determine component
    for entry in cube_coordinates:
        triangle = entry[0]    # index of triangle
        cube = entry[1]        # cube coordinate
        component = component_array [cube[0]][cube[1]][cube[2]]
        
        if component != 0:  # triangle lies in a component
        
            # append new lists for new component if necessary
            while len(component_list) < component:
                component_list.append([])
        
            # add triangle to its component list
            component_list[component - 1].append(triangle)
    
    # show component_list if wanted
    if print_list==True:
        print_component_list (component_list)
    
    return component_list



##### helpfunctions for plot_connected_2/3d_components #####

# asings to a given component a colour to distinguish different components in the plot
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
