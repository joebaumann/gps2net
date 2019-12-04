from shapely.geometry import LineString
from shapely.geometry import Point, LineString
import fiona
import networkx as nx
import math

import matplotlib.pyplot as plt
import numpy as np

from timeit import default_timer as timer
import doctest

import sys
import os


# global variable: empty Directed Graph
DG = nx.DiGraph()


def blockPrint():
    '''This method is used to disable print() messages.
    '''
    sys.stdout = open(os.devnull, 'w')


def enablePrint():
    '''This method restores print() messages. 
    '''
    sys.stdout = sys.__stdout__


# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd='\r'):
    r'''Call in a loop to create terminal progress bar.

    Parameters
    ----------
    iteration : int
        current iteration
    total : Int
        total iterations
    prefix : str, optional
        prefix string, by default ''
    suffix : str, optional
        suffix string, by default ''
    decimals : int, optional
        positive number of decimals in percent complete, by default 1
    length : int, optional
        character length of bar, by default 100
    fill : str, optional
        bar fill character, by default '█'
    printEnd : str, optional
        end character (e.g. '\\r', '\\r\\n'), by default '\\r'

    '''

    percent = ('{0:.' + str(decimals) + 'f}').format(100 *
                                                     (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)

    # enablePrint()

    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

    # blockPrint()


def distFrom(lng1, lat1, lng2, lat2):
    '''Returns the distance between two points in meters.

    Parameters
    ----------
    lng1 : float
        This float is the longitude of the source.
    lat1 : float
        This float is the latitude of the source.
    lng2 : float
        This float is the longitude of the target.
    lat2 : float
        This float is the latitude of the target.


    Notes
    -----
    This function calculates the distance between two gps points.
    The result is not 100% correct as the function does not consider the elipsis-like shape of the earth. Instead it just uses an earth radius of 6371000 meters for the calculation.

    A more accurate calculation of the distance can be done in QGIS. For more details, please visit: http://www.qgistutorials.com/en/docs/calculating_line_lengths.html
    However, this function was used to be able to run the code independently of QGIS.

    Examples
    --------

    >>> x=12
    >>> x
    12
    >>> (5<10)
    True
    >>> myDist = distFrom(-122.115, 37.115, -122.111, 37.111)
    >>> myDist
    568.8872918546489
    >>> distFrom(-122.115, 37.115, -122.111, 37.111) # doctest: +SKIP
    5


    Returns
    --------
    dist : int
        The distance between two points in meters.

    '''

    earthRadius = 6371000  # meters
    dLat = math.radians(lat2-lat1)
    dLng = math.radians(lng2-lng1)
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(math.radians(lat1)) * \
        math.cos(math.radians(lat2)) * \
        math.sin(dLng/2) * math.sin(dLng/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    dist = earthRadius * c

    return dist


def cut(line, distance, point):
    '''Cuts a line in two at a distance from its starting point.

    Parameters
    ----------
    line : Shapely 'LineString' object
        [description]
    distance : float
        [description]
    point : Shapely Point object
        [description]

    Notes
    -----
    ...
    'line.interpolate(distance)' should be equal to the parameter 'point'. However, due to precision errors 'point' is used.
    ...

    Returns
    -------
    list with LineStrings
        Two LineStrings are returned. 'point' is the end point of the first LineString and the start point of the second LineString.
    '''

    # check if point lies on line. If not, return uncut line.
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)

    # loop through all all points (nodes) of the linestring
    for i, p in enumerate(coords):
        # i is the index
        # p is the point on the line
        # pd is the distance of point p on the line
        pd = line.project(Point(p))
        # pd==distance means that point where the line should be cut is already a point on the line at index i (= so this means that there is a linesegment in the line which has the point as a starting or end point)
        if (pd == distance):
            # return two lines (whithout cutting any line segment)
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        # pd>distance means that a line segment of the line has to be cut
        # If the line is a circle (start of the line equals the end of the line) and the line has to be cut in the last line segment 'pd>distance' is not enough as pd=0.0 for the end point.
        if (pd > distance or (i != 0 and pd == 0.0 and coords[0] == coords[len(coords)-1])):
            # return the cut line
            return [
                LineString(coords[:i] + [(point.x, point.y)]),
                LineString([(point.x, point.y)] + coords[i:])]


def createGraphFromSHPInput(filepath_shp):
    '''Creates a directed graph from a shp file.

    Parameters
    ----------
    filepath_shp : str
        The path where the shp file (which contains the street data) is stored.


    Notes
    -----
    This graph is only initialized once in the beginning of the script. The DiGraph contains all linesegment of the shp file as directed edges. The start and end node of each edge are tuples containing the coordinates of the position. The weight of an edge is the air_line_distance from the start to the end of the linesegment. Each edge containes the following attributes:

    - id : int
        The id of the street (LineString) which the linesegment belongs to according to the shapefile.
    - oneway : {'B', 'F', 'T'}
        The 'oneway'-property of the street (LineString) which the linesegment belongs to according to the shapefile. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).


    Returns
    -------
    Directed Graph
        The Directed Graph which contains all linesegments of the shapefiles as edges is returned.
    '''

    GraphFromSHP = nx.DiGraph()

    # blockPrint()

    counter = 0

    print('Graph already exists')

    nr_elements_in_SHP_file = 0
    with fiona.open(filepath_shp) as input_lines:
        nr_elements_in_SHP_file = sum(1 for line in input_lines)

    print('nr_elements_in_SHP_file: ', nr_elements_in_SHP_file)

    # Initial call to print 0% progress ProgressBar
    suffix = '| SHP file items: {}/{}'.format(
        counter+1, nr_elements_in_SHP_file)
    printProgressBar(0, nr_elements_in_SHP_file,
                     prefix='Creating Graph:', suffix=suffix, length=50)

    # print('graph nr of edges: ', GraphFromSHP.number_of_edges())
    print('starting to build graph: ')

    with fiona.open(filepath_shp) as input_lines:

        for line in list(input_lines):
            counter_inner = 0

            previous_coord = (0, 0)
            id = line['id']
            oneway = line['properties']['oneway']

            # print('id: ', id)
            # print('oneway: ', oneway)

            for coord in line['geometry']['coordinates']:

                if (previous_coord != (0, 0)):
                    # print('previous_coord: ', previous_coord)
                    # print('coord: ', coord)

                    calculated_length = distFrom(
                        coord[0], coord[1], previous_coord[0], previous_coord[1])

                    if(oneway == 'B'):
                        # create two directed edges for both directions
                        GraphFromSHP.add_edge(
                            coord, previous_coord, weight=calculated_length, id=id, oneway=oneway)
                        GraphFromSHP.add_edge(
                            previous_coord, coord, weight=calculated_length, id=id, oneway=oneway)
                    elif(oneway == 'F'):
                        # create a directed edge from the from-node to the to-node
                        GraphFromSHP.add_edge(
                            previous_coord, coord, weight=calculated_length, id=id, oneway=oneway)
                    elif(oneway == 'T'):
                        # create a directed edge from the to-node to the from-node
                        GraphFromSHP.add_edge(
                            coord, previous_coord, weight=calculated_length, id=id, oneway=oneway)

                # else:
                    # print('revious_coord    ==    (0,0)   --------------')

                # print(counter_inner)
                counter_inner += 1
                previous_coord = coord

            # enablePrint()

            # Update Progress Bar
            suffix = '| SHP file items: {}/{}'.format(
                counter+1, nr_elements_in_SHP_file)
            printProgressBar(counter + 1, nr_elements_in_SHP_file,
                             prefix='Creating Graph:', suffix=suffix, length=50)

            # blockPrint()

            counter += 1

    # enablePrint()

    print('graph nr of edges AFTER CREATION: ', GraphFromSHP.number_of_edges())

    # blockPrint()

    # print('GraphFromSHP: ', (GraphFromSHP.edges(data=True)))

    # enablePrint()

    # print('graph nr of edges: ', GraphFromSHP.number_of_edges())

    # blockPrint()

    return GraphFromSHP


# %%

def getShortestPathAStar(source, target, source_line, target_line, source_line_oneway, target_line_oneway, filepath_shp, ignore_oneway=False):

    # make sure global variable is used. 'DG' is a Directed Graph
    global DG

    # blockPrint()

    path = None
    path_length = None
    path_IDs = []
    all_added_edges = []

    def air_line_distance(source, target):
        '''Heuristic function for A star algorithm: returns the air line distance from the source to the target.

        Parameters
        ----------
        source : tuple (float, float)
            This is the source point. It is a tuple of the form (lng1, lat1) where lng1 is the longitude of the source and lat1 is the latitude of the source.
        target : tuple (float, float)
            This is the source point. It is a tuple of the form (lng2, lat2) where lng2 is the longitude of the target and lat2 is the latitude of the target.

        Notes
        -----
        This function is used as the heuristic function for the A Star algorithm.
        It calculates the air line distance between two gps points.
        The result is not 100% correct as the function does not consider the elipsis-like shape of the earth. Instead it just uses an earth radius of 6371000 meters for the calculation.

        Examples
        --------

        >>> myAirLineDist = distFrom((-122.115, 37.115), (-122.111, 37.111))
        >>> myAirLineDist
        568.8872918546489


        Returns
        -------
        int
            The air line distance between two points (source and target) in meters.
        '''
        distance = distFrom(source[0], source[1], target[0], target[1])
        return distance

    def temporarily_add_edge_to_graph(startNode, endNode, edgeWeight, edgeId, direction):
        '''Temporarily adds an edge to the global graph.

        Parameters
        ----------
        startNode : tuple (float, float)
            [description]
        endNode : tuple (float, float)
            [description]
        edgeWeight : int
            The edgeweight is the distance between the startNode and the endNode in meters.
        edgeId : int
            The id of the street (LineString) which the linesegment belongs to according to the shapefile.
        direction : {'B', 'F', 'T'}
            The 'oneway'-property of the street (LineString) which the linesegment belongs to according to the shapefile.
            The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).

        Notes
        -----
        This function temporarily adds an edge to the global graph. It only adds the edge if it deosn't exist in the graph yet. Further, the new edge is added to the list 'all_added_edges' so that it can be removed again in the end.
        '''
        if (not DG.has_edge(startNode, endNode)):
            # print('edge has to be added')
            # print('adding ({},{})'.format(startNode, endNode))
            # print('')

            DG.add_edge(startNode, endNode, weight=edgeWeight,
                        id=edgeId, oneway=direction)
            all_added_edges.append((startNode, endNode))
        # else:
            # print('graph does have({},{})'.format(startNode,endNode))
            # print('edge:',DG.get_edge_data(startNode, endNode))
            # print('')

    # check if the global variable 'DG' is an empty directed graph. If yes, create a graph from the shp-file content.
    if(nx.is_empty(DG)):
        DG = (createGraphFromSHPInput(filepath_shp))

    # check if target lies exactly on the beginning/end of a line segment
    # if yes, target already exists as a node
    if (target not in target_line):
        # remove the target_line edges and add the line_segments as edges

        d_target = LineString(target_line).project(Point(target))

        cut_line_target = cut(LineString(target_line), d_target, Point(target))
        # print('cut_line_target: ', cut_line_target)
        new_line_target_after_cut = [list(x.coords) for x in cut_line_target]

        len_line_segment = len(new_line_target_after_cut[0])
        edge_to_remove_T_source = new_line_target_after_cut[0][len_line_segment-2]
        edge_to_remove_T_target = new_line_target_after_cut[1][1]

        # get the edge (including the attributes) that should be removed from the graph

        if(target_line_oneway == 'B'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_source, edge_to_remove_T_target)
            print('B:', edge_to_remove_from_graph)

        elif(target_line_oneway == 'F'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_source, edge_to_remove_T_target)
            print('F:', edge_to_remove_from_graph)

        elif(target_line_oneway == 'T'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_target, edge_to_remove_T_source)
            print('T:', edge_to_remove_from_graph)

        edge_to_remove_from_graph_id = edge_to_remove_from_graph['id']
        edge_to_remove_from_graph_oneway = edge_to_remove_from_graph['oneway']

        # add the line segments for the target

        if(edge_to_remove_from_graph_oneway == 'B' or ignore_oneway == True):

            # add edges in direction of from-node to to-node

            temporarily_add_edge_to_graph(edge_to_remove_T_source, target, distFrom(
                edge_to_remove_T_source[0], edge_to_remove_T_source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(target, edge_to_remove_T_target, distFrom(
                edge_to_remove_T_target[0], edge_to_remove_T_target[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(edge_to_remove_T_target, target, distFrom(
                edge_to_remove_T_target[0], edge_to_remove_T_target[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(target, edge_to_remove_T_source, distFrom(
                edge_to_remove_T_source[0], edge_to_remove_T_source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

        elif(edge_to_remove_from_graph_oneway == 'F'):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_to_remove_T_source, target, distFrom(
                edge_to_remove_T_source[0], edge_to_remove_T_source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(target, edge_to_remove_T_target, distFrom(
                edge_to_remove_T_target[0], edge_to_remove_T_target[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

        else:
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(target, edge_to_remove_T_source, distFrom(
                edge_to_remove_T_source[0], edge_to_remove_T_source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(edge_to_remove_T_target, target, distFrom(
                edge_to_remove_T_target[0], edge_to_remove_T_target[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

    else:
        print('no need to add additional edges target')

    # check if source lies exactly on the beginning/end of a line segment
    # if yes, source already exists as a node
    if (source not in source_line):
        # remove the source_line edges and add the line_segments as edges

        # if source and target line are equal, the already cut line has to be used

        d_source = LineString(source_line).project(Point(source))
        cut_line_source = cut(LineString(source_line), d_source, Point(source))
        new_line_source_after_cut = [list(x.coords) for x in cut_line_source]

        len_line_segment = len(new_line_source_after_cut[0])
        edge_to_remove_S_source = new_line_source_after_cut[0][len_line_segment-2]
        edge_to_remove_S_target = new_line_source_after_cut[1][1]

        # get the edge (including the attributes) that should be removed from the graph

        if(source_line_oneway == 'B'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_S_source, edge_to_remove_S_target)

        elif(source_line_oneway == 'F'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_S_source, edge_to_remove_S_target)

        elif(source_line_oneway == 'T'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_S_target, edge_to_remove_S_source)

        edge_to_remove_from_graph_id = edge_to_remove_from_graph['id']
        edge_to_remove_from_graph_oneway = edge_to_remove_from_graph['oneway']

        # if source and target line on the same line segment, add a connection between the two
        if((source_line == target_line) and (target not in target_line) and (edge_to_remove_T_source == edge_to_remove_S_source)and (edge_to_remove_T_target == edge_to_remove_S_target)):

            # check which point comes first on the line to then know which direction the line between the two should be
            d_target_same_line = LineString(source_line).project(Point(target))
            d_source_same_line = LineString(source_line).project(Point(source))

            # get the edge (including the attributes) that should be removed from the graph

            if(source_line_oneway == 'B' or ignore_oneway == True):
                temporarily_add_edge_to_graph(source, target, distFrom(
                    source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
                temporarily_add_edge_to_graph(target, source, distFrom(
                    source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

            elif(source_line_oneway == 'F'):
                # if the source_line (which equals target_line) is onway from the from-node to the to-node, the point which is closer at the beginning of the line has to be the starting point of the directed edge which is added to the graph
                if(d_source_same_line > d_target_same_line):
                    temporarily_add_edge_to_graph(target, source, distFrom(
                        source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
                else:
                    temporarily_add_edge_to_graph(source, target, distFrom(
                        source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

            elif(source_line_oneway == 'T'):
                # if the source_line (which equals target_line) is oneway from the to-node to the from-node, the point which is closer at the end of the line has to be the starting point of the directed edge which is added to the graph
                if(d_source_same_line > d_target_same_line):
                    temporarily_add_edge_to_graph(source, target, distFrom(
                        source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
                else:
                    temporarily_add_edge_to_graph(target, source, distFrom(
                        source[0], source[1], target[0], target[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

        # add the line segments for the source
        if(edge_to_remove_from_graph_oneway == 'B' or ignore_oneway == True):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_to_remove_S_source, source, distFrom(
                edge_to_remove_S_source[0], edge_to_remove_S_source[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(source, edge_to_remove_S_target, distFrom(
                edge_to_remove_S_target[0], edge_to_remove_S_target[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(edge_to_remove_S_target, source, distFrom(
                edge_to_remove_S_target[0], edge_to_remove_S_target[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(source, edge_to_remove_S_source, distFrom(
                edge_to_remove_S_source[0], edge_to_remove_S_source[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
        elif(edge_to_remove_from_graph_oneway == 'F'):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_to_remove_S_source, source, distFrom(
                edge_to_remove_S_source[0], edge_to_remove_S_source[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(source, edge_to_remove_S_target, distFrom(
                edge_to_remove_S_target[0], edge_to_remove_S_target[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
        else:
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(edge_to_remove_S_target, source, distFrom(
                edge_to_remove_S_target[0], edge_to_remove_S_target[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)
            temporarily_add_edge_to_graph(source, edge_to_remove_S_source, distFrom(
                edge_to_remove_S_source[0], edge_to_remove_S_source[1], source[0], source[1]), edge_to_remove_from_graph_id, edge_to_remove_from_graph_oneway)

    try:

        # try to find a path
        path = nx.astar_path(DG, source, target,
                             heuristic=air_line_distance, weight='weight')

        path_length = nx.astar_path_length(
            DG, source, target, heuristic=air_line_distance, weight='weight')

        # get the IDs of the traversed lines
        prevNode = (0, 0)
        for node in path:
            if (prevNode != (0, 0)):

                if(DG.get_edge_data(prevNode, node) != None):
                    # get id of edge
                    path_IDs.append(DG.get_edge_data(prevNode, node)['id'])
                else:
                    # edge was not found
                    path_IDs.append('edge not found! -- prevNode:', prevNode)
                    path_IDs.append('edge not found! -- node:', node)
            prevNode = node

    except:
        # no path was found
        pass

    # remove all the edges which where added to the graph.
    DG.remove_edges_from(all_added_edges)

    return path, path_length, path_IDs


# %%


def calculateClosestPointAndShortestPath(filepath, filepath_shp, minNumberOfLines=2, aStar=1):
    '''[summary]

    Parameters
    ----------
    filepath : str
        The path where the txt file (which contains the gps trajectories) is stored.
    filepath_shp : str
        The path where the shp file (which contains the street data) is stored.
    minNumberOfLines : int, optional
        [description], by default 2
    aStar : int, optional
        [description], by default 1

    Notes
    -----


    - outlier : Number of data points which where flagged as outliers.
    - taxi_did_not_move : x
    - no_path_found : 
    - cannot_compute_shortest_path_as_previous_point_is_0_0 : 
    - path_from_target_to_source : 
    - checked_other_solution_index : 
    - chose_other_solution_index : 
    - solution_already_lies_on_shortest_path : 
    - no_solution_lies_on_shortest_path : 
    - other_solution_is_worse : 
    - other_solution_no_path_found : 



    Returns
    -------
    calculatedSolution : list


    statistics : list
        Statistics for the txt file (containing gps trajectories) for which the solution is calculated.
        The statistics contain the following numbers (all of type int):

        - outlier : Number of data points which where flagged as outliers.
        - taxi_did_not_move : Number of times when a data point's gps position was identical to the gps position of the previous data point.
        - no_path_found : Number of times no path was found with the A Star algorithm. 
        - cannot_compute_shortest_path_as_previous_point_is_0_0 : 
        - path_from_target_to_source : 
        - checked_other_solution_index : 
        - chose_other_solution_index : 
        - solution_already_lies_on_shortest_path : 
        - no_solution_lies_on_shortest_path : 
        - other_solution_is_worse : 
        - other_solution_no_path_found : 

    '''

    # blockPrint()

    counter = 0

    calculatedSolution = []

    # initialize the dict in which all the statistics of the currently calculated solution will be stored.
    statistics = {}
    statistics['outlier'] = 0
    statistics['taxi_did_not_move'] = 0
    statistics['no_path_found'] = 0
    statistics['cannot_compute_shortest_path_as_previous_point_is_0_0'] = 0
    statistics['path_from_target_to_source'] = 0
    statistics['checked_other_solution_index'] = 0
    statistics['chose_other_solution_index'] = 0
    statistics['solution_already_lies_on_shortest_path'] = 0
    statistics['no_solution_lies_on_shortest_path'] = 0
    statistics['other_solution_is_worse'] = 0
    statistics['other_solution_no_path_found'] = 0

    previous_target = (0, 0)
    previous_source = (0, 0)
    previous_intersected_line = None
    previous_timestamp = None
    previous_intersected_line_oneway = None

    def getLocationResult(filepath_shp, x, y, passenger, timestamp, previous_point, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway, minNumberOfLines=2):

        print('')
        print('')
        print('CHECKING ---------------------------')
        print('')
        print(x)
        print(y)
        print(passenger)
        print(timestamp)
        print(previous_point)
        print(previous_intersected_line)
        print(previous_timestamp)
        print(previous_intersected_line_oneway)

        location_result = {}
        # all the paths who were not changed from source-->target to target-->source have the following property
        path_from_target_to_source = 0

        closest_intersection_x = 0
        closest_intersection_y = 0
        intersected_line = None
        intersected_line_oneway = None

        # set those key which are not necessarily set (e.g. when previous point was an outlier)

        location_result['solution_index'] = 0
        location_result['target'] = previous_point
        if (previous_intersected_line != None):
            location_result['previous_intersected_line'] = previous_intersected_line
        else:
            location_result['previous_intersected_line'] = ''

        if(previous_intersected_line_oneway != None):
            location_result['previous_intersected_line_oneway'] = previous_intersected_line_oneway
        else:
            location_result['previous_intersected_line_oneway'] = ''

        location_result['path_time'] = ''
        location_result['path'] = ''
        location_result['path_length'] = ''
        location_result['air_line_length'] = ''
        location_result['path_length/air_line_length'] = ''
        location_result['velocity_m_s'] = ''
        location_result['pathIDs'] = ''
        location_result['comment'] = ''
        location_result['path_from_target_to_source'] = ''
        location_result['closest_intersection_x'] = ''
        location_result['closest_intersection_y'] = ''
        location_result['relative_position'] = ''
        location_result['relative_position_normalized'] = ''
        location_result['intersected_line'] = ''
        location_result['intersected_line_oneway'] = ''
        location_result['linestring_adjustment_visualization'] = ''

        location_result['solution_id'] = ''

        comment = ''

        location_result['x'] = x
        location_result['y'] = y
        location_result['source'] = (x, y)
        location_result['passenger'] = passenger
        location_result['timestamp'] = timestamp
        location_result['previous_timestamp'] = previous_timestamp

        # how to get the closest line from one point

        focal_pt = Point(x, y)  # radial sweep centre point

        # load the input lines and combine them into one geometry

        # the whole san francisco shp file (which was exported as a copy)
        # TODO: write better comment like:  min nr just needed to increse areasize
        with fiona.open(filepath_shp) as input_lines:
            # define size of area to filter map lines
            areasize = 0.0005
            max_areasize = 0.001
            number_of_lines = 0

            # define number of lines which have to be looked at after filter
            while ((number_of_lines < minNumberOfLines) and (areasize < max_areasize)):

                input_shapes = list(input_lines.items(
                    bbox=((x-areasize), (y-areasize), (x+areasize), (y+areasize))))

                number_of_lines = len(input_shapes)

                # double are size so that if no lines are found, a bigger area is taken into account when filtering
                areasize = areasize*1.25

                # if the areasize already arrived at its maximuum size decrease the set minimum number of lines by one.
                # This makes sure that (if minNumberOfLines>1) a point is only market as an outlier if there is not even one line within the aresize.
                # Without that it is possible that points are market as outliers just because there are not many lines close to the point,
                # even though the point lies very close to a line.

                # if((not (number_of_lines < minNumberOfLines)) and (not (areasize < max_areasize)) and (minNumberOfLines > 1)):
                #    minNumberOfLines -= 1

        # point is oulier only if no line was found
        if(number_of_lines == 0):
            print('This point is an outlier.')
            comment += 'This point is an outlier.'

            # makes sure that next point doesn't calculate shortest path
            previous_point = (0, 0)

            # update statistics
            statistics['outlier'] += 1

        # point is not an outlier
        else:

            solution_dict = {}

            for input_line in input_shapes:

                lineID = input_line[1]['id']

                newLS = LineString(input_line[1]['geometry']['coordinates'])
                # get the distance along the LineString to a point nearest to the point
                relative_position = newLS.project(focal_pt)
                # get the distance normalized to the length of the LineString
                relative_position_normalized = newLS.project(
                    focal_pt, normalized=True)

                closest_point_on_line = newLS.interpolate(relative_position)

                inter_dict_point = {}

                inter_dict_point['closest_point_on_line'] = closest_point_on_line
                inter_dict_point['relative_position_closest_point_on_line'] = relative_position
                inter_dict_point['relative_position_normalized_closest_point_on_line'] = relative_position_normalized
                inter_dict_point['input_line'] = newLS

                inter_dict_point['oneway'] = input_line[1]['properties']['oneway']

                # distance needs to be float so that it can be sorted appropriately afterwards
                inter_dict_point['distance'] = float(
                    focal_pt.distance(closest_point_on_line))

                solution_dict[lineID] = inter_dict_point

            # sort the nested dict by the 'distance' which lies  within all the inner dicts
            # first element in the sorted dict is the point with the shortest distance to the focal_pt
            solution_dict_sorted_by_distance = sorted(
                solution_dict.items(), key=lambda kv: (kv[1]['distance']))

            print('solution NEW:')
            # get the first element as it is the one with the smallest distance (which is the closest point)
            solution_id = solution_dict_sorted_by_distance[0][0]
            solution = solution_dict_sorted_by_distance[0][1]
            location_result['all_solutions_sorted'] = solution_dict_sorted_by_distance
            location_result['solution_id'] = solution_id
            location_result['solution'] = solution

            for s in solution:
                print(str(s))
                print(str(solution[s]))
            print('')
            print('end solution')
            print('')
            print('')
            print('')

            closest_intersection_x = solution['closest_point_on_line'].x
            closest_intersection_y = solution['closest_point_on_line'].y
            relative_position = solution['relative_position_closest_point_on_line']
            relative_position_normalized = solution['relative_position_normalized_closest_point_on_line']

            intersected_line = solution['input_line']
            intersected_line_oneway = solution['oneway']

            # when appending the solution we have to parse the floats into strings
            location_result['closest_intersection_x'] = closest_intersection_x
            location_result['closest_intersection_y'] = closest_intersection_y

            # append the relative position / and the normalized relative position of the solution on the intersected line
            location_result['relative_position'] = relative_position
            location_result['relative_position_normalized'] = relative_position_normalized

            # append the underlying map line on which the final point lies AND the oneway_property
            location_result['intersected_line'] = intersected_line
            location_result['intersected_line_oneway'] = intersected_line_oneway

            # append the line which visualizes the way from the start point to the new position of the closest intersection point
            linestring_adjustment_visualization = LineString(
                [(focal_pt.x, focal_pt.y), (closest_intersection_x, closest_intersection_y)])

            location_result['linestring_adjustment_visualization'] = linestring_adjustment_visualization

            # check if previous point is set. if yes: calculate the shortest path from current point to previous point
            if(previous_point == (0, 0)):
                # no previous point set, do not calculate shortest path.
                print('cannot compute shortest path as previous point is (0,0)')
                comment += 'Cannot compute shortest path as previous point is (0,0). '

                # update statistics
                statistics['cannot_compute_shortest_path_as_previous_point_is_0_0'] += 1

            elif((closest_intersection_x, closest_intersection_y) == previous_point):
                # if the taxi did not move, path is not created
                print('source==target. Taxi did not move.')

                comment += 'source==target. Taxi did not move. '

                location_result['path_length'] = 0
                location_result['air_line_length'] = 0
                location_result['path_length/air_line_length'] = 1
                location_result['velocity_m_s'] = 0

                # update statistics
                statistics['taxi_did_not_move'] += 1

            else:

                path_time = abs(previous_timestamp-timestamp)

                location_result['path_time'] = path_time

                if(aStar == 1):

                    path = None
                    path2 = None
                    # calculate the shortest path with A STAR algorithm
                    path, path_length, pathIDs = getShortestPathAStar((closest_intersection_x, closest_intersection_y), previous_point, list(
                        intersected_line.coords), list(previous_intersected_line.coords), intersected_line_oneway, previous_intersected_line_oneway, filepath_shp)

                    # calculate air line distance between source and target
                    air_line_length = distFrom(
                        closest_intersection_x, closest_intersection_y, previous_point[0], previous_point[1])

                    # if this is true, it is possible that a taxi has to ride all around the block because it is a oneway street and source and target lie on the same line
                    if(air_line_length < 20 and intersected_line_oneway != 'B' and previous_intersected_line_oneway != 'B' and (list(intersected_line.coords) == list(previous_intersected_line.coords))):
                        print('if statement inside')

                        # calculate where the source lies on the line
                        d_source = LineString(list(intersected_line.coords)).project(
                            Point(closest_intersection_x, closest_intersection_y))

                        # calculate where the target lies on the line
                        d_target = LineString(list(previous_intersected_line.coords)).project(
                            Point(previous_point))

                        print('d_source:', d_source)
                        print('d_target:', d_target)

                        # check if there was a glipse in the gps coordinates which resulted in a situation where the source lies after the target
                        # if it is a oneway line from the FROM-node to the TO-node and source lies AFTER the target   OR   if it is a oneway line from the TO-node to the FROM-node and source lies BEFORE the target
                        if((intersected_line_oneway == 'F' and d_source > d_target)or(intersected_line_oneway == 'T' and d_source < d_target)):
                            print('if statement inside SECOND!!!')

                            # calculate the shortest path with A STAR algorithm
                            # set ignore_oneway=True : this makes sure that the source_line and the target_line both are treated as lines where driving in both directions is allowed (so feature 'oneway' is ignored)
                            path2, path_length2, pathIDs2 = getShortestPathAStar((closest_intersection_x, closest_intersection_y), previous_point, list(intersected_line.coords), list(
                                previous_intersected_line.coords), intersected_line_oneway, previous_intersected_line_oneway, filepath_shp, ignore_oneway=True)

                            print('path_length:', path_length)
                            print('path_length2:', path_length2)
                            print('path:', path)
                            print('path2:', path2)
                            print('')
                            print('check if path2 is better:')
                            print(path == None and path2 != None)
                            print('.. or ..')
                            print(path != None and path2 != None and (
                                path_length > path_length2))

                        # override path/path_lenght/pathIDs to make sure to use the solution which ignores the 'oneway' property
                        if((path == None and path2 != None) or (path != None and path2 != None and (path_length > path_length2))):
                            print('OVERRIDING THE PATH')
                            path = path2
                            path_length = path_length2
                            pathIDs = pathIDs2
                            comment += 'The oneway-property was ignored. '
                            path_from_target_to_source = 1

                            # update statistics
                            statistics['path_from_target_to_source'] += 1

                    # if a path was found, append the solution to the calculatedSolution
                    if(path != None):
                        velocity_m_s = path_length/path_time

                        location_result['path'] = LineString(path)
                        location_result['path_length'] = path_length
                        location_result['air_line_length'] = air_line_length
                        location_result['path_length/air_line_length'] = path_length / \
                            air_line_length
                        location_result['velocity_m_s'] = velocity_m_s
                        location_result['pathIDs'] = pathIDs
                        location_result['comment'] = comment

                    else:

                        # no previous point set, do not calculate shortest path.
                        comment += 'No path was found with A Star algorithm. '
                        print('No path was found with A Star algorithm.')

        location_result['path_from_target_to_source'] = path_from_target_to_source
        location_result['comment'] = comment

        return location_result, (closest_intersection_x, closest_intersection_y), previous_point, intersected_line, previous_intersected_line, timestamp, intersected_line_oneway, previous_intersected_line_oneway

    # calculate the number of lines of the file we are looking at
    # lines_in_textfile=sum(1 for line in open(filepath))
    lines_in_textfile = 0
    with open(filepath, 'r') as f:
        lines_in_textfile = sum(1 for line in f)

    print(lines_in_textfile)

    # Initial call to print 0% progress ProgressBar
    suffix = '| current file: {}/{} lines'.format(counter+1, lines_in_textfile)
    printProgressBar(0, lines_in_textfile, prefix='Progress:',
                     suffix=suffix, length=50)

    previous_location_result = {}

    with open(filepath, 'r') as f:
        mylist = f.read().splitlines()

        # append an artificaial line at the end. This will make sure that also the path from the last to the second-last element will be calculated and possible improvements will be checked.
        mylist.append('artificialline')

        for lines in mylist:

            # the last line of the text file is artificial as it was added before. this line is just needed to make sure the second last line (which is actually the last actual line of the text file is calculated correctly)
            if(lines != 'artificialline'):

                values = lines.split(' ')
                # the longitude and latitude were read from the txt file an saved as a string, here we transform it to a float
                x = float(values[1])
                y = float(values[0])
                passenger = int(values[2])
                timestamp = int(values[3])

                # get the result for the current location
                current_location_result, source, target, intersected_line, target_intersected_line, timestamp, intersected_line_oneway, target_intersected_line_oneway = getLocationResult(
                    filepath_shp, x, y, passenger, timestamp, previous_source, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway, minNumberOfLines)

                print('lets check previous_target:', previous_target)
                print('lets check target:', target)
                print('lets check previous_target2:', previous_target)

            if (previous_location_result != {}):

                # check if the result of the previous location can be improved
                # we check if either velocity>35m/s or if path_length/air_line_length > 2

                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')
                print('')

                # TODO: check why it crashes if previous_location_result['path']!='' is removed
                # TODO: why was there no path found???
                if(previous_target != (0, 0) and target != (0, 0) and previous_location_result['path'] != ''):

                    print('')
                    print('')
                    print(
                        'previous_location_result[velocity_m_s]:', previous_location_result['velocity_m_s'])
                    print('previous_location_result[path_length/air_line_length]:',
                          previous_location_result['path_length/air_line_length'])
                    print('CHECK1 OTHER SOLUTION???? because PATH=NONE',
                          (previous_location_result['path'] == ''))
                    print('CHECK2 OTHER SOLUTION???? because PREVIOUS velocity or length to big', ((float(
                        previous_location_result['velocity_m_s']) > 35.0) or ((previous_location_result['path_length/air_line_length']) > 2.0)))
                    print(
                        'current_location_result[velocity_m_s]:', current_location_result['velocity_m_s'])
                    print(current_location_result['velocity_m_s'] == '')
                    print('CHECK3 OTHER SOLUTION???? because CURRENT velocity or length to big', (current_location_result['velocity_m_s'] != '' and (
                        (float(current_location_result['velocity_m_s']) > 35.0) or ((current_location_result['path_length/air_line_length']) > 2.0))))
                    print('')

                    # check if either the path from current source to target or from previous source to previous target might be needed to improved
                    # if (current_location_result['velocity_m_s']='') means that the taxi was not moving at the previous point --> for this reason it should not be checked
                    if(current_location_result['velocity_m_s'] != '' and ((previous_location_result['path'] == '') or (float(previous_location_result['velocity_m_s']) > 35.0) or ((previous_location_result['path_length/air_line_length']) > 2.0) or (float(current_location_result['velocity_m_s']) > 35.0) or ((current_location_result['path_length/air_line_length']) > 2.0))):

                        # update statistics
                        statistics['checked_other_solution_index'] += 1

                        print('OKAY: previous_target is known')
                        print('source:', str(source))
                        print('target:', str(target))
                        print('previous_target', str(previous_target))

                        print('')
                        print('')
                        print('')
                        print(
                            'LETS calculate shortest path between current source and previous target')
                        print('')

                        # check if there is another possible solution
                        # first calculate shortest path between current source and previous target
                        path_to_previous_target, path_length_to_previous_target, pathIDs_to_previous_target = getShortestPathAStar(source, previous_target, list(
                            intersected_line.coords), list(previous_target_intersected_line.coords), intersected_line_oneway, previous_target_intersected_line_oneway, filepath_shp)
                        print('path_to_previous_target:',
                              str(path_to_previous_target))
                        print('')
                        print('path_length_to_previous_target:',
                              path_length_to_previous_target)
                        print('')
                        print('pathIDs_to_previous_target:',
                              pathIDs_to_previous_target)
                        print('')
                        print('')
                        print('all_solutions_sorted:',
                              previous_location_result['all_solutions_sorted'])
                        print('solution_id:',
                              previous_location_result['solution_id'])
                        print('solution:',
                              previous_location_result['solution'])

                        print('')
                        print('FOR LOOP:')

                        if(str(previous_location_result['solution_id']) in pathIDs_to_previous_target or previous_location_result['solution_id'] in pathIDs_to_previous_target):
                            print(
                                'Solution is already on the shortest path between the points before and after --> therefore the solution can not be improved.')

                            # update statistics
                            statistics['solution_already_lies_on_shortest_path'] += 1

                        else:
                            print('trying to find solution on shortest path.')

                            new_solution_id = (-1)
                            new_solution_index = 0
                            new_solution = {}

                            for k, v in previous_location_result['all_solutions_sorted']:
                                if(str(k) in pathIDs_to_previous_target):
                                    new_solution_id = k
                                    new_solution = v
                                    print('now break')
                                    break
                                print('k:', k)
                                print('v:', v)
                                new_solution_index += 1

                            if (new_solution_id == (-1)):
                                print(
                                    'There is no other solution which lies on the shortest path.')

                                # update statistics
                                statistics['no_solution_lies_on_shortest_path'] += 1

                            else:
                                # check if this solution is actually better
                                print('Check if this solution is actually better.')
                                print('found new_solution_id:', new_solution_id)
                                print('found new_solution:', new_solution)
                                previous_commulated_path_lengt = float(
                                    previous_location_result['path_length']) + float(current_location_result['path_length'])
                                print('previous_commulated_path_lengt',
                                      previous_commulated_path_lengt)

                                new_solution_point_x = list(
                                    new_solution['closest_point_on_line'].coords)[0][0]
                                new_solution_point_y = list(
                                    new_solution['closest_point_on_line'].coords)[0][1]
                                print('new_solution_point_x:',
                                      new_solution_point_x)
                                print('new_solution_point_y:',
                                      new_solution_point_y)

                                current_location_result_new, source_new, target_new, intersected_line_new, target_intersected_line_new, timestamp_new, intersected_line_oneway_new, target_intersected_line_oneway_new = getLocationResult(
                                    filepath_shp, x, y, passenger, timestamp, (new_solution_point_x, new_solution_point_y), new_solution['input_line'], previous_location_result['timestamp'], new_solution['oneway'], minNumberOfLines)

                                print('current_location_result_new:',
                                      current_location_result_new)
                                print('current_location_result_new -- PATH:',
                                      str(current_location_result_new['path']))
                                print('source_new:', source_new)
                                print('target_new:', target_new)
                                print('intersected_line_new:',
                                      intersected_line_new)
                                print('target_intersected_line_new:',
                                      target_intersected_line_new)
                                print('timestamp_new:', timestamp_new)
                                print('intersected_line_oneway_new:',
                                      intersected_line_oneway_new)
                                print('target_intersected_line_oneway_new:',
                                      target_intersected_line_oneway_new)

                                previous_location_result_new, previous_source_new, previous_target_new, previous_intersected_line_new, previous_target_intersected_line_new, previous_timestamp_new, previous_intersected_line_oneway_new, previous_target_intersected_line_oneway_new = getLocationResult(
                                    filepath_shp, new_solution_point_x, new_solution_point_y, previous_location_result['passenger'], previous_location_result['timestamp'], previous_location_result['target'], previous_location_result['previous_intersected_line'], previous_location_result['previous_timestamp'], previous_location_result['previous_intersected_line_oneway'], minNumberOfLines)
                                print('previous_location_result_new:',
                                      previous_location_result_new)
                                print('previous_location_result_new -- PATH:',
                                      str(previous_location_result_new['path']))
                                print('previous_source_new:',
                                      previous_source_new)
                                print('previous_target_new:',
                                      previous_target_new)
                                print('previous_intersected_line_new:',
                                      previous_intersected_line_new)
                                print('previous_target_intersected_line_new:',
                                      previous_target_intersected_line_new)
                                print('previous_timestamp_new:',
                                      previous_timestamp_new)
                                print('previous_intersected_line_oneway_new:',
                                      previous_intersected_line_oneway_new)
                                print('previous_target_intersected_line_oneway_new:',
                                      previous_target_intersected_line_oneway_new)

                                if(previous_location_result_new['path_length'] != '' and current_location_result_new['path_length'] != ''):

                                    new_commulated_path_lengt = float(
                                        previous_location_result_new['path_length']) + float(current_location_result_new['path_length'])
                                    print('')
                                    print('')
                                    print('')
                                    print('previous_commulated_path_lengt:',
                                          previous_commulated_path_lengt)
                                    print('new_commulated_path_lengt:',
                                          new_commulated_path_lengt)
                                    print('is the new solution better???', (
                                        previous_commulated_path_lengt > new_commulated_path_lengt))
                                    print('')

                                    # if the new cummulated path is shorter set it as the correct one
                                    if((previous_commulated_path_lengt > new_commulated_path_lengt)):
                                        print(
                                            'previous_location_result[solution_id]', previous_location_result['solution_id'])
                                        print('new_solution_id:',
                                              new_solution_id)

                                        print_str = 'Set the second checked solution (with solution_id: {}) as the correct one since the new cummulated path was shorter than the initial solution (with solution_id: {}). '.format(
                                            new_solution_id, previous_location_result['solution_id'])
                                        print(print_str)
                                        new_comment = 'The source point of this path was reset as this leads to much shorter path lengths. Previous solution: {} / New solution: {}. '.format(
                                            previous_source, (new_solution_point_x, new_solution_point_y))

                                        previous_location_result = previous_location_result_new
                                        previous_location_result['comment'] += new_comment
                                        previous_location_result['solution_id'] = new_solution_id
                                        previous_location_result['solution_index'] = new_solution_index

                                        current_location_result = current_location_result_new
                                        current_location_result['comment'] += 'The target point of this path was reset as this leads to much shorter path lengths.'

                                        # update statistics
                                        statistics['chose_other_solution_index'] += 1

                                    else:
                                        new_comment = 'Checked the second solution (solution_id: {} / solution_index: {} / new cummulated path_length: {}) initial solution (solution_id: {} / solution_index: {} / initial cummulated path_length: {}) was chosen due to a shorter cummulated path_length.'.format(
                                            new_solution_id, new_solution_index, new_commulated_path_lengt, previous_location_result['solution_id'], previous_location_result['solution_index'], previous_commulated_path_lengt)
                                        print(new_comment)
                                        previous_location_result['comment'] += new_comment

                                        # update statistics
                                        statistics['other_solution_is_worse'] += 1

                                else:
                                    # no path was found with new solution
                                    new_comment = 'No path found for new solution. '
                                    print(new_comment)
                                    previous_location_result['comment'] += new_comment

                                    # update statistics
                                    statistics['other_solution_no_path_found'] += 1

                else:
                    print('Tried to check other solution but the target is not known. ')

                # append the result of the previous location to the entire solution
                calculatedSolution.append(previous_location_result)

                # check if a path was found xxx
                if (previous_location_result['path'] == ''):
                    # update statistics
                    statistics['no_path_found'] += 1

            # set the result of the current location as the result of the previous location
            previous_location_result = current_location_result

            previous_source = source
            previous_target = target

            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('')
            print('previous_source', previous_source)
            print('previous_target', previous_target)
            print('')
            print('')
            print('')
            print('')

            previous_intersected_line = intersected_line
            previous_target_intersected_line = target_intersected_line
            previous_timestamp = timestamp
            previous_intersected_line_oneway = intersected_line_oneway
            previous_target_intersected_line_oneway = target_intersected_line_oneway

            # TODO: Delete the following lines as this is just for testing.
            # stop script after a few lines (for testing reasons)
            counter += 1
            print(counter)
            if (counter > 10):
                print('ende')
                #break

            # Update Progress Bar

            # enablePrint()

            suffix = '| current file: {}/{} lines'.format(
                counter, lines_in_textfile)
            x = 0
            if (x):
                suffix += ' | total: {} of {} files'.format(
                    counter, lines_in_textfile)
            printProgressBar(counter, lines_in_textfile,
                             prefix='Progress:', suffix=suffix, length=50)

            # blockPrint()

        '''
        NOT NEEDED ANYMORE, as we have artificial line now

        # the result of a location is only added to the intire solution after it is confirmed as the best solution by the next location
        # for this reason the result of the last location of the text_file has to be added to the entire solution (since there isn't a next solution that will do that)
        calculatedSolution.append(previous_location_result)
        '''

    return calculatedSolution, statistics


def getTimeDifferences(filepath, timestampPosition):

    # initialise empty list
    timeDifferences = []
    # set previousTimestamp to 0
    previousTimestamp = 0

    # open the file
    with open(filepath, 'r') as f:
        mylist = f.read().splitlines()

        # loop through the txt file
        for lines in mylist:
            values = lines.split(' ')
            # get the timestamp of the current line
            timestamp = int(values[timestampPosition])

            # if it't not the forst entry of the txt file, append the difference to the set
            if(previousTimestamp != 0):
                timeDifferences.append(previousTimestamp-timestamp)

            previousTimestamp = timestamp

    return timeDifferences


def plotAndSaveHistogram(input, minXLabel, maxXLabel, binSize, filename, title, xlabel):
    # get the max value of the imputs
    maxInput = max(input)

    # check if inputs are floats (as in the case of velocities)
    if isinstance(maxInput, float):
        # round the float to the next higher number and then convert it to int
        maxInput = int(math.ceil(maxInput))

    # create the bins (as a range from 0 to maxXLabel or maxInput depending on which one is smaller). +2 is needed to display also the max number.
    if isinstance(binSize, float):
        # create a range with floats
        myBins = np.arange(minXLabel, min(
            maxXLabel, maxInput)+(min(2*binSize, 2)), binSize)
    else:
        # create a range with ints
        myBins = range(minXLabel, min(maxXLabel, maxInput) +
                       (min(2*binSize, 2)), binSize)

    # initialize the figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # plot the histogram. 'np.clip' makes sure that the time differences which are bigger then the max X value are visualized in the last bin.
    n, bins, patches = ax.hist(x=np.clip(
        input, 0, myBins[-1]), bins=myBins, color='#0504aa', alpha=0.7, rwidth=0.85)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(xlabel)
    plt.ylabel('Frequency')
    plt.title(title)
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

    if isinstance(binSize, float):
        xlabels = bins[0:].astype('|S3')
    else:
        xlabels = bins[0:].astype(str)

    # remove the last x label
    xlabels = xlabels[:-1]
    # add a '+' to the label of the last bin to visualize that this bin contains all timeDifferences  which are equal or bigger than the xTick.
    if not isinstance(binSize, float):
        xlabels[-1] += '+'

    myXTicks = np.array(myBins)
    # set the x ticks for the figure
    plt.xticks(myXTicks, xlabels)

    # save the figure in the folder of the current taxi
    plt.savefig(filename)
    # plt.show()
    # close the current figure
    plt.close()

# %%


def getFilename(path):
    '''Returns the name of a file from a specific path (excluding directories and extension).

    Parameters
    ----------
    path : [string]
        The path of a specific data file.

    Returns
    -------
    string
        The name of the file. Namely, the last part of the path (excluding the extgension).


    Examples
    --------

    >>> myPath = 'dir1/dir2/MyFilename.txt'
    >>> filename = getFilename(myPath)
    >>> filename
    'MyFilename'

    '''
    # get the name of the file from the path (without the directories and without extension)
    filename = os.path.splitext(os.path.basename(os.path.normpath(path)))[0]
    return filename


# %%
# only selected data points from this taxi
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/new_abboip_selection_test.txt'
# all data points from this taxi
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_selection_for_SP_tests.txt'

def main():

    # blockPrint()

    def caculationForOneTXTFile(filepath_shp, dirName, filepath, aStar=1):

        # set the header of the output txt file which will contain the calculated solution.
        header = ['latitude(y);longitude(x);hasPassenger;time;closest_intersection_x;closest_intersection_y;relative_position;relative_position_normalized;intersected_line_oneway;intersected_line_as_linestring;linestring_adjustment_visualization;path_time']

        if(aStar == 1):
            header += ';path_as_linestring;path_length;air_line_length;path_length/air_line_length;velocity_m_s;pathIDs;solution_id;solution_index;path_from_target_to_source;comment'

        header += '\n'

        new_filename_solution = os.path.join(
            dirName, 'calculatedSolution') + '.txt'
        new_filename_statistics = os.path.join(dirName, 'statistics') + '.txt'
        new_filename_velocities = os.path.join(dirName, 'velocitiesPLOT.png')
        new_filename_path_length_air_line_length = os.path.join(
            dirName, 'path_length_air_line_length_PLOT.png')

        myCalculatedSolution, mySolutionStatistics = calculateClosestPointAndShortestPath(
            filepath, filepath_shp, minNumberOfLines=2, aStar=aStar)

        # this saves a new text file which includes the calculated parameters
        with open(new_filename_solution, 'w') as new_file:
            # write the header to the new txt file
            new_file.writelines(header)

            # write the calculated solution to the new txt file
            for location_result in myCalculatedSolution:
                new_file.write(str(location_result['y']))
                new_file.write(';')
                new_file.write(str(location_result['x']))
                new_file.write(';')
                new_file.write(str(location_result['passenger']))
                new_file.write(';')
                new_file.write(str(location_result['timestamp']))
                new_file.write(';')
                new_file.write(str(location_result['closest_intersection_x']))
                new_file.write(';')
                new_file.write(str(location_result['closest_intersection_y']))
                new_file.write(';')
                new_file.write(str(location_result['relative_position']))
                new_file.write(';')
                new_file.write(
                    str(location_result['relative_position_normalized']))
                new_file.write(';')
                new_file.write(str(location_result['intersected_line_oneway']))
                new_file.write(';')
                new_file.write(str(location_result['intersected_line']))
                new_file.write(';')
                new_file.write(
                    str(location_result['linestring_adjustment_visualization']))
                new_file.write(';')
                new_file.write(str(location_result['path_time']))
                new_file.write(';')
                new_file.write(str(location_result['path']))
                new_file.write(';')
                new_file.write(str(location_result['path_length']))
                new_file.write(';')
                new_file.write(str(location_result['air_line_length']))
                new_file.write(';')
                new_file.write(
                    str(location_result['path_length/air_line_length']))
                new_file.write(';')
                new_file.write(str(location_result['velocity_m_s']))
                new_file.write(';')
                new_file.write(str(location_result['pathIDs']))
                new_file.write(';')
                new_file.write(str(location_result['solution_id']))
                new_file.write(';')
                new_file.write(str(location_result['solution_index']))
                new_file.write(';')
                new_file.write(
                    str(location_result['path_from_target_to_source']))
                new_file.write(';')
                new_file.write(str(location_result['comment']))
                new_file.write('\n')

        # initialise empty lists
        velocities = []
        path_length_air_line_length = []
        velocities_none_counter = 0

        # get all velocities of the solution
        for location_result in myCalculatedSolution:
            if (location_result['velocity_m_s'] == ''):
                velocities_none_counter += 1
            else:
                velocities.append(float(location_result['velocity_m_s']))

            if (location_result['path_length/air_line_length'] != ''):
                path_length_air_line_length.append(
                    float(location_result['path_length/air_line_length']))

        # plot the timeDifferences in a histogram
        plotAndSaveHistogram(velocities, 0, 80, 5, new_filename_velocities,
                             'Histogram of velocities between gps points', 'velocity in m/s')

        # plot the the path_length/air_line_length in a histogram
        plotAndSaveHistogram(path_length_air_line_length, 1.0, 2.5, 0.1, new_filename_path_length_air_line_length,
                             'Histogram of path lengths divided by air line lengths', 'path length / air line length')

        # SAVE STATISTICS IN NEW FILE

        with open(new_filename_statistics, 'w') as new_file:
            # write statistics to the file
            new_file.write('path of shp file: ')
            new_file.write(str(filepath_shp))
            new_file.write('\n')
            new_file.write('path of old file: ')
            new_file.write(str(filepath))
            new_file.write('\n')
            new_file.write('path of new file: ')
            new_file.write(str(new_filename_solution))
            new_file.write('\n')

            # calculate the number of lines in the txt file
            with open(filepath, 'r') as f:
                num_lines = sum(1 for line in f)
            new_file.write('number of lines in txt file: ')
            new_file.write(str(num_lines))

            new_file.write('\n')
            new_file.write('\n')

            # write the statistics whcih were returned from the algorithm to the file
            for item in mySolutionStatistics.items():
                new_file.write(str(item[0]))
                new_file.write(': ')
                new_file.write(str(item[1]))
                new_file.write('\n')

            new_file.write('\n')
            new_file.write('VELOCITIES PLOT: ')
            new_file.write('\n')
            new_file.write('The Velocity plot contains velocities from {} data points. For the remaining {} data points the velocity could not be calculated (e.g. because it is an outlier).'.format(
                num_lines-velocities_none_counter, velocities_none_counter))
            new_file.write('\n')
            new_file.write('\n')

        # enablePrint()

    #filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy.txt'

    # filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small.txt'
    # filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall.txt'
    # filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_2.txt'
    # filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_3.txt'

    # filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_4.txt'

    #filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_verysmall_closestIsBest_oneline.txt'

    #filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_SMALL_to_check_outlier.txt'

    ##filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_verysmall_closestIsBest_1stSolution.txt'

    #filepath = '/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy_check_wrong_outliers.txt'

    filepath1 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy.txt'
    filepath2 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abgibo_copy.txt'
    filepath3 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abmuyawm_copy.txt'
    filepath4 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abniar_copy.txt'
    filepath5 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abnovkak_copy.txt'
    filepath6 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abtyff_copy.txt'
    filepath7 = '/Users/Joechi/Google Drive/gps2net/Data/test_data/other_taxi/new_abwecnij_copy.txt'

    #filepath_shp = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/geo_SF_lines_exported_for_testing.shp'

    filepath_shp = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/San Francisco Basemap Street Centerlines/geo_export_e5dd0539-2344-4e87-b198-d50274be8e1d.shp'

    filepaths = []
    filepaths.append(filepath6)
    filepaths.append(filepath7)
    filepaths.append(filepath1)
    filepaths.append(filepath2)
    filepaths.append(filepath3)
    filepaths.append(filepath4)
    filepaths.append(filepath5)

    # loop through all the filepaths
    for path in filepaths:

        new_filename = getFilename(path)

        dirName = os.path.join('output_files', new_filename)

        #dirName += new_filename

        # os.makedirs(dirName)

        # Create target directory & all intermediate directories if don't exists
        try:
            os.makedirs(dirName)
            print('Directory ', dirName,  ' created.')
        except FileExistsError:
            print('Directory ', dirName,  ' already exists.')

        caculationForOneTXTFile(filepath_shp, dirName, path, 1)

        # get all timestamp differences of a text file
        timeDifferences = getTimeDifferences(path, 3)

        timedifferencesFileName = os.path.join(
            dirName, 'timedifferencesPLOT.png')

        # plot the timeDifferences in a histogram
        plotAndSaveHistogram(timeDifferences, 0, 300, 25, timedifferencesFileName,
                             'Histogram of time differences between gps points', 'time difference in seconds')


# %%
if __name__ == '__main__':
    import doctest
    # running 'doctest.testmod()' test the examples in docstrings. Alternatively, the examples in docstrings can be tested by navigating to the 'docs' directory and running the following command in the terminal: python gps2net.py -v
    doctest.testmod()

    # run the main method
    main()
