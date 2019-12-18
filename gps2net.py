# %%
import doctest
import math
import os
import sys

import fiona
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from shapely.geometry import LineString, Point


# global variable: empty Directed Graph
DG = nx.DiGraph()
current_txt_file = 0
number_of_txt_files = 0


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
        prefix str, by default ''
    suffix : str, optional
        suffix str, by default ''
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

    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration >= total:
        print()


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

    >>> myDist = distFrom(-122.115, 37.115, -122.111, 37.111)
    >>> myDist
    568.8872918546489
    >>> distFrom(-122.115, 37.115, -122.111, 37.111)
    568.8872918546489


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
    'line.interpolate(distance)' should be equal to the parameter 'point'. However, due to precision errors 'point' is used.

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
    This graph is only initialized once in the beginning of the script. The DiGraph contains all linesegment of the shp file (i.e. streetsegments) as directed edges. The start and end node of each edge are tuples containing the coordinates of the position. The weight of an edge is the air_line_distance from the start to the end of the linesegment. Each edge containes the following attributes:

    - id : int
        The id of the street (LineString) which the linesegment belongs to according to the shapefile.
    - oneway : {'B', 'F', 'T'}
        The 'oneway'-property of the street (LineString) which the linesegment belongs to according to the shapefile. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).

    Depending on the 'oneway'-property of the street either one or two edges are added to the graph. For 'B', two directed edges (from-node --> to-node AND to-node --> from-node) are added. For 'F' or 'T', only one directed edge is added to the graph.

    Returns
    -------
    Directed Graph
        The Directed Graph which contains all linesegments of the shapefiles as edges is returned.
    '''

    GraphFromSHP = nx.DiGraph()

    counter = 0

    # calculate how many street segments the shp file contains
    nr_elements_in_SHP_file = 0
    with fiona.open(filepath_shp) as street_lines:
        nr_elements_in_SHP_file = sum(1 for street_segment in street_lines)

    # Initial call to print 0% progress ProgressBar
    suffix = '| SHP file street segments: {}/{}'.format(
        counter+1, nr_elements_in_SHP_file)
    printProgressBar(0, nr_elements_in_SHP_file,
                     prefix='The Graph is being created:', suffix=suffix, length=50)

    with fiona.open(filepath_shp) as street_lines:

        # loop through all streets
        for street_segment in list(street_lines):

            previous_segment = (0, 0)
            # get the attributes from the street
            id = street_segment['id']
            oneway = street_segment['properties']['oneway']

            # loop through all segments of a street
            for segment in street_segment['geometry']['coordinates']:

                if (previous_segment != (0, 0)):

                    # calculate the length of a street segment
                    calculated_length = distFrom(
                        segment[0], segment[1], previous_segment[0], previous_segment[1])

                    # add streetsegments as edges (depending on the oneway property)
                    if(oneway == 'B'):
                        # create two directed edges for both directions
                        GraphFromSHP.add_edge(
                            segment, previous_segment, weight=calculated_length, id=id, oneway=oneway)
                        GraphFromSHP.add_edge(
                            previous_segment, segment, weight=calculated_length, id=id, oneway=oneway)
                    elif(oneway == 'F'):
                        # create a directed edge from the from-node to the to-node
                        GraphFromSHP.add_edge(
                            previous_segment, segment, weight=calculated_length, id=id, oneway=oneway)
                    elif(oneway == 'T'):
                        # create a directed edge from the to-node to the from-node
                        GraphFromSHP.add_edge(
                            segment, previous_segment, weight=calculated_length, id=id, oneway=oneway)

                previous_segment = segment

            # Update Progress Bar
            suffix = '| SHP file street segments: {}/{}'.format(
                counter+1, nr_elements_in_SHP_file)
            printProgressBar(counter + 1, nr_elements_in_SHP_file,
                             prefix='The Graph is being created:', suffix=suffix, length=50)

            counter += 1

    print('The graph was created and contains {} edges.'.format(
        GraphFromSHP.number_of_edges()))

    return GraphFromSHP


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
    This function is only called in the following function: :func:`~gps2net.getShortestPathAStar`

    This function is used as the heuristic function for the A Star algorithm.
    It calculates the air line distance between two gps points.
    The result is not 100% correct as the function does not consider the elipsis-like shape of the earth. Instead it just uses an earth radius of 6371000 meters for the calculation.

    Examples
    --------

    >>> myAirLineDist = air_line_distance((-122.115, 37.115), (-122.111, 37.111))
    >>> myAirLineDist
    568.8872918546489


    Returns
    -------
    int
        The air line distance between two points (source and target) in meters.
    '''
    distance = distFrom(source[0], source[1], target[0], target[1])
    return distance


def getShortestPathAStar(source, target, source_line, target_line, source_line_oneway, target_line_oneway, filepath_shp, ignore_oneway=False):
    """Calculates the shortest – most likely – path between two GPS positions (on a street segment) based on the streets of the underlying network.

    Parameters
    ----------
    source : tuple of float
        GPS position (float, float). The start node of the path.
    target : tuple of float
        GPS position (float, float). The end node of the path.
    source_line : list of coordinates
        The source_line is a list of GPS coordinates. These coordinates represent the street on which the source point lies.
    target_line : list of coordinates
        The target_line is a list of GPS coordinates. These coordinates represent the street on which the target point lies.
    source_line_oneway : {'B', 'F', 'T'}
        The 'oneway'-property of the source_line.
        The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).
    target_line_oneway : {'B', 'F', 'T'}
        The 'oneway'-property of the target_line.
        The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).
    filepath_shp : str
        The path where the shp file (which contains the street data) is stored.
    ignore_oneway : bool, optional
        By default False. If True the oneway-property will be ignored when adding edges to the graph. This is set to True when the algorithm detects a GPS glipse (the vehicle seems to drive a tiny bit backwards on a oneway street which is not possible) which resulted in a wrong path. In this case the path might actually be the path from the target to the source.

    Returns
    -------
    path : list of coordinates
        The shortest – most likely – path between two data point based on the streets of the underlying network.
    path_length : float
        The length of the path in meters. This is an approximation – see :func:`~gps2net.distFrom`.
    path_IDs : list
        The path IDs of all street segments which are traversed on the path.
    """

    # make sure global variable is used. 'DG' is a Directed Graph
    global DG

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
        This function is only called in the following function: :func:`~gps2net.getShortestPathAStar`

        This function temporarily adds an edge to the global graph. It only adds the edge if it deosn't exist in the graph yet. Further, the new edge is added to the list 'all_added_edges' so that it can be removed again in the end.
        '''
        if (not DG.has_edge(startNode, endNode)):
            # edge has to be added
            # print('adding ({},{})'.format(startNode, endNode))

            DG.add_edge(startNode, endNode, weight=edgeWeight,
                        id=edgeId, oneway=direction)
            all_added_edges.append((startNode, endNode))
        else:
            # edge is already in graph
            # print('graph does have({},{})'.format(startNode,endNode))
            # print('edge:',DG.get_edge_data(startNode, endNode))
            pass

    path = None
    path_length = None
    path_IDs = []
    all_added_edges = []

    # check if the global variable 'DG' is an empty directed graph. If yes, create a graph from the shp file content.
    if(nx.is_empty(DG)):
        DG = (createGraphFromSHPInput(filepath_shp))

    # check if target lies exactly on the beginning/end of a line segment
    # if yes, target already exists as a node in the graph
    if (target not in target_line):

        # add the line_segments as edges

        d_target = LineString(target_line).project(Point(target))

        cut_line_target = cut(LineString(target_line), d_target, Point(target))
        new_line_target_after_cut = [list(x.coords) for x in cut_line_target]

        len_line_segment = len(new_line_target_after_cut[0])
        edge_target_start = new_line_target_after_cut[0][len_line_segment-2]
        edge_target_end = new_line_target_after_cut[1][1]

        # get the edge (including the attributes)

        if(target_line_oneway == 'B'):
            graph_edge = DG.get_edge_data(
                edge_target_start, edge_target_end)

        elif(target_line_oneway == 'F'):
            graph_edge = DG.get_edge_data(
                edge_target_start, edge_target_end)

        elif(target_line_oneway == 'T'):
            graph_edge = DG.get_edge_data(
                edge_target_end, edge_target_start)

        graph_edge_id = graph_edge['id']
        graph_edge_oneway = graph_edge['oneway']

        # add the line segments for the target as edges to the graph

        if(graph_edge_oneway == 'B' or ignore_oneway == True):

            # add edges in both directions
            temporarily_add_edge_to_graph(edge_target_start, target, distFrom(
                edge_target_start[0], edge_target_start[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(target, edge_target_end, distFrom(
                edge_target_end[0], edge_target_end[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(edge_target_end, target, distFrom(
                edge_target_end[0], edge_target_end[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(target, edge_target_start, distFrom(
                edge_target_start[0], edge_target_start[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

        elif(graph_edge_oneway == 'F'):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_target_start, target, distFrom(
                edge_target_start[0], edge_target_start[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(target, edge_target_end, distFrom(
                edge_target_end[0], edge_target_end[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

        else:
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(target, edge_target_start, distFrom(
                edge_target_start[0], edge_target_start[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(edge_target_end, target, distFrom(
                edge_target_end[0], edge_target_end[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

    else:
        # no need to add additional edges target
        pass

    # check if source lies exactly on the beginning/end of a line segment
    # if yes, source already exists as a node
    if (source not in source_line):
        # add the line_segments as edges

        # if source and target line are equal, the already cut line has to be used

        d_source = LineString(source_line).project(Point(source))
        cut_line_source = cut(LineString(source_line), d_source, Point(source))
        new_line_source_after_cut = [list(x.coords) for x in cut_line_source]

        len_line_segment = len(new_line_source_after_cut[0])
        edge_source_start = new_line_source_after_cut[0][len_line_segment-2]
        edge_source_end = new_line_source_after_cut[1][1]

        # get the edge (including the attributes)

        if(source_line_oneway == 'B'):
            graph_edge = DG.get_edge_data(
                edge_source_start, edge_source_end)

        elif(source_line_oneway == 'F'):
            graph_edge = DG.get_edge_data(
                edge_source_start, edge_source_end)

        elif(source_line_oneway == 'T'):
            graph_edge = DG.get_edge_data(
                edge_source_end, edge_source_start)

        graph_edge_id = graph_edge['id']
        graph_edge_oneway = graph_edge['oneway']

        # if source and target line on the same line segment, add a connection between the two
        if((source_line == target_line) and (target not in target_line) and (edge_target_start == edge_source_start) and (edge_target_end == edge_source_end)):

            # check which point comes first on the line to then know which direction the line between the two should be
            d_target_same_line = LineString(source_line).project(Point(target))
            d_source_same_line = LineString(source_line).project(Point(source))

            # get the edge (including the attributes)

            # add adges in both directions
            if(source_line_oneway == 'B' or ignore_oneway == True):
                temporarily_add_edge_to_graph(source, target, distFrom(
                    source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
                temporarily_add_edge_to_graph(target, source, distFrom(
                    source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

            # add adges in only one direction
            elif(source_line_oneway == 'F'):
                # if the source_line (which equals target_line) is oneway from the from-node to the to-node, the point which is closer to the beginning of the line has to be the starting point of the directed edge which is added to the graph
                if(d_source_same_line > d_target_same_line):
                    temporarily_add_edge_to_graph(target, source, distFrom(
                        source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
                else:
                    temporarily_add_edge_to_graph(source, target, distFrom(
                        source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

            # add adges in only one direction
            elif(source_line_oneway == 'T'):
                # if the source_line (which equals target_line) is oneway from the to-node to the from-node, the point which is closer to the end of the line has to be the starting point of the directed edge which is added to the graph
                if(d_source_same_line > d_target_same_line):
                    temporarily_add_edge_to_graph(source, target, distFrom(
                        source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)
                else:
                    temporarily_add_edge_to_graph(target, source, distFrom(
                        source[0], source[1], target[0], target[1]), graph_edge_id, graph_edge_oneway)

        # add the line segments for the source as edges to the graph

        if(graph_edge_oneway == 'B' or ignore_oneway == True):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_source_start, source, distFrom(
                edge_source_start[0], edge_source_start[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(source, edge_source_end, distFrom(
                edge_source_end[0], edge_source_end[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(edge_source_end, source, distFrom(
                edge_source_end[0], edge_source_end[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(source, edge_source_start, distFrom(
                edge_source_start[0], edge_source_start[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)

        elif(graph_edge_oneway == 'F'):
            # add edges in direction of from-node to to-node
            temporarily_add_edge_to_graph(edge_source_start, source, distFrom(
                edge_source_start[0], edge_source_start[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(source, edge_source_end, distFrom(
                edge_source_end[0], edge_source_end[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)

        else:
            # add edges in direction of to-node to from-node
            temporarily_add_edge_to_graph(edge_source_end, source, distFrom(
                edge_source_end[0], edge_source_end[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)
            temporarily_add_edge_to_graph(source, edge_source_start, distFrom(
                edge_source_start[0], edge_source_start[1], source[0], source[1]), graph_edge_id, graph_edge_oneway)

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


def calculateMostLikelyPointAndPaths(filepath, filepath_shp, minNumberOfLines=2, criticalVelocity=35.0, criticalPathLength=2.0):
    '''Maps GPS positions to the most likely points (based on the underlying street network) and obtains the most likely paths between those points based on the underlying street network.

    Parameters
    ----------
    filepath : str
        The path where the txt file (which contains the taxi mobility trace) is stored.
    filepath_shp : str
        The path where the shp file (which contains the street data) is stored.
    minNumberOfLines : int, optional
        By default 2.
    criticalVelocity : float, optional
        The critical velocity is measured in m/s.
        The critical velocity defines the critical threashhold. If for any found solution the calculated velocity is above this threashhold the algorithm checks if there is another solution which would lead to a lower velocity.
        By default 35.0 m/s. This value was chosen since cars normally do not drive that fast.
    criticalPathLength : float, optional
        The critical PathLength defines the critical threashhold. If for any found solution the pathLength devided by the airLineLength above this threashhold the algorithm checks if there is another solution which would lead to a shorter path.
        By default 2.0. This value was chosen since the data with which the algorithm was evaluated showed that values higher than 2.0 are often due due to erroneously mapped points (e.g. if a gps point is close to a crossing and is then mapped to the wrong street). By checking other solutiones when the threashold is exeeded, more likely paths can be found and thereby the accuracy of the algorithm is improved.


    Notes
    -----
    Every txt input file contains the mobility trace of a taxi. The format of each mobility trace file is the following - each line contains [latitude, longitude, occupancy, time], e.g.: [37.75134 -122.39488 0 1213084687], where latitude and longitude are in decimal degrees, occupancy shows if a cab has a fare (1 = occupied, 0 = free) and time is in UNIX epoch format.
    The output 'linestring_adjustment_visualization' can be used to visualize the mapping of points, e.g. by adding it as a layer in GQIS (Layer > Add Layer > Add delimited text layer --> Choose WKT as geometry definition and 'linestring_adjustment_visualization' as geometry field).


    Returns
    -------
    calculatedSolution : list of dictionaries
        Every line of the txt file (which was the input containing the GPS positions) yields a dictionary. This dictionary is appended to the calculatedSolution list. In the end every dictionary will be written to a new txt file as a separate line. This will yield a txt file (similar to the input file but) with additional information. The following values are taken from the input file which contains the GPS trajectories:

        - y : float
            Latitude of the GPS position in decimal degrees.

        - x : float
            Longitude of the GPS position in decimal degrees.

        - passenger : {0, 1}
            Occupancy shows if a cab has a fare (1 = occupied, 0 = free).

        - timestamp : int
            Time is in UNIX epoch format.

        In addition to those four values, every dictionary contains the following parameters which where computed by the algorithm (key : type):

        - closest_intersection_x : float
            Longitude of the mapped position based on underlying network topology.

        - closest_intersection_y : float
            Latitude of the mapped position based on underlying network topology.

        - relative_position : float
            Distance along the street (where the initial point was mapped to).

        - relative_position_normalized : float between 0 and 1
            Normalized distance along the street (where the initial point was mapped to).

        - intersected_line_oneway : {'B', 'F', 'T'}
            The 'oneway'-property of the street (LineString) which the GPS point was mapped to. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).

        - intersected_line : list of coordinates
            The source_line is a list of GPS coordinates. These coordinates represent the street on which the source point lies.

        - linestring_adjustment_visualization : list with two elements (GPS_position, mapped_point)
            The first element is the GPS position from the input file (x, y). The second element is the point on a street where the initial GPS position was mapped to (closest_intersection_x, closest_intersection_y). This list can be used to visualize the mapping of points.

        - path_time : int
            The time from the current to the next GPS position in seconds.

        - path : list of coordinates
            The shortest – most likely – path between two data point based on the streets of the underlying network.
            --> For more details on how the path is calculated, see :func:`~gps2net.getShortestPathAStar`

        - path_length : float
            The length of the path in meters. This is an approximation – see :func:`~gps2net.distFrom`.
            --> For more details on how the path is calculated, see :func:`~gps2net.getShortestPathAStar`

        - air_line_length : float
            The air line length of the current position to the next position (from start to end of 'path') in meters. This is an approximation – see :func:`~gps2net.distFrom`.

        - path_length/air_line_length : float
            'path_length' devided by 'air_line_length'

        - velocity_m_s : float
            The average velocity along the path in m/s (path_length devided by path_time).

        - pathIDs : list
            The path IDs of all street segments which are traversed on the path.
            --> For more details on how the path is calculated, see :func:`~gps2net.getShortestPathAStar`

        - solution_id : int
            The id of street where the mapped point lies on.

        - solution_index : int
            The index of the chosen solution. All possible solutions are in a list which is sorted by the distance of the solution to the GPS position. 0 means that the point was mapped to the closest intersection with a street. If this solution does not seem to be probable, other (further away) solution where checked and if one of these solutions yielded a better outcome (e.g. shorter paths) than the index of this solution is taken. Consequently, solution_index=3 means that the third-closest solution was chosen.

        - path_from_target_to_source : {0, 1}
            1 means that the on the 'oneway'-property was ignored. In this case the path might actually be the path from the target to the source. This is done when the algorithm detects a GPS glipse (the vehicle seems to drive a tiny bit backwards on a oneway street which is not possible) which resulted in a wrong path.

        - taxi_did_not_move : {0, 1}
            0 means that the taxi did move. 1 means that the taxi did not move.

        - second_best_solution_yields_more_found_paths : {0, 1}
            1 means that no path existed for initial solution and the second best solution lead to more found paths. In these cases the initial solution was replaced by the second best solution.

        - NO_PATH_FOUND : {0, 1}
            1 means that no path could be found even though the taxi moved and both source and target point are both no outliers.

        - outlier : {0, 1}
            1 means that the GPS position is marked as an outlier. This is done when no street could be found in a specified area.

        - comment : str
            Comment which summarizes most important assumptions/results/concerns for a specific GPS position.

    statistics : list
        Statistics for the txt file (containing GPS positions) for which the solution is calculated.
        The statistics contain the following numbers (all of type int):

        - outlier : int
            Number of data points which where flagged as outliers.
        - taxi_did_not_move : int
            Number of times when a data point's GPS position was identical to the GPS position of the previous data point.
        - no_path_found : int
            Number of times no path was found even though the taxi moved and both source and target point are both no outliers.
        - cannot_compute_shortest_path_as_previous_point_is_outlier : int
            Number of times the path could not be computed since the previous point was an outlier.
        - path_from_target_to_source : int
            Number of times the 'oneway'-property was ignored. In this case the path might actually be the path from the target to the source. This is done when the algorithm detects a GPS glipse (the vehicle seems to drive a tiny bit backwards on a oneway street which is not possible) which resulted in a wrong path.
        - checked_other_solution_index : int
            Number of times another solution was checked because a path either lead to a velocity of more than 35 m/s or a path length which is more than double the air line length.
        - chose_other_solution_index : int
            Number of times the other solution actually lead to a shorter cummulated path (path to previous point plus path to next point) than the initial solution
        - solution_already_lies_on_shortest_path : int
            Number of times when velocity was more than 35 m/s or path length was is more than double the air line length BUT the initial solution was already on the shortest path from the previous point to the next point (in this case, the solution is optimal which is why no other solution was checked).
        - no_solution_lies_on_shortest_path : int
            Number of times when velocity was more than 35 m/s or path length was is more than double the air line length BUT there was no other solution which lied on the shortest path from the previous point to the next point (in this case, no other solution could b checked).
        - other_solution_is_worse : int
            Number of times when another solution was checked, but since this new solution was worse than the initial solution it wasn't taken.
        - other_solution_no_path_found : int
            Number of times when another solution was checked, but since there was no valid path for this new solution it wasn't taken.

        - checked_if_path_exists_for_second_best_solution_index : int
            Number of times no path existed for initial solution. In these cases the secont best solution (if it existed) was checked.
        - second_best_solution_yields_more_found_paths : int
            Number of times no path existed for initial solution and the second best solution lead to more found paths. In these cases the initial solution was replaced by the second best solution.

    '''

    def getLocationResult(filepath_shp, x, y, passenger, timestamp, previous_point, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway, minNumberOfLines=2):
        """Calculates all additional parameters for one GPS position (which corresponds to one line in the txt file).

        Parameters
        ----------
        filepath_shp : str
            The path where the shp file (which contains the street data) is stored.
        y : float
            Latitude of the GPS position in decimal degrees.
        x : float
            Longitude of the GPS position in decimal degrees.
        passenger : {0, 1}
            Occupancy shows if a cab has a fare (1 = occupied, 0 = free).
        timestamp : int
            Time is in UNIX epoch format.
        previous_point : tuple of float
            GPS position of the mapped previous point based on the underlying street structure.
        previous_intersected_line : list of coordinates
            The previous_intersected_line is a list of GPS coordinates. These coordinates represent the street on which the mapped previous point lies.
        previous_timestamp : int
            Time of the previous point (in UNIX epoch format).
        previous_intersected_line_oneway : [list of coordinates
            The 'oneway'-property of the previous_intersected_line. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).
        minNumberOfLines : int, optional
            The min number of lines is just needed to increse areasize. This means that the areasize is increased as long as minNumberOfLines is not found and max areasize is not exeeded.
            A GPS position is only marked as an outlier if NOT EVEN A SINGLE STREET could be found in the area
            By default 2.

        Notes
        -----
        This function is only called in the following function: :func:`~gps2net.calculateMostLikelyPointAndPaths`

        The additional parameters are calculated in this function. However, the solution is not saved. When the next GPS position is looked at, the result is validated and only if then it still found to be the most likely result is saved as such.
        However, if then another result (namely, mapping the GPS position to another street which yields a result which lies further away from the GPS location) is found to be more likely, the initial result is overridden.

        Returns
        -------
        location_result : dict
            This dict contains the main solution. More specific, the additional attribute which will be written to the output file (after validation).
        (closest_intersection_x, closest_intersection_y) : tuple of float
            GPS position of the mapped point based on the underlying street structure.
        previous_point : tuple of float
            GPS position of the mapped previous point based on the underlying street structure.
        intersected_line : list of coordinates
            The intersected_line is a list of GPS coordinates. These coordinates represent the street on which the mapped point (closest_intersection_x, closest_intersection_y) lies.
        previous_intersected_line : list of coordinates
            The previous_intersected_line is a list of GPS coordinates. These coordinates represent the street on which the mapped previous point lies.
        timestamp : int
            Time is in UNIX epoch format.
        intersected_line_oneway : {'B', 'F', 'T'}
            The 'oneway'-property of the intersected_line. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).
        previous_intersected_line_oneway : list of coordinates
            The 'oneway'-property of the previous_intersected_line. The 'oneway'-property indicates if a street is bi-directional (B), or one way heading from the from-node to the to-node (F), or one way heading from the to-node to the from-node (T).
        """

        location_result = {}
        # all the paths who were not changed from source-->target to target-->source have the following property
        path_from_target_to_source = 0

        closest_intersection_x = 0
        closest_intersection_y = 0
        intersected_line = None
        intersected_line_oneway = None

        # initiaize all keys (as some of them are not necessarily set, e.g. when previous point was an outlier)

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
        location_result['taxi_did_not_move'] = 0
        location_result['second_best_solution_yields_more_found_paths'] = 0
        location_result['NO_PATH_FOUND'] = 0
        location_result['outlier'] = 0

        comment = ''

        location_result['x'] = x
        location_result['y'] = y
        location_result['source'] = (x, y)
        location_result['passenger'] = passenger
        location_result['timestamp'] = timestamp
        location_result['previous_timestamp'] = previous_timestamp

        # how to get the closest line from one point

        gps_point = Point(x, y)  # radial sweep centre point

        with fiona.open(filepath_shp) as street_lines:
            # define size of area to filter streets. The areasize is measured in decimal degrees.
            # an areasize of 0.001 is (in the area of San Francisco) approximately 88 meters for longitudes and approximately 111 meters for latitude. This means that if the areasize s 0.001, streets are searched in an area of approximately 196 (longitude) times 222 meters (latitute), since the aresize is added in all four directions of a GPS position. Or an even broader approximation:
            # an areasize of 0.001 means that an area of approx. 200x200 meters is covered.
            # an areasize of 0.0005 means that an area of approx. 50x50 meters is covered.
            areasize = 0.0005
            max_areasize = 0.001

            # the min number of lines is just needed to increse areasize. This means that the areasize is increased as long as minNumberOfLines is not found and max areasize is not exeeded.
            # a GPS position is only marked as an outlier if NOT EVEN A SINGLE STREET could be found in the area
            number_of_streets = 0

            while ((number_of_streets < minNumberOfLines) and (areasize < max_areasize)):

                input_shapes = list(street_lines.items(
                    bbox=((x-areasize), (y-areasize), (x+areasize), (y+areasize))))

                number_of_streets = len(input_shapes)

                # double are size so that if no lines are found, a bigger area is taken into account when filtering
                areasize = areasize*1.25

        # point is oulier only if not even a single street was found
        if(number_of_streets == 0):
            comment += 'This point is an outlier. '

            # makes sure that next point doesn't calculate shortest path
            previous_point = (0, 0)

            # set the outlier property of the solution to 1
            location_result['outlier'] = 1

        # point is not an outlier
        else:

            solution_dict = {}

            for input_line in input_shapes:

                lineID = input_line[1]['id']

                newLS = LineString(input_line[1]['geometry']['coordinates'])
                # get the distance along the LineString to a point nearest to the point
                relative_position = newLS.project(gps_point)
                # get the distance normalized to the length of the LineString
                relative_position_normalized = newLS.project(
                    gps_point, normalized=True)

                closest_point_on_line = newLS.interpolate(relative_position)

                inter_dict_point = {}

                inter_dict_point['closest_point_on_line'] = closest_point_on_line
                inter_dict_point['relative_position_closest_point_on_line'] = relative_position
                inter_dict_point['relative_position_normalized_closest_point_on_line'] = relative_position_normalized
                inter_dict_point['input_line'] = newLS

                inter_dict_point['oneway'] = input_line[1]['properties']['oneway']

                # distance needs to be float so that it can be sorted appropriately afterwards
                inter_dict_point['distance'] = float(
                    gps_point.distance(closest_point_on_line))

                solution_dict[lineID] = inter_dict_point

            # sort the nested dict by the 'distance' which lies within all the inner dicts
            # the first element in the sorted dict is the point with the shortest distance to the gps_point
            solution_dict_sorted_by_distance = sorted(
                solution_dict.items(), key=lambda kv: (kv[1]['distance']))

            # get the first element as it is the one with the smallest distance (which is the closest point)
            solution_id = solution_dict_sorted_by_distance[0][0]
            solution = solution_dict_sorted_by_distance[0][1]
            location_result['all_solutions_sorted'] = solution_dict_sorted_by_distance
            location_result['solution_id'] = solution_id
            location_result['solution'] = solution

            closest_intersection_x = solution['closest_point_on_line'].x
            closest_intersection_y = solution['closest_point_on_line'].y
            relative_position = solution['relative_position_closest_point_on_line']
            relative_position_normalized = solution['relative_position_normalized_closest_point_on_line']

            intersected_line = solution['input_line']
            intersected_line_oneway = solution['oneway']

            # append x- and y-position of the mapped point, i.e. the closest point which lies on a street of the underlying network
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
                [(gps_point.x, gps_point.y), (closest_intersection_x, closest_intersection_y)])

            location_result['linestring_adjustment_visualization'] = linestring_adjustment_visualization

            # check if previous point is set. if yes: calculate the shortest path from current point to previous point
            if(previous_point == (0, 0)):
                # no previous point set, do not calculate shortest path.
                comment += 'Cannot compute shortest path as previous point is (0,0). '

                # update statistics
                statistics['cannot_compute_shortest_path_as_previous_point_is_outlier'] += 1

            elif((closest_intersection_x, closest_intersection_y) == previous_point):
                # if the taxi did not move, path is not created

                comment += 'source==target --> Taxi did not move. '

                location_result['path_length'] = 0
                location_result['air_line_length'] = 0
                location_result['path_length/air_line_length'] = 1
                location_result['velocity_m_s'] = 0
                location_result['taxi_did_not_move'] = 1

                # update statistics
                statistics['taxi_did_not_move'] += 1

            else:

                path_time = abs(previous_timestamp-timestamp)

                location_result['path_time'] = path_time

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

                    # calculate where the source lies on the line
                    d_source = LineString(list(intersected_line.coords)).project(
                        Point(closest_intersection_x, closest_intersection_y))

                    # calculate where the target lies on the line
                    d_target = LineString(list(previous_intersected_line.coords)).project(
                        Point(previous_point))

                    # check if there was a glipse in the gps coordinates which resulted in a situation where the source lies after the target on a oneway street
                    # if it is a oneway line from the FROM-node to the TO-node and source lies AFTER the target   OR   if it is a oneway line from the TO-node to the FROM-node and source lies BEFORE the target
                    if((intersected_line_oneway == 'F' and d_source > d_target)or(intersected_line_oneway == 'T' and d_source < d_target)):

                        # calculate the shortest path with A STAR algorithm
                        # set ignore_oneway=True : this makes sure that the source_line and the target_line both are treated as lines where driving in both directions is allowed (so feature 'oneway' is ignored)
                        path2, path_length2, pathIDs2 = getShortestPathAStar((closest_intersection_x, closest_intersection_y), previous_point, list(intersected_line.coords), list(
                            previous_intersected_line.coords), intersected_line_oneway, previous_intersected_line_oneway, filepath_shp, ignore_oneway=True)

                    # if the new solution is more likely to be correct, override path/path_lenght/pathIDs to make sure to use the new solution, i.e. the one which ignores the 'oneway' property
                    if((path == None and path2 != None) or (path != None and path2 != None and (path_length > path_length2))):
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
                    printComment = 'No path was found with A Star algorithm from {} to {}. '.format(
                        (closest_intersection_x, closest_intersection_y), previous_point)

        location_result['path_from_target_to_source'] = path_from_target_to_source
        location_result['comment'] = comment

        return location_result, (closest_intersection_x, closest_intersection_y), previous_point, intersected_line, previous_intersected_line, timestamp, intersected_line_oneway, previous_intersected_line_oneway

    global current_txt_file
    global number_of_txt_files

    counter = 0

    calculatedSolution = []

    # initialize the dict in which all the statistics of the currently calculated solution will be stored.
    statistics = {}
    statistics['outlier'] = 0
    statistics['taxi_did_not_move'] = 0
    statistics['no_path_found'] = 0
    statistics['cannot_compute_shortest_path_as_previous_point_is_outlier'] = 0
    statistics['path_from_target_to_source'] = 0
    statistics['checked_other_solution_index'] = 0
    statistics['chose_other_solution_index'] = 0
    statistics['solution_already_lies_on_shortest_path'] = 0
    statistics['no_solution_lies_on_shortest_path'] = 0
    statistics['other_solution_is_worse'] = 0
    statistics['other_solution_no_path_found'] = 0

    statistics['checked_if_path_exists_for_second_best_solution_index'] = 0
    statistics['second_best_solution_yields_more_found_paths'] = 0

    previous_target = (0, 0)
    previous_source = (0, 0)
    previous_intersected_line = None
    previous_timestamp = None
    previous_intersected_line_oneway = None

    # calculate the number of lines of the file we are looking at
    lines_in_textfile = 0
    with open(filepath, 'r') as f:
        lines_in_textfile = sum(1 for line in f)

    # Initial call to print 0% progress ProgressBar
    suffix = '| current file: {}/{} lines'.format(counter+1, lines_in_textfile)
    suffix += ' | total: {} of {} files'.format(
        current_txt_file, number_of_txt_files)

    previous_location_result = {}

    with open(filepath, 'r') as f:
        mylist = f.read().splitlines()

        # append an artificial line at the end. This will make sure that also the path from the last to the second-last element will be calculated and possible improvements will be checked.
        mylist.append('artificialline')

        for txt_line in mylist:

            # the last line of the text file is artificial as it was added before. this line is just needed to make sure the second last line (which is actually the last actual line of the text file is calculated correctly)
            if(txt_line != 'artificialline'):

                values = txt_line.split(' ')
                # the longitude, latitude, passenger and timestamp properties were read from the txt file an saved as a string, here we transform it to a float / int
                x = float(values[1])
                y = float(values[0])
                passenger = int(values[2])
                timestamp = int(values[3])

                # get the result for the current location
                current_location_result, source, target, intersected_line, target_intersected_line, timestamp, intersected_line_oneway, target_intersected_line_oneway = getLocationResult(
                    filepath_shp, x, y, passenger, timestamp, previous_source, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway, minNumberOfLines)

            if (previous_location_result != {}):

                # check if the result of the previous location can be improved
                # we check if either velocity > criticalVelocity or if path_length/air_line_length > criticalRPathLength
                # by default the criticalVelocity is 35.0 m/s
                # by default the criticalPathLength is 2.0. This means that other solutions are checked if the calculated path length is >= the double of the airLineLength

                # no point is an outlier and all paths were found --> check if the solution can be improved
                if(previous_target != (0, 0) and target != (0, 0) and previous_location_result['path'] != ''):

                    # check if either the path from current source to target or from previous source to previous target might be needed to improved
                    # if (current_location_result['velocity_m_s']=='') means that the taxi was not moving at the previous point --> for this reason it should not be checked

                    if(current_location_result['velocity_m_s'] != '' and ((float(previous_location_result['velocity_m_s']) > criticalVelocity) or ((previous_location_result['path_length/air_line_length']) > criticalPathLength) or (float(current_location_result['velocity_m_s']) > criticalVelocity) or ((current_location_result['path_length/air_line_length']) > criticalPathLength))):

                        # update statistics
                        statistics['checked_other_solution_index'] += 1

                        # the goal is to check if there is another possible solution and then compare this new solution with the current solution
                        # Let's calculate shortest path between current source and previous target

                        # first calculate shortest path between current source and previous target <-- this will be the baseline so to say; We are trying to improve this solution, namely we try to shorten this path.
                        path_to_previous_target, path_length_to_previous_target, pathIDs_to_previous_target = getShortestPathAStar(source, previous_target, list(
                            intersected_line.coords), list(previous_target_intersected_line.coords), intersected_line_oneway, previous_target_intersected_line_oneway, filepath_shp)

                        # check if the initial solution already lies on a street segment which is part of the hortest path between the points before and after
                        if(str(previous_location_result['solution_id']) in pathIDs_to_previous_target or previous_location_result['solution_id'] in pathIDs_to_previous_target):
                            blockPrint()
                            print(
                                'Solution is already on the shortest path between the points before and after --> therefore the solution can not be improved.')
                            print('')
                            enablePrint()

                            # update statistics
                            statistics['solution_already_lies_on_shortest_path'] += 1

                        else:
                            # Try to find solution on shortest path.

                            new_solution_id = (-1)
                            new_solution_index = 0
                            new_solution = {}

                            # try to find another solution point which lies on the shortest path between the current source and the previous target
                            for k, v in previous_location_result['all_solutions_sorted']:
                                if(str(k) in pathIDs_to_previous_target):
                                    new_solution_id = k
                                    new_solution = v
                                new_solution_index += 1

                            if (new_solution_id == (-1)):
                                # There is no other solution which lies on the shortest path between the current source and the previous target.

                                # update statistics
                                statistics['no_solution_lies_on_shortest_path'] += 1

                            else:
                                # there is another solution point which lies on the shortest path between the current source and the previous target --> check if this solution is actually better

                                previous_commulated_path_lengt = float(
                                    previous_location_result['path_length']) + float(current_location_result['path_length'])

                                new_solution_point_x = list(
                                    new_solution['closest_point_on_line'].coords)[0][0]
                                new_solution_point_y = list(
                                    new_solution['closest_point_on_line'].coords)[0][1]

                                # calculate the solution for the data point whe currently look at with the new solution point.
                                current_location_result_new, source_new, target_new, intersected_line_new, target_intersected_line_new, timestamp_new, intersected_line_oneway_new, target_intersected_line_oneway_new = getLocationResult(
                                    filepath_shp, x, y, passenger, timestamp, (new_solution_point_x, new_solution_point_y), new_solution['input_line'], previous_location_result['timestamp'], new_solution['oneway'], minNumberOfLines)

                                # calculate the solution for the data point whe previously looked at with the new solution point.
                                previous_location_result_new, previous_source_new, previous_target_new, previous_intersected_line_new, previous_target_intersected_line_new, previous_timestamp_new, previous_intersected_line_oneway_new, previous_target_intersected_line_oneway_new = getLocationResult(
                                    filepath_shp, new_solution_point_x, new_solution_point_y, previous_location_result['passenger'], previous_location_result['timestamp'], previous_location_result['target'], previous_location_result['previous_intersected_line'], previous_location_result['previous_timestamp'], previous_location_result['previous_intersected_line_oneway'], minNumberOfLines)

                                # check if with the new solution both paths (from new solution to previous and to next point) exist.
                                if(previous_location_result_new['path_length'] != '' and current_location_result_new['path_length'] != ''):

                                    new_commulated_path_lengt = float(
                                        previous_location_result_new['path_length']) + float(current_location_result_new['path_length'])

                                    # if the new cummulated path is shorter set it as the correct one
                                    if((previous_commulated_path_lengt > new_commulated_path_lengt)):
                                        blockPrint()
                                        print_str = 'Set the second checked solution (with solution_id: {}) as the correct one since the new cummulated path was shorter than the initial solution (with solution_id: {}). '.format(
                                            new_solution_id, previous_location_result['solution_id'])
                                        print(print_str)
                                        print('')
                                        enablePrint()
                                        new_comment = 'The source point of this path was reset as this leads to much shorter path lengths. Previous solution: {} / New solution: {}. '.format(
                                            previous_source, (new_solution_point_x, new_solution_point_y))

                                        # save the initial GPS position coordinates
                                        myLat = previous_location_result['y']
                                        myLng = previous_location_result['x']

                                        # update the previous location result
                                        previous_location_result = previous_location_result_new

                                        # previous_location_result_new saved the updated x and y position. For this reason the initial GPS position coordinates have to be set again.
                                        previous_location_result['y'] = myLat
                                        previous_location_result['x'] = myLng

                                        # previous_location_result_new used the updated x and y position for linestring_adjustment_visualization. For this reason the initial GPS position coordinates have to be set again.
                                        # append the line which visualizes the way from the start point to the new position of the closest intersection point
                                        adjustment_visualization_new = LineString(
                                            [(myLng, myLat), (previous_location_result['closest_intersection_x'], previous_location_result['closest_intersection_y'])])

                                        previous_location_result['linestring_adjustment_visualization'] = adjustment_visualization_new

                                        previous_location_result['comment'] += new_comment
                                        previous_location_result['solution_id'] = new_solution_id
                                        previous_location_result['solution_index'] = new_solution_index

                                        current_location_result = current_location_result_new
                                        current_location_result['comment'] += 'The target point of this path was reset as this leads to much shorter path lengths.'

                                        # update statistics
                                        statistics['chose_other_solution_index'] += 1

                                    else:
                                        new_comment = 'Checked the second solution (solution_id: {} / solution_index: {} / new cummulated path_length: {}) but initial solution (solution_id: {} / solution_index: {} / initial cummulated path_length: {}) was chosen due to a shorter cummulated path_length.'.format(
                                            new_solution_id, new_solution_index, new_commulated_path_lengt, previous_location_result['solution_id'], previous_location_result['solution_index'], previous_commulated_path_lengt)
                                        blockPrint()
                                        print(new_comment)
                                        print('')
                                        enablePrint()
                                        previous_location_result['comment'] += new_comment

                                        # update statistics
                                        statistics['other_solution_is_worse'] += 1

                                else:
                                    # no path was found with new solution
                                    new_comment = 'No path found for new solution. '
                                    previous_location_result['comment'] += new_comment

                                    # update statistics
                                    statistics['other_solution_no_path_found'] += 1

                # no path was found, either from current source to target (current_location_result['path']) or from previous source to previous target (previous_location_result['path']), even though source and target are known --> check if a path can be found with another solution
                elif(previous_target != (0, 0) and ((previous_location_result['path'] == '' and previous_location_result['taxi_did_not_move'] != 1 and previous_location_result['outlier'] != 1) or (current_location_result['path'] == '' and current_location_result['taxi_did_not_move'] != 1 and current_location_result['outlier'] != 1)) and len(previous_location_result['all_solutions_sorted']) > 1):

                    # check how many paths could not be found --> this value is either 1 or 2 and will be our baseline so to say to check whether the new solution will yield more found paths
                    not_found_paths = 0
                    if (previous_location_result['path'] == '' and previous_location_result['taxi_did_not_move'] != 1 and previous_location_result['outlier'] != 1):
                        not_found_paths += 1
                    if (current_location_result['path'] == '' and current_location_result['taxi_did_not_move'] != 1 and current_location_result['outlier'] != 1):
                        not_found_paths += 1

                    # update statistics
                    statistics['checked_if_path_exists_for_second_best_solution_index'] += 1

                    # get the second best solution
                    new_solution_index = 1
                    new_solution_id = previous_location_result['all_solutions_sorted'][1][0]
                    new_solution = previous_location_result['all_solutions_sorted'][1][1]

                    new_solution_point_x = list(
                        new_solution['closest_point_on_line'].coords)[0][0]
                    new_solution_point_y = list(
                        new_solution['closest_point_on_line'].coords)[0][1]

                    # calculte current location result with new solution
                    current_location_result_new, source_new, target_new, intersected_line_new, target_intersected_line_new, timestamp_new, intersected_line_oneway_new, target_intersected_line_oneway_new = getLocationResult(
                        filepath_shp, x, y, passenger, timestamp, (new_solution_point_x, new_solution_point_y), new_solution['input_line'], previous_location_result['timestamp'], new_solution['oneway'], minNumberOfLines)

                    # calculte previous location result with new solution
                    previous_location_result_new, previous_source_new, previous_target_new, previous_intersected_line_new, previous_target_intersected_line_new, previous_timestamp_new, previous_intersected_line_oneway_new, previous_target_intersected_line_oneway_new = getLocationResult(
                        filepath_shp, new_solution_point_x, new_solution_point_y, previous_location_result['passenger'], previous_location_result['timestamp'], previous_location_result['target'], previous_location_result['previous_intersected_line'], previous_location_result['previous_timestamp'], previous_location_result['previous_intersected_line_oneway'], minNumberOfLines)

                    # check how many paths could not be found WITH NEW SOLUTION
                    not_found_paths_NEW_SOLUTION = 0
                    if (previous_location_result_new['path'] == ''):
                        not_found_paths_NEW_SOLUTION += 1
                    if (current_location_result_new['path'] == ''):
                        not_found_paths_NEW_SOLUTION += 1

                    # check if the new solution is actually better than the old one, namely if more paths could be found.
                    if(not_found_paths_NEW_SOLUTION < not_found_paths):

                        print_str = 'Set the second checked solution (with solution_id: {}) as the correct one since it yields more found paths than the initial solution (with solution_id: {}). '.format(
                            new_solution_id, previous_location_result['solution_id'])
                        blockPrint()
                        print(print_str)
                        print('')
                        enablePrint()
                        new_comment = 'The source point of this path was reset as this yields more found paths. Previous solution: {} / New solution: {}. '.format(
                            previous_source, (new_solution_point_x, new_solution_point_y))

                        # save the initial GPS position coordinates
                        myLat = previous_location_result['y']
                        myLng = previous_location_result['x']

                        # update the previous location result
                        previous_location_result = previous_location_result_new

                        # previous_location_result_new saved the updated x and y position. For this reason the initial GPS position coordinates have to be set again.
                        previous_location_result['y'] = myLat
                        previous_location_result['x'] = myLng

                        previous_location_result['comment'] += new_comment
                        previous_location_result['solution_id'] = new_solution_id
                        previous_location_result['solution_index'] = new_solution_index

                        current_location_result = current_location_result_new
                        current_location_result['comment'] += 'The target point of this path was reset as this yields more found paths.'

                        # update statistics
                        statistics['second_best_solution_yields_more_found_paths'] += 1
                        # update statistic in output txt file
                        previous_location_result['second_best_solution_yields_more_found_paths'] = 1

                    else:
                        # new solution does not yield more found paths than initial solution
                        blockPrint()
                        print(
                            'New solution does not yield more found paths than initial solution.')
                        new_comment_second_best_solution = 'Checked the second best solution (solution_id: {} / solution_index: {}) but it did not lead to more found paths. Therefore, initial solution (solution_id: {} / solution_index: {}) was chosen. '.format(
                            new_solution_id, new_solution_index, previous_location_result['solution_id'], previous_location_result['solution_index'])
                        print(new_comment_second_best_solution)
                        print('')
                        enablePrint()
                        previous_location_result['comment'] += new_comment_second_best_solution

                else:
                    blockPrint()
                    print('Tried to check other solution but the target is not known. ')
                    print('')
                    enablePrint()

                # append the result of the previous location to the entire solution
                calculatedSolution.append(previous_location_result)
                if (previous_location_result['outlier'] == 1):

                    # update statistics
                    statistics['outlier'] += 1

                # check if a path was found
                if (previous_location_result['path'] == '' and previous_target != (0, 0) and previous_location_result['taxi_did_not_move'] != 1):
                    # update statistics
                    statistics['no_path_found'] += 1
                    # add comment
                    new_comment = 'No path found from {} to {}. '.format(
                        (previous_location_result['closest_intersection_x'], previous_location_result['closest_intersection_y']), previous_location_result['target'])
                    previous_location_result['comment'] += new_comment
                    # update statistic in output txt file
                    previous_location_result['NO_PATH_FOUND'] = 1

            # set the result of the current location as the result of the previous location
            previous_location_result = current_location_result

            previous_source = source
            previous_target = target

            previous_intersected_line = intersected_line
            previous_target_intersected_line = target_intersected_line
            previous_timestamp = timestamp
            previous_intersected_line_oneway = intersected_line_oneway
            previous_target_intersected_line_oneway = target_intersected_line_oneway

            # The following lines are just for testing. Uncomment 'break' in order to only run the algorithm for a certain amount of lines per file. The counter counts the amount of lines of the current txt file which were already calculated.
            counter += 1
            if (counter > 50):
                # break
                pass

            # Update Progress Bar
            if(txt_line != 'artificialline'):
                suffix = '| current file: {}/{} lines'.format(
                    counter, lines_in_textfile)
                suffix += ' | total: {} of {} files'.format(
                    current_txt_file, number_of_txt_files)
                printProgressBar(counter, lines_in_textfile,
                                 prefix='Progress:', suffix=suffix, length=50)

    return calculatedSolution, statistics


def getTimeDifferences(filepath, timestampPosition):
    """Get the time differences of taxi mobility traces in a txt file.

    Parameters
    ----------
    filepath : str
        The path where the txt file (which contains the taxi mobility trace) is stored.
    timestampPosition : int
        The position on which the timestamp is (assuming that each line contains one measured GPS position and that values in each line are seperated by a space).

    Notes
    -----
    This method assumes that each line contains one measured GPS position and that values in each line are seperated by a space. Further, it assumes that time is in UNIX epoch format.
    Examples
    --------
    >>> myPath = '/Users/Joechi/Google Drive/gps2net/Data/testData/testTaxi.txt'
    >>> myTimeDifferences = getTimeDifferences(myPath, 3)
    >>> myTimeDifferences
    [28, 119]

    Returns
    -------
    list
        Time differences in seconds.
    """

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
    """Plots a histogram with a specified max x-value. Values bigger that the max x-value are clipped to the interval edge.

    Parameters
    ----------
    input : list of int or float
        Values which should be plottet in a histogram.
    minXLabel : int or float
        The min x-value of the plot.
    maxXLabel : [type]
        The max x-value of the plot.
    binSize : int
        The size of the plotted bins.
    filename : str
        Path and filename of the plot --> this specifies where and under which name the plot should be saved.
    title : str
        Title of the plot.
    xlabel : str
        X axis label of the plot.

    Notes
    -----
    The plot is saved in PNG format. 
    """
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

    try:
        maxfreq = n.max()
        # Set a clean upper y-axis limit.
        plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq %
                 10 else maxfreq + 10)

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

    except ValueError:  # raised if `n` is empty.
        pass

    # save the figure in the folder of the current taxi
    plt.savefig(filename)
    # the figure could be shown with 'plt.show()'. However, this is not necessary as it is saved anyway.
    # close the current figure
    plt.close()


def getPathFromUnmappedGpsPositions(filepath, new_filename):
    """Creates a new file which includes path, path length, path time, and velocity without considering the underlying street network.

    Parameters
    ----------
    filepath : str
        The path where the txt file (which contains the taxi mobility trace) is stored.
    new_filename : str
        The filename of the new solution.

    Notes
    -----
    This function created a txt file which (in addition to the input values) contains additional parameters which are NOT based on the underlying street network, such as:

    - path : Linestring
        Linear path with start=current_position and end=next_position.
    - path length : float
        Distance between current and next position in meters.
    - path time : float
        Time between measurement at current and next position in seconds.
    - velocity : float
        Velocity in meters/second.
    """

    header = [
        "latitude(y);longitude(x);hasPassenger;time;path;path_length;path_time;velocity_m_s\n"]

    path = ''
    path_length = ''
    path_time = ''
    velocity_m_s = ''
    time_previous = None

    with open(filepath, "r") as f:
        mylist = f.read().splitlines()

        # this saves a new text file which includes the calculated parameters
        with open(new_filename, 'w') as new_file:

            # write the header to the new txt file
            new_file.writelines(header)

            for line in mylist:

                # read the values from the txt file and transform the strings to float/int

                values = line.split(" ")
                x = float(values[1])
                y = float(values[0])
                time_current = int(values[3])
                lines_withoutspaces = str(line).replace(" ", ";")

                # no values are added to the first line of the file
                if (time_previous != None):

                    # check if data is sorted by time
                    # we start with newest data point, so if previous data time was earlier in time, it is not sorted
                    if (time_previous < time_current):
                        print('')
                        print('Algorithm cannot run since DATA IS NOT SORTED!')
                        print('Time of current position ({}) is before time of previous position ({}).'.format(
                            time_previous, time_current))
                        print(
                            'Please adjust your input file. The newest gps signal should be the first line.')
                        print('')
                        sys.exit('Execution of the algorithm stopped.')

                    else:
                        start_pt = Point(x_previous, y_previous)
                        end_pt = Point(x, y)
                        path = LineString(
                            [(start_pt.x, start_pt.y), (end_pt.x, end_pt.y)])
                        path_length = distFrom(x_previous, y_previous, x, y)
                        path_time = abs(time_previous-time_current)

                        velocity_m_s = path_length/path_time

                # write the initial values to the new file
                new_file.write(lines_withoutspaces)

                # write the new values to the new file
                new_file.write(';')
                new_file.write(str(path))
                new_file.write(';')
                new_file.write(str(path_length))
                new_file.write(';')
                new_file.write(str(path_time))
                new_file.write(';')
                new_file.write(str(velocity_m_s))
                new_file.write('\n')

                x_previous = x
                y_previous = y
                time_previous = time_current


def getFilename(path):
    '''Returns the name of a file from a specific path (excluding directories and extension).

    Parameters
    ----------
    path : str
        The path of a specific data file.

    Returns
    -------
    str
        The name of the file. Namely, the last part of the path (excluding the extgension).


    Examples
    --------

    >>> myPath = 'dir1/dir2/ThisIsMyFilename.txt'
    >>> filename = getFilename(myPath)
    >>> filename
    'ThisIsMyFilename'

    '''
    # get the name of the file from the path (without the directories and without extension)
    filename = os.path.splitext(os.path.basename(os.path.normpath(path)))[0]
    return filename


def caculationForOneTXTFile(filepath_shp, new_filename_solution, new_filename_statistics, new_filename_velocities, new_filename_path_length_air_line_length, filepath):
    """Calculates the most likely solution for one entire txt file based on the underlying street network and saves the solution.

    Parameters
    ----------
    filepath_shp : str
        The path where the shp file (which contains the street data) is stored.
    new_filename_solution : str
        The filename of the new solution.
    new_filename_statistics : str
        The filename of the new solution.
    new_filename_velocities : str
        The filename of the velocities histogram.
    new_filename_path_length_air_line_length : str
        The filename of the path_length_air_line_length histogram.
    filepath : str
        The path where the txt file (which contains the taxi mobility trace) is stored.
    """

    # set the header of the output txt file which will contain the calculated solution.
    header = ['latitude(y);longitude(x);hasPassenger;time;closest_intersection_x;closest_intersection_y;relative_position;relative_position_normalized;intersected_line_oneway;intersected_line_as_linestring;linestring_adjustment_visualization;path_time;path_as_linestring;path_length;air_line_length;path_length/air_line_length;velocity_m_s;pathIDs;solution_id;solution_index;path_from_target_to_source;taxi_did_not_move;second_best_solution_yields_more_found_paths;NO_PATH_FOUND;outlier;comment\n']

    myCalculatedSolution, mySolutionStatistics = calculateMostLikelyPointAndPaths(
        filepath, filepath_shp, minNumberOfLines=2, criticalVelocity=35.0, criticalPathLength=2.0)

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
            new_file.write(str(location_result['taxi_did_not_move']))
            new_file.write(';')
            new_file.write(
                str(location_result['second_best_solution_yields_more_found_paths']))
            new_file.write(';')
            new_file.write(str(location_result['NO_PATH_FOUND']))
            new_file.write(';')
            new_file.write(str(location_result['outlier']))
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


def main():

    global number_of_txt_files
    global current_txt_file

    filepath_shp = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/San Francisco Basemap Street Centerlines/geo_export_e5dd0539-2344-4e87-b198-d50274be8e1d.shp'

    filepaths = []

    filepath1 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi1.txt'
    filepath2 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi2.txt'
    filepath3 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi3.txt'
    filepath4 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi4.txt'
    filepath5 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi5.txt'
    filepath6 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi6.txt'
    filepath7 = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/cabspottingdata/taxi7.txt'

    filepaths.append(filepath1)
    filepaths.append(filepath2)
    filepaths.append(filepath3)
    filepaths.append(filepath4)
    filepaths.append(filepath5)
    filepaths.append(filepath6)
    filepaths.append(filepath7)

    number_of_txt_files = len(filepaths)

    # loop through all the filepaths
    for path in filepaths:

        current_txt_file += 1

        new_filename = getFilename(path)

        dirName = os.path.join('output_files', new_filename)

        # Create target directory & all intermediate directories if don't exists
        try:
            os.makedirs(dirName)
            print('')
            print('Directory ', dirName,  ' created.')
        except FileExistsError:
            print('')
            print('Directory ', dirName,  ' already exists.')

        new_filename_simple_solution = os.path.join(
            dirName, 'pathFromUnmappedGpsPositions') + '.txt'

        new_filename_solution = os.path.join(
            dirName, 'calculatedSolution') + '.txt'
        new_filename_statistics = os.path.join(dirName, 'statistics') + '.txt'
        new_filename_velocities = os.path.join(dirName, 'velocitiesPLOT.png')
        new_filename_path_length_air_line_length = os.path.join(
            dirName, 'path_length_air_line_length_PLOT.png')

        # calculate and save the simple solution (path from exact gps positions without considering the underlying street network)
        getPathFromUnmappedGpsPositions(path, new_filename_simple_solution)

        # calculate and save the full solution (most likely paths) based on the underlying street network
        caculationForOneTXTFile(filepath_shp, new_filename_solution, new_filename_statistics,
                                new_filename_velocities, new_filename_path_length_air_line_length, path)

        # get all timestamp differences of a text file
        timeDifferences = getTimeDifferences(path, 3)

        timedifferencesFileName = os.path.join(
            dirName, 'timedifferencesPLOT.png')

        # plot the timeDifferences in a histogram
        plotAndSaveHistogram(timeDifferences, 0, 300, 25, timedifferencesFileName,
                             'Histogram of time differences between gps points', 'time difference in seconds')

        print('')
        print('The following files where created:')
        print('- ' + new_filename_simple_solution)
        print('- ' + new_filename_solution)
        print('- ' + new_filename_statistics)
        print('- ' + new_filename_velocities)
        print('- ' + new_filename_path_length_air_line_length)


if __name__ == '__main__':
    import doctest
    # running 'doctest.testmod()' test the examples in docstrings. Alternatively, the examples in docstrings can be tested by commenting out 'main()' and then navigating to the 'docs' directory and running the following command in the terminal (this will also show the test results for passed tests): python gps2net.py -v
    testResults = doctest.testmod()
    print(testResults)

    # run the main method
    main()
