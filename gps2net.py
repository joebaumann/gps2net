# %%
# imports

from shapely.geometry import LineString
from matplotlib import pyplot
from shapely.geometry import Point, LineString, MultiLineString, mapping, shape
from shapely.ops import cascaded_union
from shapely import affinity
import fiona
import networkx as nx
import math


import time
from timeit import default_timer as timer


import sys
import os

# Disable


def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore


def enablePrint():
    sys.stdout = sys.__stdout__


# %%

# global variable: empty Directed Graph

# TODO: uncomment lines (just done for testing reasons so that graph is not recreated every time)
DG = nx.DiGraph()
# print("GraphFromSHP nr of edges START: ", DG.number_of_edges())

# constants used to measure used time of specific code fragements
timeCreatingGraph = 0
timeCreatingGraph2 = 0
timeAStar = 0
time_aStar_path = 0
time_aStar_length = 0

astar_copyGraph = 0
astar_addTarget = 0
astar_addSource = 0
astar_numberOfEdges = 0
astar_getEdgeData = 0
astar_addEdge = 0

no_path_AStar = 0

aStar_write = 0
timeClosestPoint = 0


# %%

# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 *
                                                     (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)

    # enablePrint()

    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

    # blockPrint()


# print("done")

# %%

# returns distance in meters

# APPROXIMATION
# QGIS could calculate accurate distance

def distFrom(lng1, lat1, lng2, lat2):
    earthRadius = 6371000  # meters
    dLat = math.radians(lat2-lat1)
    dLng = math.radians(lng2-lng1)
    a = math.sin(dLat/2) * math.sin(dLat/2) + math.cos(math.radians(lat1)) * \
        math.cos(math.radians(lat2)) * \
        math.sin(dLng/2) * math.sin(dLng/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    dist = earthRadius * c

    return dist


# %%
# how to cut a line
def cut(line, distance, point):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            # print("point.x: ", point.x)
            # print("cp.x: ", cp.x)
            # print("point.y: ", point.y)
            # print("cp.y: ", cp.y)
            return [
                LineString(coords[:i] + [(point.x, point.y)]),
                LineString([(point.x, point.y)] + coords[i:])]

# %%


def createGraphFromSHPInput(filepath_shp):

    GraphFromSHP = nx.DiGraph()

    # blockPrint()

    global timeCreatingGraph
    global timeCreatingGraph2

    counter = 0

    startGraph2 = timer()
    print("Graph already exists")

    endGraph2 = timer()
    timeCreatingGraph2 += (endGraph2-startGraph2)

    startGraph = timer()

    nr_elements_in_SHP_file = 0
    with fiona.open(filepath_shp) as input_lines:
        for i in input_lines:
            nr_elements_in_SHP_file += 1
    print("nr_elements_in_SHP_file: ", nr_elements_in_SHP_file)

    # Initial call to print 0% progress ProgressBar
    suffix = '| SHP file items: {}/{}'.format(
        counter+1, nr_elements_in_SHP_file)
    printProgressBar(0, nr_elements_in_SHP_file,
                     prefix='Creating Graph:', suffix=suffix, length=50)

    # print("graph nr of edges: ", GraphFromSHP.number_of_edges())
    print("starting to build graph: ")

    with fiona.open(filepath_shp) as input_lines:

        for line in list(input_lines):
            counter_inner = 0

            previous_coord = (0, 0)
            id = line['id']
            oneway = line['properties']['oneway']

            # print("id: ", id)
            # print("oneway: ", oneway)

            for coord in line['geometry']['coordinates']:

                if (previous_coord != (0, 0)):
                    # print("previous_coord: ", previous_coord)
                    # print("coord: ", coord)

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
                    # print("revious_coord    ==    (0,0)   --------------")

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

    print("graph nr of edges AFTER CREATION: ", GraphFromSHP.number_of_edges())

    # blockPrint()

    endGraph = timer()
    timeCreatingGraph += (endGraph-startGraph)

    # print("GraphFromSHP: ", (GraphFromSHP.edges(data=True)))

    # enablePrint()

    # print("graph nr of edges: ", GraphFromSHP.number_of_edges())

    # print("timeCreatingGraph: ", timeCreatingGraph)

    # blockPrint()

    return GraphFromSHP


# %%

def getShortestPathAStar(source, target, source_line, target_line, source_line_oneway, target_line_oneway, filepath_shp, ignore_oneway=False):

    global DG

    # blockPrint()

    global no_path_AStar
    global time_aStar_path
    global time_aStar_length

    global astar_copyGraph
    global astar_addTarget
    global astar_addSource
    global astar_numberOfEdges
    global astar_getEdgeData
    global astar_addEdge

    path = None
    path_length = None
    path_IDs = []
    all_added_edges = []

    # heuristic function for A star: returns the air line distance from source to the target

    def air_line_distance(source, target):
        distance = distFrom(source[0], source[1], target[0], target[1])
        return distance

    # function to temporarily add an edge to the graph. Only adds the edge if it deosn't exist yet. Further, new edge is added to list so that it can be removed again in the end
    def temporarily_add_edge_to_graph(startNode, endNode, edgeWeight, edgeId, direction):
        if (not DG.has_edge(startNode, endNode)):
            # print("adding ({},{})".format(startNode, endNode))

            DG.add_edge(startNode, endNode, weight=edgeWeight,
                        id=edgeId, oneway=direction)
            all_added_edges.append((startNode, endNode))
        # else:
            # print("graph does have({},{})".format(startNode,endNode))
            # print("edge:",DG.get_edge_data(startNode, endNode))

    astar_copyGraph_start = timer()

    # nrOfEdges=DG.number_of_edges()

    # if(nrOfEdges==0):
    if(nx.is_empty(DG)):

        DG = (createGraphFromSHPInput(filepath_shp))

    # DG=(createGraphFromSHPInput(filepath_shp)).copy()
    # print(DG.number_of_edges())

    astar_copyGraph_end = timer()
    astar_copyGraph += (astar_copyGraph_end-astar_copyGraph_start)

    astar_addTarget_start = timer()

    # check if target lies exactly on the beginning/end of a line segment
    # if yes, target already exists as a node
    if (target not in target_line):
        # remove the target_line edges and add the line_segments as edges

        d_target = LineString(target_line).project(Point(target))

        cut_line_target = cut(LineString(target_line), d_target, Point(target))
        # print("cut_line_target: ", cut_line_target)
        new_line_target_after_cut = [list(x.coords) for x in cut_line_target]

        len_line_segment = len(new_line_target_after_cut[0])
        edge_to_remove_T_source = new_line_target_after_cut[0][len_line_segment-2]
        edge_to_remove_T_target = new_line_target_after_cut[1][1]

        astar_getEdgeData_start = timer()

        # get the edge (including the attributes) that should be removed from the graph

        if(target_line_oneway == 'B'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_source, edge_to_remove_T_target)
            print("B:", edge_to_remove_from_graph)

        elif(target_line_oneway == 'F'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_source, edge_to_remove_T_target)
            print("F:", edge_to_remove_from_graph)

        elif(target_line_oneway == 'T'):
            edge_to_remove_from_graph = DG.get_edge_data(
                edge_to_remove_T_target, edge_to_remove_T_source)
            print("T:", edge_to_remove_from_graph)

        astar_getEdgeData_end = timer()
        astar_getEdgeData += (astar_getEdgeData_end-astar_getEdgeData_start)

        edge_to_remove_from_graph_id = edge_to_remove_from_graph['id']
        edge_to_remove_from_graph_oneway = edge_to_remove_from_graph['oneway']

        """

        # i think that it is not necessary to remove an edge

        # remove the edge from the graph
        if(remove_reverse==1):
            DG.remove_edge(edge_to_remove_source, edge_to_remove_target)
        else:
            DG.remove_edge(edge_to_remove_target, edge_to_remove_source)

        print("edge_to_remove_from_graph after removed: ",
              DG.get_edge_data(edge_to_remove_source,edge_to_remove_target))
        """

        astar_addEdge_start = timer()

        # add the line segments for the target
        print("source_line_oneway:", source_line_oneway)
        print("source_line_oneway:", source_line_oneway)
        print("check if ignore oneway is true:", ignore_oneway)

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

        astar_addEdge_end = timer()
        astar_addEdge += (astar_addEdge_end-astar_addEdge_start)

        astar_getEdgeData_start = timer()

        astar_getEdgeData_end = timer()
        astar_getEdgeData += (astar_getEdgeData_end-astar_getEdgeData_start)

    else:
        print("no need to add additional edges target")

    astar_addTarget_end = timer()
    astar_addTarget += (astar_addTarget_end-astar_addTarget_start)

    astar_addSource_start = timer()

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

        astar_getEdgeData_start = timer()

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

        astar_getEdgeData_end = timer()
        astar_getEdgeData += (astar_getEdgeData_end-astar_getEdgeData_start)

        edge_to_remove_from_graph_id = edge_to_remove_from_graph['id']
        edge_to_remove_from_graph_oneway = edge_to_remove_from_graph['oneway']

        astar_addEdge_start = timer()

        # if source and target lie on the same line segment, add a connection between the two
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

        astar_addEdge_end = timer()
        astar_addEdge += (astar_addEdge_end-astar_addEdge_start)

        astar_getEdgeData_start = timer()

        astar_getEdgeData_end = timer()
        astar_getEdgeData += (astar_getEdgeData_end-astar_getEdgeData_start)

    astar_addSource_end = timer()
    astar_addSource += (astar_addSource_end-astar_addSource_start)

    astar_numberOfEdges_start = timer()

    astar_numberOfEdges_end = timer()
    astar_numberOfEdges += (astar_numberOfEdges_end-astar_numberOfEdges_start)

    try:
        # try to find a path
        time_aStar_path_start = timer()

        path = nx.astar_path(DG, source, target,
                             heuristic=air_line_distance, weight='weight')

        time_aStar_path_end = timer()
        time_aStar_path += (time_aStar_path_end-time_aStar_path_start)

        time_aStar_length_start = timer()

        path_length = nx.astar_path_length(
            DG, source, target, heuristic=air_line_distance, weight='weight')

        # get the IDs of the traversed lines
        print('get the IDs of the traversed lines')
        prevNode = (0, 0)
        for node in path:
            print('prevNode:', prevNode)
            if (prevNode != (0, 0)):
                print('node', node)

                if(DG.get_edge_data(prevNode, node) != None):
                    print(DG.get_edge_data(prevNode, node)['id'])
                    if(DG.get_edge_data(prevNode, node)['id'] in path_IDs):
                        #sys.exit('id already exists')
                        print('id already exists')
                        print("path_IDs1:", path_IDs)
                    else:
                        path_IDs.append(DG.get_edge_data(prevNode, node)['id'])
                else:
                    path_IDs.append('edge not found! -- prevNode:', prevNode)
                    path_IDs.append('edge not found! -- node:', node)
                print("")
            prevNode = node
        print("path_IDs2:", path_IDs)

        time_aStar_length_end = timer()
        time_aStar_length += (time_aStar_length_end-time_aStar_length_start)

    except:
        no_path_AStar += 1
        """
        print("neighbors: ", list(DG.neighbors(source)))
        print("predecessors: ", list(DG.predecessors(source)))
        print("successors: ", list(DG.successors(source)))
        print("data: ", list(DG[source]))
        print("")
        print("---")
        print("")
        print("neighbors: ", list(DG.neighbors(target)))
        print("predecessors: ", list(DG.predecessors(target)))
        print("successors: ", list(DG.successors(target)))
        print("data: ", list(DG[target]))
        """

    # remove all the edges which where added to the graph.
    print("GraphFromSHP nr of edges 1: ", DG.number_of_edges())
    DG.remove_edges_from(all_added_edges)
    print("GraphFromSHP nr of edges 1: ", DG.number_of_edges())

    return path, path_length, path_IDs


# %%


def calculateClosestPointAndShortestPath(filepath, filepath_shp, minNumberOfLines=2, aStar=1):

    global timeClosestPoint

    timeClosestPointStart = timer()

    # blockPrint()

    """



    """

    # header = ["latitude(y);longitude(x);hasPassenger;time;closest_intersection_x;closest_intersection_y;intersected_line_as_linestring;linestring_adjustment_visualization;path_time;path_as_linestring;path_length;air_line_length;path_length/air_line_length;velocity_m_s;A_path_as_linestring;A_path_length;A_air_line_length;A_path_length/air_line_length;A_velocity_m_s\n"]
    header = ['latitude(y);longitude(x);hasPassenger;time;closest_intersection_x;closest_intersection_y;intersected_line_as_linestring;linestring_adjustment_visualization;path_time']

    if(aStar == 1):
        header += ';path_as_linestring;path_length;air_line_length;path_length/air_line_length;velocity_m_s;pathIDs;solution_id;solution_index;path_from_target_to_source;comment'

    header += '\n'

    counter = 0

    lines_for_new_text_file = []

    previous_target = (0, 0)
    previous_source = (0, 0)
    previous_intersected_line = None
    previous_timestamp = None
    previous_intersected_line_oneway = None

    def getLocationResult(filepath_shp, x, y, passenger, timestamp, previous_point, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway):

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

        global timeAStar
        global time_aStar_path
        global time_aStar_length
        global aStar_write

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
            location_result['previous_intersected_line'] = 'previous_intersected_line'

        if(previous_intersected_line_oneway != None):
            location_result['previous_intersected_line_oneway'] = previous_intersected_line_oneway
        else:
            location_result['previous_intersected_line_oneway'] = 'previous_intersected_line_oneway'

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
        location_result['intersected_line'] = ''
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

        # if point is oulier
        if(areasize >= max_areasize):
            print("This point is an outlier.")
            comment += 'This point is an outlier.'

            # makes sure that next point doesn't calculate shortest path
            previous_point = (0, 0)

        # point is not an outlier
        else:

            solution_dict = {}

            for input_line in input_shapes:

                lineID = input_line[1]['id']

                newLS = LineString(input_line[1]['geometry']['coordinates'])
                closest_point_on_line = newLS.interpolate(
                    newLS.project(focal_pt))

                inter_dict_point = {}

                inter_dict_point["closest_point_on_line"] = closest_point_on_line
                inter_dict_point["input_line"] = newLS

                inter_dict_point["oneway"] = input_line[1]['properties']['oneway']

                # distance needs to be float so that it can be sorted appropriately afterwards
                inter_dict_point['distance'] = float(
                    focal_pt.distance(closest_point_on_line))

                solution_dict[lineID] = inter_dict_point

            # sort the nested dict by the 'distance' which lies  within all the inner dicts
            # first element in the sorted dict is the point with the shortest distance to the focal_pt
            solution_dict_sorted_by_distance = sorted(
                solution_dict.items(), key=lambda kv: (kv[1]['distance']))

            print("solution NEW:")
            # get the first element as it is the one with the smallest distance (which is the closest point)
            solution_id = solution_dict_sorted_by_distance[0][0]
            solution = solution_dict_sorted_by_distance[0][1]
            location_result['all_solutions_sorted'] = solution_dict_sorted_by_distance
            location_result['solution_id'] = solution_id
            location_result['solution'] = solution

            for s in solution:
                print(str(s))
                print(str(solution[s]))
            print("")
            print("end solution")
            print("")
            print("")
            print("")

            # sys.exit('id already exists')

            closest_intersection_x = solution['closest_point_on_line'].x
            closest_intersection_y = solution['closest_point_on_line'].y

            intersected_line = solution['input_line']
            intersected_line_oneway = solution['oneway']

            # when appending the solution we have to parse the floats into strings
            location_result['closest_intersection_x'] = closest_intersection_x
            location_result['closest_intersection_y'] = closest_intersection_y

            # append the underlying map line on which the final point lies
            location_result['intersected_line'] = intersected_line

            # append the line which visualizes the way from the start point to the new position of the closest intersection point
            linestring_adjustment_visualization = LineString(
                [(focal_pt.x, focal_pt.y), (closest_intersection_x, closest_intersection_y)])

            location_result['linestring_adjustment_visualization'] = linestring_adjustment_visualization

            # check if previous point is set. if yes: calculate the shortest path from current point to previous point
            if(previous_point == (0, 0)):
                # no previous point set, do not calculate shortest path.
                print("cannot compute shortest path as previous point is (0,0)")
                comment += 'Cannot compute shortest path as previous point is (0,0). '

            elif((closest_intersection_x, closest_intersection_y) == previous_point):
                # if the taxi did not move, path is not created
                print("source==target. Taxi did not move.")

                comment += 'source==target. Taxi did not move. '

            else:

                path_time = abs(previous_timestamp-timestamp)

                location_result['path_time'] = path_time

                if(aStar == 1):

                    aStarStart = timer()

                    path = None
                    path2 = None
                    # calculate the shortest path with A STAR algorithm
                    path, path_length, pathIDs = getShortestPathAStar((closest_intersection_x, closest_intersection_y), previous_point, list(
                        intersected_line.coords), list(previous_intersected_line.coords), intersected_line_oneway, previous_intersected_line_oneway, filepath_shp)

                    print("")
                    print("path:", path)
                    print("")

                    aStarEnd = timer()
                    timeAStar += (aStarEnd-aStarStart)

                    aStar_write_Start = timer()

                    # calculate air line distance between source and target
                    air_line_length = distFrom(
                        closest_intersection_x, closest_intersection_y, previous_point[0], previous_point[1])

                    """
                    print("air_line_length s to t:", air_line_length)
                    print("path_length:", path_length)

                    print("if statement")
                    print(air_line_length < 20)
                    print(intersected_line_oneway != 'B')
                    print(previous_intersected_line_oneway != 'B')
                    print(list(intersected_line.coords) == list(
                        previous_intersected_line.coords))
                    """

                    # if this is true, it is possible that a taxi has to ride all around the block because it is a oneway street and source and target lie on the same line
                    if(air_line_length < 20 and intersected_line_oneway != 'B' and previous_intersected_line_oneway != 'B' and (list(intersected_line.coords) == list(previous_intersected_line.coords))):
                        print("if statement inside")

                        # calculate where the source lies on the line
                        d_source = LineString(list(intersected_line.coords)).project(
                            Point(closest_intersection_x, closest_intersection_y))

                        # calculate where the target lies on the line
                        d_target = LineString(list(previous_intersected_line.coords)).project(
                            Point(previous_point))

                        print("d_source:", d_source)
                        print("d_target:", d_target)

                        # check if there was a glipse in the gps coordinates which resulted in a situation where the source lies after the target
                        # if it is a oneway line from the FROM-node to the TO-node and source lies AFTER the target   OR   if it is a oneway line from the TO-node to the FROM-node and source lies BEFORE the target
                        if((intersected_line_oneway == 'F' and d_source > d_target)or(intersected_line_oneway == 'T' and d_source < d_target)):
                            print("if statement inside SECOND!!!")

                            # calculate the shortest path with A STAR algorithm
                            # set ignore_oneway=True : this makes sure that the source_line and the target_line both are treated as lines where driving in both directions is allowed (so feature 'oneway' is ignored)
                            path2, path_length2, pathIDs2 = getShortestPathAStar((closest_intersection_x, closest_intersection_y), previous_point, list(intersected_line.coords), list(
                                previous_intersected_line.coords), intersected_line_oneway, previous_intersected_line_oneway, filepath_shp, ignore_oneway=True)

                            print("path_length:", path_length)
                            print("path_length2:", path_length2)
                            print("path:", path)
                            print("path2:", path2)
                            print("")
                            print("check if path2 is better:")
                            print(path == None and path2 != None)
                            print(".. or ..")
                            print(path != None and path2 != None and (
                                path_length > path_length2))

                        # override path/path_lenght/pathIDs to make sure to use the solution which ignores the 'oneway' property
                        if((path == None and path2 != None) or (path != None and path2 != None and (path_length > path_length2))):
                            print("OVERRIDING THE PATH")
                            path = path2
                            path_length = path_length2
                            pathIDs = pathIDs2
                            comment += ' the oneway-property was ignored'
                            path_from_target_to_source = 1

                    # if a path was found, append the solution to the lines_for_new_text_file
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
                        print("No path was found with A Star algorithm.")

                    aStar_write_End = timer()
                    aStar_write += (aStar_write_End-aStar_write_Start)

        location_result['path_from_target_to_source'] = path_from_target_to_source
        location_result['comment'] = comment

        return location_result, (closest_intersection_x, closest_intersection_y), previous_point, intersected_line, previous_intersected_line, timestamp, intersected_line_oneway, previous_intersected_line_oneway

    # calculate the number of lines of the file we are looking at
    # lines_in_textfile=sum(1 for line in open(filepath))
    lines_in_textfile = 0
    with open(filepath, "r") as f:
        for line in f:
            lines_in_textfile += 1
    print(lines_in_textfile)

    # Initial call to print 0% progress ProgressBar
    suffix = '| current file: {}/{} lines'.format(counter+1, lines_in_textfile)
    printProgressBar(0, lines_in_textfile, prefix='Progress:',
                     suffix=suffix, length=50)

    previous_location_result = {}

    with open(filepath, "r") as f:
        mylist = f.read().splitlines()

        # append an artificaial line at the end. This will make sure that also the path from the last to the second-last element will be calculated and possible improvements will be checked.
        mylist.append('artificialline')

        for lines in mylist:

            # the last line of the text file is artificial as it was added before. this line is just needed to make sure the second last line (which is actually the last actual line of the text file is calculated correctly)
            if(lines != 'artificialline'):

                values = lines.split(" ")
                # the longitude and latitude were read from the txt file an saved as a string, here we transform it to a float
                x = float(values[1])
                y = float(values[0])
                passenger = int(values[2])
                timestamp = int(values[3])

                # get the result for the current location
                current_location_result, source, target, intersected_line, target_intersected_line, timestamp, intersected_line_oneway, target_intersected_line_oneway = getLocationResult(filepath_shp, x, y, passenger, timestamp, previous_source, previous_intersected_line, previous_timestamp, previous_intersected_line_oneway)

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
                        print("all_solutions_sorted:",
                              previous_location_result['all_solutions_sorted'])
                        print("solution_id:",
                              previous_location_result['solution_id'])
                        print("solution:",
                              previous_location_result['solution'])

                        print('')
                        print('FOR LOOP:')

                        if(str(previous_location_result['solution_id']) in pathIDs_to_previous_target or previous_location_result['solution_id'] in pathIDs_to_previous_target):
                            print(
                                'Solution is already on the shortest path between the points before and after --> therefore the solution can not be improved.')

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
                                print("k:", k)
                                print("v:", v)
                                new_solution_index += 1

                            if (new_solution_id == (-1)):
                                print(
                                    'There is no other solution which lies on the shortest path.')
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

                                current_location_result_new, source_new, target_new, intersected_line_new, target_intersected_line_new, timestamp_new, intersected_line_oneway_new, target_intersected_line_oneway_new = getLocationResult(filepath_shp, x, y, passenger, timestamp, (new_solution_point_x, new_solution_point_y), new_solution['input_line'], previous_location_result['timestamp'], new_solution['oneway'])

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

                                previous_location_result_new, previous_source_new, previous_target_new, previous_intersected_line_new, previous_target_intersected_line_new, previous_timestamp_new, previous_intersected_line_oneway_new, previous_target_intersected_line_oneway_new = getLocationResult(filepath_shp, new_solution_point_x, new_solution_point_y, previous_location_result['passenger'], previous_location_result['timestamp'], previous_location_result['target'], previous_location_result['previous_intersected_line'], previous_location_result['previous_timestamp'], previous_location_result['previous_intersected_line_oneway'])
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
                                    else:
                                        new_comment = 'Checked the second solution (with solution_id: {} and solution_index: {}) but the new cummulated path_length ({}) is longer that the initial commulated path_length ({}) and solution_id: {}).'.format(
                                            new_solution_id, new_solution_index, new_commulated_path_lengt, previous_commulated_path_lengt, previous_location_result['solution_id'])
                                        print(new_comment)
                                        previous_location_result['comment'] += new_comment

                                else:
                                    # no path was found with new solution
                                    new_comment = 'No path found for new solution. '
                                    print(new_comment)
                                    previous_location_result['comment'] += new_comment

                        #sys.exit('WE STOP HERE')

                else:
                    print('Tried to check other solution but the target is not known. ')

                # append the result of the previous location to the entire solution
                lines_for_new_text_file.append(previous_location_result)

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

            # stop algorithm after a few lines (for testing reasons)
            counter += 1
            print(counter)
            if (counter > 10):
                print("ende")
                break

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

        """
        NOT NEEDED ANYMORE, as we have artificial line now

        # the result of a location is only added to the intire solution after it is confirmed as the best solution by the next location
        # for this reason the result of the last location of the text_file has to be added to the entire solution (since there isn't a next solution that will do that)
        lines_for_new_text_file.append(previous_location_result)
        """

    timeClosestPointEnd = timer()

    timeClosestPoint += (timeClosestPointEnd-timeClosestPointStart)

    return header, lines_for_new_text_file


# %%
# only selected data points from this taxi
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/new_abboip_selection_test.txt'
# all data points from this taxi
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_selection_for_SP_tests.txt'

# blockPrint()
startTotal = timer()

#filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy.txt'

# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small.txt'
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall.txt'
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_2.txt'
# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_3.txt'

# filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_small_verysmall_4.txt'

#filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_verysmall_closestIsBest_oneline.txt'

#filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_SMALL_to_check_outlier.txt'

##filepath = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/new_abboip_copy_verysmall_closestIsBest_1stSolution.txt'


filepath = '/Users/Joechi/Google Drive/gps2net/Data/test_data/just_one_taxi/new_abboip_copy.txt'



#filepath_shp = '/Users/Joechi/Google Drive/HS19 – PathPy/2_Taxi data/Tests/Exports/AllPointsForOneTaxi/geo_SF_lines_exported_for_testing.shp'

filepath_shp = '/Users/Joechi/Google Drive/gps2net/Data/taxi_san_francisco/San Francisco Basemap Street Centerlines/geo_export_e5dd0539-2344-4e87-b198-d50274be8e1d.shp'

filepathIndex = filepath.rfind('/')
filepathIndex2 = filepath.rfind('.')
new_filename = filepath[filepathIndex+1:filepathIndex2]
new_filename += '_testfile.txt'


myHeader, myCalculatedSolution = calculateClosestPointAndShortestPath(
    filepath, filepath_shp, minNumberOfLines=2, aStar=1)

print('')
print('')
print('')
print('')
print('')
print('')
print('')
#print('NEW returned calculated solution: ')
# print(myCalculatedSolution)
print('')
print('')
print('')
print('')
print('')
print('')


# print(myHeader)
# print("data:")
# print(myCalculatedSolution)
# print("data end.")


# this saves a new text file which includes the calculated parameters

# with open("CreatedFiles3/ImprovedAlgorithm/"+new_filename, "w") as new_file:
#with open("2_Taxi data/CreatedFiles3/ImprovedAlgorithm2/"+new_filename, "w") as new_file:
with open("output_files/"+new_filename, "w") as new_file:
    new_file.writelines(myHeader)

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
        new_file.write(str(location_result['path_length/air_line_length']))
        new_file.write(';')
        new_file.write(str(location_result['velocity_m_s']))
        new_file.write(';')
        new_file.write(str(location_result['pathIDs']))
        new_file.write(';')
        new_file.write(str(location_result['solution_id']))
        new_file.write(';')
        new_file.write(str(location_result['solution_index']))
        new_file.write(';')
        new_file.write(str(location_result['path_from_target_to_source']))
        new_file.write(';')
        new_file.write(str(location_result['comment']))
        new_file.write('\n')
        # new_file.writelines(myCalculatedSolution)

endTotal = timer()


# enablePrint()

print("no_path_AStar :", no_path_AStar)
print("---")
print("timeTotal :", endTotal-startTotal)
print("timeCreatingGraph :", timeCreatingGraph)
print("timeCreatingGraph2 :", timeCreatingGraph2)
print("timeAStar :", timeAStar)
print("time_aStar_path :", time_aStar_path)
print("time_aStar_length :", time_aStar_length)

print("")
print("XXXXXXXX")
print("")

print("astar_copyGraph: ", astar_copyGraph)
print("astar_addTarget: ", astar_addTarget)
print("astar_addSource: ", astar_addSource)
print("astar_numberOfEdges: ", astar_numberOfEdges)
print("astar_getEdgeData: ", astar_getEdgeData)
print("astar_addEdge: ", astar_addEdge)

print("")

print("-------")
print("")


print("aStar_write :", aStar_write)
print("timeClosestPoint Total :", timeClosestPoint)
print("timeClosestPoint without graphs :",
      timeClosestPoint-timeAStar)
print("")
print("-------")
print("")


print("GraphFromSHP nr of edges END: ", DG.number_of_edges())


# %%
