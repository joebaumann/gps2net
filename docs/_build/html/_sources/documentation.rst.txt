*****************
Algorithm
*****************


In the following, a documentation for all functions of the gps2net algorithm is provided. Further, code examples illustrate thir usages.

Usage
=========

**gps2net** is a python module which addresses the problem outlined above. It contains a suitable algorithm to extract path data from GPS trajectories considering the underlying network topology.


Prerequisits
------------

- In order to run the algorithm, the filepaths of the txt files containing GPS data have to be specified in the main() function.
- The GPS data in the txt files must be sorted (first line is newest).
- The path of the shp-file (which provides the underlying network structure) has to be set.

Running the Algorithm
---------------------

After the filepaths are specified in the main() function, the algorithm can be run::

    python gps2net.py


Solution Output
--------------------

The output is saved in the *output_files* directory.
For every line of the input file (which contains the values) a dictionary (named after the input file) is created in *output_files*. Every created directory contains the following solution files:

- **pathFromUnmappedGpsPositions.txt**
    A txt file which (in addition to the input values) contains additional parameters which are NOT based on the underlying street network, such as:
    
    - path : Linestring
        Linear path with start=current_position and end=next_position.
    - path length : float
        Distance between current and next position in meters.
    - path time : float
        Time between measurement at current and next position in seconds.
    - velocity : float
        Velocity in meters/second.
- **calculatedSolution.txt**
    A txt file which (in addition to the input values) contains additional parameters which are all based on the underlying street network, such as:

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
- **statistics.txt**
    A txt file which lists various statistics, such as e.g. the following numbers:

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
- **timedifferencesPLOT.png**
    A figure which plots a histogram of all time differences between GPS measurements.
- **velocitiesPLOT.png**
    A figure which plots a histogram of the velocities at all GPS positions.
- **path_length_air_line_length_PLOT.png**
    A figure which plots a histogram of the length of the most likely path between two positions divided by the air line length between those two points. The histogram plots an overview of this measurement for at all GPS positions.


Visualizing the Solution Output
-------------------------------

The solution can be visualized using `QGIS`_.

- In order to visualize Points (such as e.g. the measured or the mapped GPS positions), add it as a layer in GQIS (*Layer* > *Add Layer* > Add delimited text layer --> Choose *Point Coordinates* as geometry definition and 'longitude(x)' as X field and 'latitude(y)' as Y field respectively).
- In order to visualize Linestrings (such as e.g. path or linestring_adjustment_visualization), add it as a layer in GQIS (*Layer* > *Add Layer* > Add delimited text layer --> Choose *WKT* as geometry definition and 'linestring_adjustment_visualization' as geometry field).

.. _`QGIS`: https://www.qgis.org/


Running Doctests
----------------

In order to run the doctests make sure to first comment out the call of the main() function gps2net module::

    # main()

Then navigate to the *docs* directory and run the following command in the terminal::
    
    python gps2net.py -v


Code Documentation: Functions
=============================

.. automodule:: gps2net
    :members: