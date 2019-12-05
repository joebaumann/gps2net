**************
gps2net script
**************


In the following, a documentation for all functions of the gps2net script is provided. Further, code examples illustrate thir usages.

Usage
=========

**gps2net** is a python module which addresses the problem outlined above. It contains a suitable algorithm to extract path data from GPS trajectories considering the underlying network topology.


Prerequisits
------------

In order to run the script, the filepaths of the txt files containing GPS data have to be specified in the main() function. In addition to that, the path of the shp-file (which provides the underlying network structure) has to be set.

Running the script
------------------

After the filepaths are specified in the main() function, the script can be run::

    python gps2net.py


Solution Output
--------------------

The output is saved in the *output_files* directory.
For every line of the input file (which contains the GPS positions) a dictionary (named after the input file) is created in *output_files*. Every created directory contains the following solution files:

- calculatedSolution.txt :: A txt which (in addition to the input GPS positions) contains additional parameters such as::
    - mapped position based on underlying network topology
    - path to next GPS position based on underlying network topology
    - path length
    - path time
    - average velocity
    - ...
- statistics.txt : A txt file which lists various statistics.
- timedifferencesPLOT.png : A figure which plots a histogram of all time differences between GPS measurements.
- velocitiesPLOT.png : A figure which plots a histogram of the velocities at all GPS positions.
- path_length_air_line_length_PLOT.png : A figure which plots a histogram of the length of the most likely path between two positions divided by the air line length between those two points. The histogram plots an overview of this measurement for at all GPS positions.


Running doctests
----------------

In order to run the doctests navigate to the *docs* directory and run the following command in the terminal::
    
    python gps2net.py -v


Functions
=========

.. automodule:: gps2net
    :members: