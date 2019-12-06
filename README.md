# gps2net: Extraction of path data from GPS trajectories considering the underlying network topology

<table>
      <tr><td>Date</td><td>January 2020</td></tr>
      <tr><td>Author</td><td>Joachim Baumann</td></tr>
      <tr><td>Contact</td><td>joachimbaumann1@gmail.com</td></tr>
      <tr><td>Documentation</td><td>https://joehonig.github.io/gps2net/</td></tr>
</table>

## Note

> If you have any problems, found bugs in the code or have comments or questions, please feel free to send an e-mail to [Joachim Baumann](mailto:joachimbaumann1@gmail.com).

## Background
A growing community of researchers in networks of scientific studies has recognised that through the combination of spatial-temporal-resolved sequential data (e.g., such as recorded GPS<sup id="a1">[1](#f1)</sup> trajectories) and the underlying topology (e.g., a road network) higher-order, non-dyadic graph representations (Rosvall *et al.*, 2014<sup id="a2">[2](#f2)</sup>; Scholtes, 2017<sup id="a3">[3](#f3)</sup>; Scholtes *et al.*, 2016<sup id="a4">[4](#f4)</sup>, 2014<sup id="a5">[5](#f5)</sup>). Different from standard (first-order) network models, such higher-order models do not merely assume that paths in networks are transitive, but instead, explicitly model the (causal) pathways that result from the timing and chronological ordering of spatial-temporal resolved observations.

In order to perform such analyses, path data is required for the underlying network. Since GPS data only represents the position (latitudinal, longitudinal) at a certain point in time, they have to be converted into suitable trajectories. In addition, the accuracy of the recorded position varies, which makes it difficult to assign it to the edges in the network. This can occur especially in urban areas, where GPS signal is influenced by buildings and the density of edges is relatively high, which complicates the assignment process. Furthermore, it can happen that due to the high speeds or gaps in the position recording, no complete path can be observed, i.e., the position jumps from one edge to another. However, these edges are not topologically connected. In this respect, suitable algorithms are needed to converge GPS data into paths, taking into account the previously mentioned constraints.

## Data

In order to validate the algorithm, the data-set of mobility traces of taxi cabs in San Francisco, USA,<sup id="a6">[6](#f6)</sup> in combination with the open-source street map of San Francisco<sup id="a7">[7](#f7)</sup> was used.

The following data was used:

**GPS data:**

Dataset of mobility traces of taxi cabs in San Francisco, USA.

- File: cabspottingdata.tar.gz
- URL: https://crawdad.org/epfl/mobility/20090224/

*Note: Registration and log in is necesary in order to get the data*

**GIS Road Map data:**

- File: San Francisco Basemap Street Centerlines.zip
- URL: https://data.sfgov.org/Geographic-Locations-and-Boundaries/San-Francisco-Basemap-Street-Centerlines/7hfy-8sz8/data


## Solution

>**gps2net** is a python module which addresses the problem outlined above. It contains a suitable algorithm to extract path data from GPS trajectories considering the underlying network topology.

###### Prerequisites

In order to run the script, the filepaths of the txt files containing GPS data have to be specified in the main() function. In addition to that, the path of the shp-file (which provides the underlying network structure) has to be set.

###### Running the script

After the filepaths are specified in the main() function, the script can be run:
```
python gps2net.py
```

###### Solution Output

The output is saved in the *output_files* directory.
For every line of the input file (which contains the GPS positions) a dictionary (named after the input file) is created in *output_files*. Every created directory contains the following solution files:
- calculatedSolution.txt : A txt which (in addition to the input GPS positions) contains additional parameters such as:
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


## Running doctests

In order to run the doctests navigate to the *docs* directory and run the following command in the terminal:
```
python gps2net.py -v
```

## Footnotes & References

<b id="f1">1</b> In the remainder of this text, the colloquial term GPS is used to describe the global navigation satellite system. [↩](#a1)

<b id="f2">2</b> Rosvall, M.; Esquivel, A. V.; Lancichinetti, A.; West, J. D.; Lambiotte, R. (2014). Memory in network flows and its effects on spreading dynamics and community detection. *Nature communications* **5**, 4630. [↩](#a2)

<b id="f3">3</b> Scholtes, I. (2017). When is a network a network?: Multi-order graphical model selection in pathways and temporal networks. In: *Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*. ACM, pp. 1037-1046. [↩](#a3)

<b id="f4">4</b> Scholtes, I.; Pfitzner, R.; Garas, A.; Tessone, C. J.; Schweitzer, F. (2014). Causality-driven slow-down and speed-up of diffusion in non-Markovian temporal networks. *Nature communications* **5**, 5024. [↩](#a4)

<b id="f5">5</b> Scholtes, I.; Wider, N.; Garas, A. (2016). Higher-order aggregate networks in the analysis of temporal networks: path structures and centralities. *The European Physical Journal B* **89(3)**, 61. [↩](#a5)

<b id="f6">6</b> Piorkowski, M.; Sarafijanovic‑Djukic, N.; Grossglauser, M. (2009). CRAWDAD dataset epfl/mobility (v. 2009‑02‑24), downloaded from https://crawdad.org/epfl/mobility/20090224, https://doi.org/10.15783/C7J010. [↩](#a6)

<b id="f7">7</b> https://datasf.org/ [↩](#a7)
