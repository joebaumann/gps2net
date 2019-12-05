************
Introduction
************

Background
==========

Today, we live in an information society in which the creation, distribution, use, integration and manipulation of information is an important economic, political and cultural activity. New technologies have resulted in an information explosion and are profoundly changing all aspects of our lives. For example, satellite navigation systems such as GPS, GLONASS, Galileo, BeiDou revolutionized our perception of space. Global Positioning System (GPS) [#]_ permits the user to locate the position of any place anywhere in the world to the accuracy of fewer than 10 meters. Nowadays, GPS is an essential utility for a vast number of commercial, strategic and scientific applications ranging from air, satellite system (GIS), its applications increase manifold. It is utilised for the management of truck, taxi, ambulance, and other patrol fleets, in remote sensing and environment monioring, in disaster relief, in earth and atmospheric sciences and in several new applications which are emerging rapidly ([GeospatialWorld2010]_).


Goal
==========

A growing community of researchers in networks of scientific studies has recognised that through the combination of spatial-temporal-resolved sequential data (e.g., such as recorded GPS trajectories) and the underlying topology (e.g., a road network) higher-order, non-dyadic graph representations ([Rosvall2014]_; [Scholtes2014]_; [Scholtes2016]_, [Scholtes2017]_). Different from standard (first-order) network models, such higher-order models do not merely assume that paths in networks are transitive, but instead, explicitly model the (causal) pathways that result from the timing and chronological ordering of spatial-temporal resolved observations.

In order to perform such analyses, path data is required for the underlying network. Since GPS data only represents the position (latitudinal, longitudinal) at a certain point in time, they have to be converted into suitable trajectories. In addition, the accuracy of the recorded position varies, which makes it difficult to assign it to the edges in the network. This can occur especially in urban areas, where GPS signal is influenced by buildings and the density of edges is relatively high, which complicates the assignment process. Furthermore, it can happen that due to the high speeds or gaps in the position recording, no complete path can be observed, i.e., the position kumps from one edge to another. However, these edges are not topologically connected. In this respect, suitable algorithms are needed to converge GPS data into paths, taking into account the previously mentioned contraints.


Solution
========

**gps2net** is a python module which addresses the problem outlined above. It contains a suitable algorithm to extract path data from GPS trajectories considering the underlying network topology.


.. rubric:: Footnotes

.. [#] In the remainder of this text, the colloquial term GPS is used to describe the global navigation satellite system.