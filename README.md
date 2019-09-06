# Bachelor Project 2019
This is the code used for my bachelor thesis *Fast Simulations of Star Cluster Evolution*, 2019

Author: Lina Warntoft 

Supervisor: Florent Renaud

Department of Astronomy and Theoretical Physics, Lund University

(link: https://lup.lub.lu.se/student-papers/search/publication/8971716)

## Abstract
The evolution of globular star clusters tells the history of the formation of galaxies, how stellar streams might have formed and
how the future of the Universe might look like. Although the general theoretical background of the evolution is well understood, it
is difficult to conduct accurate and realistic computer simulations of objects containing thousands of interacting bodies. This
project introduces a new method tackling this issue, treating star clusters as point objects while they are losing stars as a result
of tidal stripping. This is all conducted based on the concept of potential escapers. Combining the computational effectiveness of a
simplified approximation of a star cluster with the accuracy of an N-body simulation, this project introduces a new option for
globular star cluster simulation.

The code written for this model allows to introduce any star cluster into any external potential. The star cluster is then described
using five main coupled differential equations: number of stars, half mass radius, total energy, core radius and density profile,
which, when integrated using a 4th order Runge-Kutta scheme, allows to model the full evolution of the cluster in less than a minute
on a single computer. Tracer particles for potential escapers are implemented so that the tidal stripping works in constant as well
as varying tidal fields. The results of this project show that the star loss over time is comparable to the ones conducted using much
slower direct N-body simulations, e.g. the step-like loss of stars in elliptical tidal fields, while at the same time having a run-time
of about a few minutes. Further development of this method might contribute to the understanding of stellar streams and build-up of
stellar populations in galaxies.

## The code and how to use it

### StarCluster1

### StarCluster2

### StarCluster3

#### Tides
The files contained within the folder *Tides* are written by my supervisor Florent Renaud.
