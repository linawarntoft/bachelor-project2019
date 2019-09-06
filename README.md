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
The code is written in c++ for runtime speed efficiency. For more detailed information about the code, please read the full thesis.
There are comments in the code to help you understand it as you go.

### StarCluster1
The numerical method is based on the one presented in Alexander and Gieles [2012] and reproduced for our project. For using the code:

1. Enter the relevant data in an `input.txt` in the format found in the folder *InputFiles*
2. Create a results folder named *Results*. Your data will be found here
3. Compile all the files in the *StarCluster1* directory using the `g++ *.cpp` command in your Linux terminal
4. Run the program using `./a.out` in the terminal

Remember to compile after every change to `input.txt`.

### StarCluster2
The numerical method is built upon the one used in *StarCluster1* and was based on Gieles et al. [2014]. For using the code:

Follow the same steps as presented for *StarCluster1*. Note that there are more parameters needed in `input.txt` for this version.

### StarCluster3
The numerical method used for *StarCluster3* is the same as *StarCluster2*. The difference here is that we are now considering potentially escaping stars as different bodies in our simulation. They are created using the object `PotEsc`. For using this code:

Follow the same steps as presented for *StarCluster2*. Note that there are more parameters needed in `input.txt` for this version. You may also choose to not enter any initial coordinates and velocities in the input-file.


#### Tides
The files contained within the folder *Tides* are written by my supervisor Florent Renaud. For more detailes on the theory behind his code, read Renaud et al. [2011] and Renaud and Gieles [2015]. For a shorter description, read part 2.3 of my thesis.

## Contact

For more information or questions about the code, please contact me at lwarntoft@gmail.com.

