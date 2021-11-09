
# A framework of functions for PIA approximation and interpolation on B-Splines in MATLAB

This repository contains the code I implemented for my degree thesis in Computer Science for Management at the University of Bologna, which I completed in 2021. 

**Thesis title**: Progressive Iteration Approximation: a class of methods for spline interpolation and approximation  
**Supervised by**: prof. Giulio Casciola  
**Candidate**: Giacomo Orsi

The thesis (in Italian) is available on the file `tesi.pdf`

## Abstract
The aim of this thesis is to introduce some iterative-geometric interpolation and approximation methods for spline curves and to analyse their differences. In particular, a class of methods called "Progressive Iteration Approximation" is analysed. These methods allow the control points of interpolating or approximating curves to be determined iteratively, with the advantage of handling very large data sets efficiently. In addition, a summary is presented of some numerical results obtained in MATLAB environment following the study, implementation and experimentation of the algorithms analysed.


## Codice 
This folder contains the code for the implementation for B-Spline curves in MATLAB of the algorithms
- PIA (Progressive Iteration Approximation) (https://doi.org/10.1016/j.cad.2013.08.012)
- WPIA (Weighted Progressive Iteration Approximation) (https://doi.org/10.1016/j.cagd.2009.11.001)
- PIA Adaptive (https://doi.org/10.1016/j.cagd.2012.03.005)
- LSPIA (Least Squares Progressive Iteration Approximation) (https://doi.org/10.1016/j.cad.2013.08.012)
- LSPIA Progressive (https://doi.org/10.1016/j.cad.2013.08.012)


PIA and WPIA are both contained within the folder `pia_bspline`. 

Each folder contains the file `<algorithmName>_body.m` and  `<algorithmName>_example.m`

At the beginning of the file  `body` it is possible to set all the execution parameters of the environment. 

In the files `example` there are some example parameters to execute `body`. 

In the folder `svg` there are several svg images composed of only one  *path*, and therefore suitable to be approximated/interpolated with only one curve. 

In the folder `librerie` all the used libraries are listed. 

To use the scripts, it is suggested to open the entire folder in MATLAB and execute the script `load_libraries.m` that adds all the libraries and the svg pictures to the working path. 

In the folder `figure_tesi` there is a file `.txt` with the instructions to replicate the figures of the thesis (in Italian). 