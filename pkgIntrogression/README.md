# pkgIntrogression
Introgression simulations complementing paper in prep

## System requirements
pkgIntrogression can be installed on unix-flavored systems, and requires the following key elements:

* C++11
* Open MP support
* boost libaries
* shiny package
* plotly package

## Installation
The easiest way to install pkgIntrogression is using devtools::install_github() from R:
```
library(devtools)
install_github("freekdh/Introgression/pkgIntrogression")
```
You need to have boost (e.g. `sudo apt-get install libboost-dev`) and openmp libraries installed.
