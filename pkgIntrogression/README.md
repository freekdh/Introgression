# pkgIntrogression
Introgression simulations complementing published paper

## System requirements
pkgIntrogression can be installed on unix-flavored systems, and requires the following key elements:

* C++11
* Open MP support
* boost libaries
* shiny library
* plotly library

## Installation
The easiest way to install velocyto.R is using devtools::install_github() from R:
```
library(devtools)
install_github("freekdh/Introgression/pkgIntrogression")
```
You need to have boost (e.g. `sudo apt-get install libboost-dev`) and openmp libraries installed. You can see detailed installation commands in the dockers/debian9/Dockerfile. 

### Dockers
If you are having trouble installing the package on your system, you can build a docker instance that can be used on a wide range of systems and cloud environments. To install docker framework on your system see [installation instruction](https://github.com/wsargent/docker-cheat-sheet#installation). After installing the docker system, use the following commands to build a velocyto.R docker instance:
```bash
cd velocyto.R/dockers/debian9
docker build -t velocyto .
docker run --name velocyto -it velocyto
```
