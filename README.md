This repository is an incomplete Fortran version of the Python [Constrained Spectral Approximation Method](https://github.com/ray-chew/spec_appx).

The Fortran implementation was meant to be the orographic gravity-wave source for the [Multi-scale Gravity-wave Model](https://arxiv.org/abs/2309.11257) that is coupled to the German Weather Service's [Icosahedral Nonhydrostatic (ICON) model](https://icon-model.org/). This codebase builds, runs, and outputs the spectral approximation of topographic data distributed in a geodesic grid cell. However, not all of the latest algorithmic and scientific features have been implemented in the Fortran code. 

  * A simple plotter for the output data can be found in `notebooks/plotter.ipynb`.
  * Namelist files are in `scripts`.

Note that this codebase uses features from the [Fortran 2008 standard](https://en.wikipedia.org/wiki/Fortran#Fortran_2008).

## Requirements
* gfortran-11
* BLAS
* LAPACK
* [fortran-stdlib](https://github.com/fortran-lang/stdlib)
* OpenMP
* [NetCDF-Fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/)
* CMake

Refer to `CMakeLists.txt` for more details.

## Build
```bat
cd trifourtopo
mkdir build
cd build
cmake ..
make
```

The `make` command will create two binaries in `build/bin`:
* `linker`
* `four_topo`

## Run
This codebase uses the [USGS GMTED 2010 topographic dataset](https://www.usgs.gov/coastal-changes-and-impacts/gmted2010) and an [ICON grid](http://icon-downloads.mpimet.mpg.de/).

Before executing `four_topo` for the spectral analysis, you **must** link the topographic dataset to the corresponding grid by setting the path in `linker.nml` and executing `build/bin/linker`. The program extends the NetCDF grid data with a lookup table to the corresponding data slice in the topographic dataset.

The binary `four_topo` outputs a NetCDF4 file comprising the computed spectral approximation over all the cells provided in the grid file. The output may be visualised with `notebooks/plotter.ipynb`.

The notebook `nc_compactifie.ipynb` creates a regional grid and topographic dataset and may help with debugging and testing.
