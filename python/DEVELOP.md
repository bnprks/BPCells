## Building from source

1. Install required C++ dependencies, setting any required environment variables for
   compilation:
    - Eigen3:
        - Typically available via package managers, e.g. `apt-get install libeigen3-dev`
        - Installation path must be passed via standard environment variables, e.g. 
        `export CPATH="${CPATH:+$CPATH:}/usr/include/eigen3"` (note that typical installations
        prefix with the eigen3 directory, necessitating this CPATH modification)
    - HDF5:
        - Typically available via package managers, e.g. `apt-get install libhdf5-dev`
        - Any non-standard installation paths must be passed via the standard environment variables `CPATH` for include and `LIBRARY_PATH` for lib, 
        e.g. `export CPATH="${CPATH:+$CPATH:}/usr/include/hdf5/serial"` and `export LIBRARY_PATH=${LIBRARY_PATH:+$LIBRARY_PATH:}/usr/lib/hdf5/serial`
    - Highway:
        - Available [here](https://github.com/google/highway/), and also for Debian available via
        `apt-get install libhwy-dev`. Requires cmake to build from source
        - `bash ./scripts/install_highway_cmake.sh $(pwd)/highway` will fetch source and install, taking one optional
        command-line argument to pass installation prefix
        - Non-standard installation paths can be specified by setting `CPATH` and `LIBRARY_PATH` as with HDF5, e.g. `export CPATH=${CPATH:+$CPATH:}$(pwd)/highway/include; export LIBRARY_PATH=${LIBRARY_PATH:+$LIBRARY_PATH:}$(pwd)/highway/lib`
2. Run `pip install .`

On Linux/Mac machines the `scripts/install_{eigen,hdf5,highway}.sh` will perform builds from
source of the respective libraries. These scripts accept one optional parameter with the 
root directory to install into (default `/usr/local`). They also can optionally use the
`LIB_CACHE` environment variable to help cache outputs across CI jobs.

## Running tests
BPCells uses pytest, so just running `pytest` will run man unit tests.
Certain slow-running tests that require downloading data will only run if
the environment variable `BPCELLS_PYTEST_DATA` is set, pointing to the
directory that should be used to cache file downloads.

### Test dependencies

`pip install pytest h5py anndata`


### Documentation dependencies

`pip install sphinx myst_nb pydata_sphinx_theme jupytext`

## Example `.envrc`

One easy way to get set up with environment variables and a virtual environment for development is
with [`direnv`](https://direnv.net/). If you install direnv, you can make a `.envrc` in this folder
and it will automatically enable a virtualenv and set environment variables whenever you switch
into this folder.

The exact right setup will depend on your OS, but this works well for Ubuntu
```bash
layout python3

# Eigen3 installation path (depends on your OS setup; this works for Ubuntu/Debian)
export CPATH="${CPATH:+$CPATH:}/usr/include/eigen3"

# HDF5 installation path (depends on your OS setup; this works for Ubuntu/Debian)
export CPATH="${CPATH:+$CPATH:}/usr/include/hdf5/serial"
export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}/usr/lib/x86_64-linux-gnu/hdf5/serial"

# Export highway lib path for building
export CPATH="${CPATH:+$CPATH:}$(pwd)/highway/include"
export LIBRARY_PATH="${LIBRARY_PATH:+$LIBRARY_PATH:}$(pwd)/highway/lib"

# Use ccache to speed up re-compilation cycles
export CXX="ccache g++"
export CC="ccache gcc"

# Directory for real test data
export BPCELLS_PYTEST_DATA_CACHE="$(pwd)/tests/data"
```