## BPCells python bindings

The current BPCells python bindings provide limited access to BPCells functionality. This
mainly focuses on integer matrix slicing / import, and pseudobulk fragment insertion counts.
The API is subject to change somewhat once a full-featured port of the BPCells R functionality is
completed.

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
        - Any non-standard installation paths must be passed via the standard environment variables, 
        e.g. `CPATH` and `LIB_PATH`
    - Highway:
        - Available [here](https://github.com/google/highway/), and also for Debian available via
        `apt-get install libhwy-dev`. Requires cmake to build from source
        - `bash ./scripts/install_highway_cmake.sh $(pwd)/highway` will fetch source and install, taking one optional
        command-line argument to pass installation prefix
        - Non-standard installation paths can be specified by setting the environment variables
        `HWY_INCLUDE_DIR` and `HWY_LIB_DIR`, e.g. `export HWY_INCLUDE_DIR=$(pwd)/highway/include; export HWY_LIB_DIR=$(pwd)/highway/lib`
2. Run `pip install .`

## Running tests
BPCells uses pytest, so just running `pytest` will run man unit tests.
Certain slow-running tests that require downloading data will only run if
the environment variable `BPCELLS_PYTEST_DATA` is set, pointing to the
directory that should be used to cache file downloads.