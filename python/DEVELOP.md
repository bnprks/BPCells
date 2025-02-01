## Building from source

The core dependencies required are: Eigen3, HDF5, and [Highway](https://github.com/google/highway/).

### Setup on Linux/Mac

Build dependencies from source (slower but always consistent):
```shell
mkdir -p build/deps-root
python scripts/install_deps.py build/deps-root
```

Then install via pip with appropriate environment variables set:
```shell
export CPATH="$(pwd)/build/deps-root/include"
export LIBRARY_PATH="$(pwd)/build/deps-root/lib:$(pwd)/build/deps-root/lib64"

pip install .
```

After having built the dependencies once, just running the `export` and `pip` commands is sufficient to re-build.

Re-building can also be sped up by installing `ccache` (`brew install ccache`, `sudo apt install ccache`, etc.), then setting:
```shell
export CXX="ccache g++"
export CC="ccache gcc"
```

### Setup on Windows

Recommend using vcpkg to install dependencies. ([Installation instructions](https://learn.microsoft.com/en-us/vcpkg/get_started/get-started?pivots=shell-cmd#1---set-up-vcpkg))

```shell
vcpkg install hdf5 eigen3 highway[contrib] zlib 
```

Then set environment variables and run pip (assuming vcpkg root is at `C:\vcpkg`)

PowerShell
```powershell
$env:CPATH = 'C:\vcpkg\installed\x64-windows\include'
$env:LIBRARY_PATH = 'C:\vcpkg\installed\x64-windows\lib'

pip install .
```

Command Prompt
```shell
set CPATH=C:\vcpkg\installed\x64-windows\include
set LIBRARY_PATH=C:\vcpkg\installed\x64-windows\lib

pip install .
```

## Running tests
BPCells uses pytest, so just running `pytest` will run man unit tests.
Certain slow-running tests that require downloading data will only run if
the environment variable `BPCELLS_PYTEST_DATA_CACHE` is set, pointing to the
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

# Export path to the pre-built dependencies
export CPATH="$(pwd)/build/deps-root/include"
export LIBRARY_PATH="$(pwd)/build/deps-root/lib:$(pwd)/build/deps-root/lib64"

# Use ccache to speed up re-compilation cycles
export CXX="ccache g++"
export CC="ccache gcc"

# Directory for real test data
export BPCELLS_PYTEST_DATA_CACHE="$(pwd)/tests/data"
```

## CI Deploy setup

BPCells uses cibuildwheel to deploy directly to PyPI, via `.github/workflows/pypi.yml`. These workflows
pass specific `CIBW_*` environment variables for Mac, Windows, and Linux in order to achieve a build process
that works on the Github Actions runners.

The main gotchas to be aware of:
- The final PyPI publish step needs manual approval before it runs, so set a timer to come back to approve
  the deploy.
- Because BPCells uses `setuptools_scm` to detect version number from git tag, PyPI will reject changes that
  don't come from a nice tagged version (and end up with a version like `bpcells-0.1.dev1+gf9636f4`). An example
  failed run is [here](https://github.com/bnprks/BPCells/actions/runs/10544359618). You can
  double-check the active version number from the `sdist` job progress, which takes around 10 seconds to run
  and displays the full version number at the end of its "Build a source tarball" log.
- There is some caching used to try to avoid re-compiling dependencies repeatedly. If a dependency upgrade
  isn't coming through as expected, you may need to remove the caching step or bump the cache key name.
