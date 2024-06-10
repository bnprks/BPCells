import glob
import os
import os.path

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, ParallelCompile

from setuptools import setup
from setuptools.command.install import install
import setuptools_scm


# Handle improperly checked-out symlinks from git by manual copying
# This happens for certain Windows users depending on git config
if not os.path.islink("src/bpcells-cpp") and not os.path.isdir("src/bpcells-cpp"):
    import shutil
    contents = open("src/bpcells-cpp", "rb").read()
    os.remove("src/bpcells-cpp")
    shutil.copytree("../r/src/bpcells-cpp", "src/bpcells-cpp")
    open("src/bpcells-cpp/old-bpcells-cpp", "wb").write(contents)

class InstallHook(install):
    def run(self):
        # Build HWY if needed
        
        install.run(self)

        # Undo the manual copying of symlinks
        if not os.path.islink("src/bpcells-cpp") and os.path.isfile("src/bpcells-cpp/old-bpcells-cpp"):
            import shutil
            contents = open("src/bpcells-cpp/old-bpcells-cpp", "rb").read()
            shutil.rmtree("src/bpcells-cpp")
            open("src/bpcells-cpp", "wb").write(contents)

extra_args = []
if "BPCELLS_DEBUG" in os.environ:
    extra_args.append("-g")

if "BPCELLS_NO_OPT" in os.environ:
    extra_args.append("-O0")
else:
    extra_args.append("-O2")

# Plan: take an environment variable for the libs
# - HWY_INCLUDE_DIR: Include path for hwy library
# - HWY_LIB_DIR: Lib path for hwy library
include_dirs = ["src/bpcells-cpp"]
library_dirs = []
if "HWY_INCLUDE_DIR" in os.environ:
    include_dirs.append(os.environ["HWY_INCLUDE_DIR"])
if "HWY_LIB_DIR" in os.environ:
    library_dirs.append(os.environ["HWY_LIB_DIR"])

def get_version():
    return setuptools_scm.get_version(root="..", relative_to=__file__)

ext_modules = [
    Pybind11Extension("bpcells.cpp",
        # Sort input files to ensure reproducible builds (https://github.com/pybind/python_example/pull/53)
        sorted(glob.glob("src/**/*.cpp", recursive=True)),
        # Example: passing in the version to the compiled code
        define_macros = [
            ('VERSION_INFO', get_version()), 
            ('EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS', None)
        ],
        extra_compile_args = [
            '-std=c++17',
            '-Wno-unused-but-set-variable',
        ] + extra_args,
        libraries=["hdf5", "hwy"],
        include_dirs = include_dirs,
        library_dirs = library_dirs,
        ),
]

ParallelCompile("BPCELLS_NUM_BUILD_JOBS").install()

setup(
    cmdclass={"install":InstallHook},
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
)

