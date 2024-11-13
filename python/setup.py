import glob
import os
import os.path
import platform

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, ParallelCompile

from setuptools import setup
from setuptools.command.build import build
import setuptools_scm


# Handle improperly checked-out symlinks from git by manual copying
# This happens for certain Windows users depending on git config
import shutil
for path in ["bpcells-cpp", "vendor"]:
    if not os.path.islink(f"src/{path}") and not os.path.isdir(f"src/{path}"):
        contents = open(f"src/{path}", "rb").read()
        os.remove(f"src/{path}")
        shutil.copytree(f"../r/src/{path}", f"src/{path}")
        open(f"src/{path}/old-{path}", "wb").write(contents)

class BuildHook(build):
    def run(self): 
        try:
            build.run(self)
        finally:
            # Undo the manual copying of symlinks
            for path in ["bpcells-cpp", "vendor"]:
                if not os.path.islink(f"src/{path}") and os.path.isfile(f"src/{path}/old-{path}"):
                    contents = open(f"src/{path}/old-{path}", "rb").read()
                    shutil.rmtree(f"src/{path}")
                    open(f"src/{path}", "wb").write(contents)

extra_compile_args = []
if "BPCELLS_DEBUG" in os.environ:
    extra_compile_args.append("-g")

if "BPCELLS_NO_OPT" in os.environ:
    extra_compile_args.append("-O0")
else:
    extra_compile_args.append("-O2")

# Plan: take an environment variable for the libs
# - HWY_INCLUDE_DIR: Include path for hwy library
# - HWY_LIB_DIR: Lib path for hwy library
libraries = ["hdf5", "hwy"]
include_dirs = ["src/bpcells-cpp", "src/vendor"]
library_dirs = []

if platform.system() == "Windows":
    extra_compile_args = ["/std:c++17"]
    if "BPCELLS_DEBUG" in os.environ:
        extra_compile_args.append("/Zi")

    if "BPCELLS_NO_OPT" in os.environ:
        extra_compile_args.append("/Od")
    else:
        extra_compile_args.append("/O2")

    libraries.append("zlib")

    # Mimic the unix compiler environment variable settings
    if "CPATH" in os.environ:
        include_dirs.extend(os.environ["CPATH"].split(";"))
    if "LIBRARY_PATH" in os.environ:
        library_dirs.extend(os.environ["LIBRARY_PATH"].split(";"))
else:
    extra_compile_args = [
        "-std=c++17",
        '-Wno-unused-but-set-variable',
    ] 
    if "BPCELLS_DEBUG" in os.environ:
        extra_compile_args.append("-g")

    if "BPCELLS_NO_OPT" in os.environ:
        extra_compile_args.append("-O0")
    else:
        extra_compile_args.append("-O2")

    libraries.append("z")

    # Mimic the unix compiler environment variable settings
    if "CPATH" in os.environ:
        include_dirs.extend(os.environ["CPATH"].split(":"))
    if "LIBRARY_PATH" in os.environ:
        library_dirs.extend(os.environ["LIBRARY_PATH"].split(":"))

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
        extra_compile_args = extra_compile_args,
        libraries=libraries,
        include_dirs = include_dirs,
        library_dirs = library_dirs,
        ),
]

ParallelCompile("BPCELLS_NUM_BUILD_JOBS").install()

setup(
    cmdclass={"build":BuildHook},
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
)

