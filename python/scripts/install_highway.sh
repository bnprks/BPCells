set -euo
set -x


if [ "$#" -eq 1 ]; then
    if [ -f $1/lib/libhwy.a ]; then
        # Quit early if the output directory already exists (assume it's all good)
        exit 0
    fi
    INSTALL_DIR="$1"
else
    # cmake on linux defaults to /usr/local https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html
    INSTALL_DIR="/usr/local"
fi

# Check for cached files
if [ ! -z "${LIB_CACHE+x}" ] && [ -f "${LIB_CACHE}/hwy/lib64/libhwy.a" ]; then
    echo "Copying hwy from cache"
    mkdir -p "$INSTALL_DIR"/include "$INSTALL_DIR"/lib
    cp -r "${LIB_CACHE}"/hwy/include/* "$INSTALL_DIR"/include
    cp "${LIB_CACHE}"/hwy/lib*/libhwy.a "$INSTALL_DIR"/lib 
    exit 0;
fi

# Cache is empty, so fill cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    MAIN_INSTALL_DIR="$INSTALL_DIR"
    INSTALL_DIR="$LIB_CACHE/hwy"
    mkdir -p "$INSTALL_DIR"
fi

curl -L https://github.com/google/highway/archive/refs/tags/1.1.0.tar.gz > highway.tar.gz
tar -xzf highway.tar.gz
rm highway.tar.gz

cd highway-1.1.0
mkdir -p build-dir
cd build-dir

cmake .. \
    -DHWY_ENABLE_TESTS:BOOL=OFF \
    -DHWY_ENABLE_EXAMPLES:BOOL=OFF \
    -DHWY_ENABLE_CONTRIB:BOOL=OFF \
    -DCMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR"

make -j4 install

cd ../../

# Manually copy the math-inl contrib folder. (This is much faster than having to build all the vectorized sort algorithms)
mkdir -p "$INSTALL_DIR/include/hwy/contrib/math"
cp -r highway-1.1.0/hwy/contrib/math/math-inl.h "$INSTALL_DIR/include/hwy/contrib/math"

# Fill cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    mkdir -p "$MAIN_INSTALL_DIR"/include "$MAIN_INSTALL_DIR"/lib
    cp -r "${LIB_CACHE}"/hwy/include/* "$MAIN_INSTALL_DIR"/include

    # Can start in lib or lib64
    cp "${LIB_CACHE}"/hwy/lib*/libhwy.a "$MAIN_INSTALL_DIR"/lib 
fi

rm -rf highway-1.1.0