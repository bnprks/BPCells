set -euo
set -x


if [ "$#" -eq 1 ]; then
    if [ -f $1/lib/libhwy.a ]; then
        # Quit early if the output directory already exists (assume it's all good)
        exit 0
    fi
    SET_INSTALL_DIR="-DCMAKE_INSTALL_PREFIX:PATH=$1"
    INSTALL_DIR="$1"
else
    SET_INSTALL_DIR=""
    # cmake on linux defaults to /usr/local https://cmake.org/cmake/help/latest/variable/CMAKE_INSTALL_PREFIX.html
    INSTALL_DIR="/usr/local"
fi


wget https://github.com/google/highway/archive/refs/tags/1.1.0.tar.gz
tar -xzf 1.1.0.tar.gz
rm 1.1.0.tar.gz

cd highway-1.1.0
mkdir -p build-dir
cd build-dir

cmake .. \
    -DHWY_ENABLE_TESTS:BOOL=OFF \
    -DHWY_ENABLE_EXAMPLES:BOOL=OFF \
    -DHWY_ENABLE_CONTRIB:BOOL=OFF \
    $SET_INSTALL_DIR

make -j4 install

cd ../../

# Manually copy the math-inl contrib folder. (This is much faster than having to build all the vectorized sort algorithms)
mkdir -p "$INSTALL_DIR/include/hwy/contrib/math"
cp -r highway-1.1.0/hwy/contrib/math/math-inl.h "$INSTALL_DIR/include/hwy/contrib/math"

rm -r highway-1.1.0