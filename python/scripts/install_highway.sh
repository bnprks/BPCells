set -euo pipefail
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

rm -rf highway-1.1.0