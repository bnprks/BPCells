set -euo
set -x


if [ "$#" -eq 1 ]; then
    if [ -f $1/lib/libhwy.a ]; then
        # Quit early if the output directory already exists (assume it's all good)
        exit 0
    fi
    SET_INSTALL_DIR="-DCMAKE_INSTALL_PREFIX:PATH=$1"
else
    SET_INSTALL_DIR=""
fi


wget https://github.com/google/highway/archive/refs/tags/1.0.7.tar.gz
tar -xzf 1.0.7.tar.gz
rm 1.0.7.tar.gz

cd highway-1.0.7
mkdir -p build-dir
cd build-dir

cmake .. \
    -DHWY_ENABLE_TESTS:BOOL=OFF \
    -DHWY_ENABLE_EXAMPLES:BOOL=OFF \
    $SET_INSTALL_DIR

make -j4 install

cd ../../
rm -r highway-1.0.7