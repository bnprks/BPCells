set -euo
set -x


if [ "$#" -eq 1 ]; then
    INSTALL_DIR="$1"
else
    INSTALL_DIR="/usr/local"
fi

# Check for cached files
if [ ! -z "${LIB_CACHE+x}" ] && [ -f "${LIB_CACHE}/hdf5/lib/libhdf5.so" ]; then
    echo "Copying hdf5 from cache"
    mkdir -p "$INSTALL_DIR"/include "$INSTALL_DIR"/lib 
    cp "${LIB_CACHE}"/hdf5/include/* "$INSTALL_DIR"/include
    cp "${LIB_CACHE}"/hdf5/lib/* "$INSTALL_DIR"/lib 
    exit 0;
fi

# Cache is empty, so fill cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    MAIN_INSTALL_DIR="$INSTALL_DIR"
    INSTALL_DIR="$LIB_CACHE/hdf5"
    mkdir -p "$INSTALL_DIR"
fi

curl -L https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.4.3.tar.gz > hdf5.tar.gz
tar -xzf hdf5.tar.gz
rm hdf5.tar.gz

cd hdf5-hdf5_1.14.4.3 
./configure --prefix="$INSTALL_DIR" --enable-tools=no --with-szlib=no --enable-tests=no
make -j4
make install
cd ..

# Fill cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    mkdir -p "$MAIN_INSTALL_DIR"/include "$MAIN_INSTALL_DIR"/lib
    cp "${LIB_CACHE}"/hdf5/include/* "$MAIN_INSTALL_DIR"/include
    cp "${LIB_CACHE}"/hdf5/lib/* "$MAIN_INSTALL_DIR"/lib 
fi

rm -rf hdf5-hdf5_1.14.4.3 