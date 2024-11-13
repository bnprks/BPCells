set -euo pipefail
set -x


if [ "$#" -eq 1 ]; then
    INSTALL_DIR="$1"
else
    INSTALL_DIR="/usr/local"
fi

curl -L https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.4.3.tar.gz > hdf5.tar.gz
tar -xzf hdf5.tar.gz
rm hdf5.tar.gz

cd hdf5-hdf5_1.14.4.3 
./configure --prefix="$INSTALL_DIR" --enable-tools=no --with-szlib=no --enable-tests=no
make -j4
make install
cd ..

rm -rf hdf5-hdf5_1.14.4.3 