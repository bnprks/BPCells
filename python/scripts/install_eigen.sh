set -euo
set -x


if [ "$#" -eq 1 ]; then
    INSTALL_DIR="$1"
else
    INSTALL_DIR="/usr/local"
fi

# Check for cached files
if [ ! -z "${LIB_CACHE+x}" ] && [ -f "${LIB_CACHE}/eigen/include/Eigen/Core" ]; then
    echo "Copying Eigen from cache"
    mkdir -p "$INSTALL_DIR/include"
    cp -r "${LIB_CACHE}/eigen/include/Eigen" "$INSTALL_DIR/include"
    exit 0;
fi

curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz > eigen.tar.gz
tar -xzf eigen.tar.gz
rm eigen.tar.gz

mkdir -p "$INSTALL_DIR"/include
cp -r eigen-3.4.0/Eigen $INSTALL_DIR/include

# Put files in cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    mkdir -p "${LIB_CACHE}/eigen/include"
    cp -r "$INSTALL_DIR/include/Eigen" "${LIB_CACHE}/eigen/include/"
fi

rm -rf eigen-3.4.0