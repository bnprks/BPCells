set -euo pipefail
set -x


if [ "$#" -eq 1 ]; then
    INSTALL_DIR="$1"
else
    INSTALL_DIR="/usr/local"
fi

curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz > eigen.tar.gz
tar -xzf eigen.tar.gz
rm eigen.tar.gz

mkdir -p "$INSTALL_DIR"/include
cp -r eigen-3.4.0/Eigen $INSTALL_DIR/include

rm -rf eigen-3.4.0