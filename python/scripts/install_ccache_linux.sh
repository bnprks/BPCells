set -euo
set -x


if [ "$#" -eq 1 ]; then
    INSTALL_DIR="$1"
else
    INSTALL_DIR="/usr/local"
fi

# Check for cached files
if [ ! -z "${LIB_CACHE+x}" ] && [ -f "${LIB_CACHE}/ccache" ]; then
    echo "Copying ccache from cache"
    mkdir -p "$INSTALL_DIR/bin"
    cp -r "${LIB_CACHE}/ccache" "$INSTALL_DIR/bin"
    exit 0;
fi

curl -L https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-linux-x86_64.tar.xz > ccache.tar.xz
tar -xJf ccache.tar.xz
rm ccache.tar.xz

mkdir -p $INSTALL_DIR/bin
cp ccache-*/ccache $INSTALL_DIR/bin

# Put files in cache
if [ ! -z "${LIB_CACHE+x}" ]; then
    mkdir -p "${LIB_CACHE}"
    cp "$INSTALL_DIR/bin/ccache" "${LIB_CACHE}/ccache"
fi


rm -rf ccache-*