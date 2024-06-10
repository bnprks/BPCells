
INCLUDE_FILE=$1
TEMPDIR=$(mktemp -d -t includetest.XXXXX)

echo '
#include <'"$INCLUDE_FILE"'>
int main() {
    return 0;
}
' > $TEMPDIR/test.cc

$CXX $TEMPDIR/test.cc -o $TEMPDIR/test.out 2> /dev/null && SUCCESS="yes"


rm -rf $TEMPDIR

if [ -z $SUCCESS ]; then
    exit 1
else
    exit 0
fi