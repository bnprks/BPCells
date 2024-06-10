
MACRO=$1
TEMPDIR=$(mktemp -d -t macrotest.$1.XXXXX)

echo '
int main() {
    #if !defined('"$MACRO"')
    static_assert(false, "'"$MACRO"' is not defined");
    #endif
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