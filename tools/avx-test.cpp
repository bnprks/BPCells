#include "../src/lib/sleef_wrapper.h"

// AVX parameters on Windows can cause seg faults with GCC:
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54412
// This is enough to reproduce the bug at config time

void foo(BPCells::vec_float x) {}

int main() {
    BPCells::vec_float r = BPCells::splat_float(0.0f);
    foo(r);
}