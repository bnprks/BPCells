#include <hwy/highway.h>
#include <hwy/aligned_allocator.h>

int main() {
    auto alignedAlloc = hwy::AllocateAligned<float>(1024);
    return 0;
}