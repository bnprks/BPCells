#include <iostream>
#include "../src/bpcells-cpp/utils/filesystem_compat.h"

int main() {
    BPCells::std_fs::path p("test/path.txt");
    std::cout << p << std::endl;
}