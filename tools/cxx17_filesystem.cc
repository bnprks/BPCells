#include <iostream>
#include "../src/utils/filesystem_compat.h"

int main() {
    BPCells::std_fs::path p("test/path.txt");
    std::cout << p << std::endl;
}