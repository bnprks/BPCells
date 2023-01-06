#include <iostream>
#include <filesystem>

int main() {
    std::filesystem::path p("test/path.txt");
    std::cout << p << std::endl;
}