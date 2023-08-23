#pragma once

// Include std::filesystem, with proper support for the experimental/filesystem fallback headers
// Also make the namespace BPCells::std_fs an alias for std::filesystem or std::experimental::filesystem as needed

// Auto-detection of header (works for gcc + clang I think)
#ifdef __has_include
#if !__has_include(<filesystem>)
#define BPCELLS_CPP17_EXPERIMENTAL_FILESYSTEM
#endif
#endif

// Manual option to force experimental/filesystem
#ifdef BPCELLS_CPP17_EXPERIMENTAL_FILESYSTEM
#include <experimental/filesystem>
namespace BPCells {
    namespace std_fs = std::experimental::filesystem;
}
#else
#include <filesystem>
namespace BPCells {
    namespace std_fs = std::filesystem;
}
#endif