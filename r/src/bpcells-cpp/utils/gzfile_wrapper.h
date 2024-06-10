// Copyright 2023 BPCells contributors
// 
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

#pragma once

#include <utility>

#include <zlib.h>

#include "filesystem_compat.h"

namespace BPCells {

// This wrapper helps for opening gzFiles, and it also ensures that a gzFile is properly closed via
// destructor
class gzFileWrapper {
  private:
    gzFile f = NULL;

  public:
    gzFileWrapper() = default;

    gzFileWrapper(const std::string path, std::string mode, uint32_t buffer_size = 1 << 20) {
        std::string str_path(path);

        if (mode == "r") {
            f = gzopen(path.c_str(), "r");
        } else if (mode == "w") {
            // Create directory if it doesn't already exist
            std_fs::path fpath(path);
            if (fpath.has_parent_path() && !std_fs::exists(fpath.parent_path())) {
                std_fs::create_directories(fpath.parent_path());
            }

            size_t extension_idx = str_path.rfind(".");
            if (extension_idx != std::string::npos && str_path.substr(extension_idx) == ".gz") {
                f = gzopen(path.c_str(), "wb1");
            } else {
                f = gzopen(path.c_str(), "wT");
            }
        } else {
            throw std::runtime_error("gzFileWrapper: invalid mode string: " + mode);
        }

        if (f == NULL) {
            int errnum;
            throw std::runtime_error(
                std::string("gzFileWrapper: Could not open file: ") + gzerror(f, &errnum)
            );
        }

        // Note default of large 1MB buffer to speed up reading
        gzbuffer(f, buffer_size);
    }

    ~gzFileWrapper() {
        if (f != NULL) {
            gzclose(f);
        }
    }

    // No copy assignment or construction
    gzFileWrapper(const gzFileWrapper &other) = delete;
    gzFileWrapper &operator=(const gzFileWrapper &other) = delete;

    // Normal move assignment and construction
    gzFileWrapper(gzFileWrapper &&other) noexcept : f(std::exchange(other.f, nullptr)) {}
    gzFileWrapper &operator=(gzFileWrapper &&other) noexcept {
        std::swap(f, other.f);
        return *this;
    }

    inline gzFile operator*() const noexcept { return f; }
};

} // namespace BPCells