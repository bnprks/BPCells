#pragma once
#include "../utils/radix_sort.h"
#include "MatrixIterator.h"

namespace BPCells {

// Rank-transform each column of a matrix. Offset the ranking such that 0 values
// are assigned a rank of 0. Tied ranks are averaged together.
template <class T> class ColwiseRank : public MatrixLoader<double> {
  private:
    std::unique_ptr<MatrixLoader<T>> loader;
    uint32_t current_col = UINT32_MAX - 1;
    bool take_ownership = true;

    std::vector<uint32_t> row_data, row_buf;
    std::vector<T> val_data, val_buf;
    std::vector<double> ranks;
    int64_t tie_statistic;

    uint32_t idx = 0;
    uint32_t cap = 0;
    uint32_t load_size;

  public:
    ColwiseRank(std::unique_ptr<MatrixLoader<T>> &&loader, uint32_t load_size = 1024)
        : loader(std::move(loader))
        , load_size(load_size) {
        row_data.resize(this->loader->rows());
        row_buf.resize(this->loader->rows());
        val_data.resize(this->loader->rows());
        val_buf.resize(this->loader->rows());
        ranks.resize(this->loader->rows());
    }

    ~ColwiseRank() {
        if (!take_ownership) loader.release();
    }
    ColwiseRank() = default;
    ColwiseRank(ColwiseRank&&) = default;
    ColwiseRank& operator=(ColwiseRank&&) = default;

    // Set the object so that the inner loader will be preserved
    // rather than calling the destructor when this loader is destructed
    void preserve_input_loader() {take_ownership = false;}

    // Get this column's sum of t^3 - t
    // where t is the number of tied entries at a given rank
    int64_t tieStatistic() {return tie_statistic;}

    uint32_t rows() const override { return loader->rows(); }
    uint32_t cols() const override { return loader->cols(); }

    const char *rowNames(uint32_t row) override { return loader->rowNames(row); }
    const char *colNames(uint32_t col) override { return loader->colNames(col); }

    void restart() override {
        idx = 0;
        cap = 0;
        this->loader->restart();
    }

    void seekCol(uint32_t col) override {
        idx = 0;
        cap = 0;
        this->loader->seekCol(col);
    }

    bool nextCol() override {
        idx = 0;
        cap = 0;
        return this->loader->nextCol();
    }

    uint32_t currentCol() const override { return loader->currentCol(); }

    bool load() override {
        if (idx == 0 && cap == 0) {
            tie_statistic = 0;
            uint32_t negative_count = 0;
            uint32_t explicit_zeros_count = 0;
            // Load all the values for the column and sort
            while (this->loader->load()) {
                uint32_t loaded = this->loader->capacity();

                uint32_t *row_ptr = this->loader->rowData();
                T *val_ptr = this->loader->valData();

                std::memmove(row_data.data() + cap, row_ptr, sizeof(uint32_t) * loaded);
                for (uint32_t i = 0; i < loaded; i++) {
                    val_data[cap + i] = val_ptr[i];
                    negative_count += val_ptr[i] < 0;
                    explicit_zeros_count += val_ptr[i] == 0;
                }
                cap += loaded;
            }

            // Calculate the offset to apply to make rank of 0 map to 0
            uint32_t implicit_zeros_count = rows() - cap;
            double zero_rank =
                negative_count + (1 + explicit_zeros_count + implicit_zeros_count) / 2.0;
            if (explicit_zeros_count == 0 && implicit_zeros_count == 0) zero_rank = 0;

            // Sort the entries by value
            lsdRadixSortArrays<T, uint32_t>(cap, val_data, row_data, val_buf, row_buf);

            // Loop over negatives
            size_t i = 0;
            size_t ties = 1;
            for (; i < cap && val_data[i] < 0; i += ties) {
                ties = 1;
                while (i + ties < cap && val_data[i] == val_data[i + ties])
                    ties++;
                if (ties > 1) tie_statistic += ties * ties * ties - ties;
                for (size_t j = 0; j < ties; j++) {
                    ranks[i + j] = i + 1 + (ties - 1) / 2.0 - zero_rank;
                }
            }
            // Loop over explicit zeros
            for (; i < cap && val_data[i] == 0; i += 1) {
                ranks[i] = 0;
            }
            double total_zeros = implicit_zeros_count + explicit_zeros_count;
            tie_statistic += total_zeros * total_zeros * total_zeros - total_zeros;
            // Loop over positives
            for (; i < cap; i += ties) {
                ties = 1;
                while (i + ties < cap && val_data[i] == val_data[i + ties])
                    ties++;
                if (ties > 1) tie_statistic += ties * ties * ties - ties;
                for (size_t j = 0; j < ties; j++) {
                    ranks[i + j] = i + 1 + (ties - 1) / 2.0 + implicit_zeros_count - zero_rank;
                }
            }
        } else {
            idx += std::min(cap - idx, load_size);
        }
        return idx < cap;
    }

    uint32_t capacity() const override { return std::min(cap - idx, load_size); }
    uint32_t *rowData() override { return row_data.data() + idx; }
    double *valData() override { return ranks.data() + idx; }
};

} // end namespace BPCells