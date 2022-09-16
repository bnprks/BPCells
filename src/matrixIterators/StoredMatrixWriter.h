#pragma once

#include "../arrayIO/array_interfaces.h"
#include "../arrayIO/bp128.h"
#include "MatrixIterator.h"
#include "OrderRows.h"
#include "StoredMatrix.h"

namespace BPCells {

// Class for writing a column-major MatrixLoader as either a column-major
// or row-major disk format. Row-major outputs will be the transpose of the
// column-major inputs
template <typename T> class StoredMatrixWriter : public MatrixWriter<T> {
  private:
    UIntWriter row, col_ptr, shape;
    NumWriter<T> val;
    std::unique_ptr<StringWriter> row_names, col_names, storage_order;
    bool row_major;

  public:
    StoredMatrixWriter(
        UIntWriter &&row,
        NumWriter<T> &&val,
        UIntWriter &&col_ptr,
        UIntWriter &&shape,
        std::unique_ptr<StringWriter> &&row_names,
        std::unique_ptr<StringWriter> &&col_names,
        std::unique_ptr<StringWriter> &&storage_order,
        bool row_major = false
    )
        : row(std::move(row))
        , col_ptr(std::move(col_ptr))
        , shape(std::move(shape))
        , val(std::move(val))
        , row_names(std::move(row_names))
        , col_names(std::move(col_names))
        , storage_order(std::move(storage_order))
        , row_major(row_major) {}

    static StoredMatrixWriter createUnpacked(WriterBuilder &wb, bool row_major = false) {
        wb.writeVersion(StoredMatrix<T>::versionString(false));
        return StoredMatrixWriter(
            wb.createUIntWriter("index"),
            wb.create<T>("val"),
            wb.createUIntWriter("idxptr"),
            wb.createUIntWriter("shape"),
            wb.createStringWriter("row_names"),
            wb.createStringWriter("col_names"),
            wb.createStringWriter("storage_order"),
            row_major
        );
    }

    static StoredMatrixWriter
    createPacked(WriterBuilder &wb, bool row_major = false, uint32_t buffer_size = 1024) {
        wb.writeVersion(StoredMatrix<T>::versionString(true));
        NumWriter<T> val;

        if constexpr (std::is_same_v<T, uint32_t>) {
            val = UIntWriter(
                std::make_unique<BP128_FOR_UIntWriter>(
                    wb.createUIntWriter("val_data"), wb.createUIntWriter("val_idx")
                ),
                buffer_size
            );
        } else {
            val = wb.create<T>("val");
        }

        return StoredMatrixWriter(
            UIntWriter(
                std::make_unique<BP128_D1Z_UIntWriter>(
                    wb.createUIntWriter("index_data"),
                    wb.createUIntWriter("index_idx"),
                    wb.createUIntWriter("index_starts")
                ),
                buffer_size
            ),
            std::move(val),
            wb.createUIntWriter("idxptr"),
            wb.createUIntWriter("shape"),
            wb.createStringWriter("row_names"),
            wb.createStringWriter("col_names"),
            wb.createStringWriter("storage_order"),
            row_major
        );
    }

    void write(MatrixLoader<T> &mat_in, void (*checkInterrupt)(void) = NULL) override {
        // Ensure that we write matrices sorted by row
        OrderRows<T> mat(mat_in);
        uint32_t col = 0;
        uint32_t idx = 0; // Index of for col_ptr array

        uint32_t max_capacity = std::max(row.maxCapacity(), val.maxCapacity());

        col_ptr.write_one(idx);

        while (mat.nextCol()) {
            if (checkInterrupt != NULL) checkInterrupt();
            if (mat.currentCol() < col)
                throw std::runtime_error("StoredMatrixWriter encountered out-of-order columns");
            while (col < mat.currentCol()) {
                col_ptr.write_one(idx);
                col += 1;
            }
            while (mat.load()) {
                uint32_t i = 0;
                while (i < mat.capacity()) {
                    uint32_t capacity = std::min(max_capacity, mat.capacity() - i);

                    row.ensureCapacity(capacity);
                    val.ensureCapacity(capacity);

                    std::memmove(row.data(), mat.rowData() + i, capacity * sizeof(uint32_t));
                    std::memmove(val.data(), mat.valData() + i, capacity * sizeof(T));

                    row.advance(capacity);
                    val.advance(capacity);
                    idx += capacity;
                    i += capacity;
                }

                if (checkInterrupt != NULL) checkInterrupt();
            }
        }
        if (row_major) {
            shape.write_one(mat.cols());
            shape.write_one(mat.rows());
        } else {
            shape.write_one(mat.rows());
            shape.write_one(mat.cols());
        }
        col_ptr.write_one(idx);

        row.finalize();
        val.finalize();
        col_ptr.finalize();
        shape.finalize();

        // Get row and col names. This probably incurs a few extra copies,
        // but it shouldn't matter since writing the actual matrix should dominate cost
        std::vector<std::string> col_names(0);
        for (int i = 0;; i++) {
            const char *col_name = mat.colNames(i);
            if (col_name == NULL) break;
            col_names.push_back(std::string(col_name));
        }

        std::vector<std::string> row_names(0);
        for (int i = 0;; i++) {
            const char *row_name = mat.rowNames(i);
            if (row_name == NULL) break;
            row_names.push_back(std::string(row_name));
        }

        if (row_major) std::swap(row_names, col_names);
        this->col_names->write(VecStringReader(col_names));
        this->row_names->write(VecStringReader(row_names));

        std::vector<std::string> storage_order;
        storage_order.push_back(row_major ? "row" : "col");
        this->storage_order->write(VecStringReader(storage_order));
    }
};

} // end namespace BPCells