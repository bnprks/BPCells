#include "StoredMatrix.h"

namespace BPCells {

StoredMatrixWriter::StoredMatrixWriter(UIntWriter &&row, UIntWriter &&val, UIntWriter &&col_ptr, 
    UIntWriter &&row_count, 
    std::unique_ptr<StringWriter> &&row_names, std::unique_ptr<StringWriter> &&col_names) :
    row(std::move(row)), val(std::move(val)), col_ptr(std::move(col_ptr)), row_count(std::move(row_count)),
    row_names(std::move(row_names)), col_names(std::move(col_names)) {}

StoredMatrixWriter StoredMatrixWriter::createUnpacked(WriterBuilder &wb) {
    wb.writeVersion("unpacked-uint-matrix-v1");
    return StoredMatrixWriter(
        wb.createUIntWriter("row"),
        wb.createUIntWriter("val"),
        wb.createUIntWriter("col_ptr"),
        wb.createUIntWriter("row_count"),
        wb.createStringWriter("row_names"),
        wb.createStringWriter("col_names")
    );
}
StoredMatrixWriter StoredMatrixWriter::createPacked(WriterBuilder &wb, uint32_t buffer_size) {
    wb.writeVersion("packed-uint-matrix-v1");
    return StoredMatrixWriter(
        UIntWriter(
            std::make_unique<BP128_D1Z_UIntWriter>(
                wb.createUIntWriter("row_data"),
                wb.createUIntWriter("row_idx"),
                wb.createUIntWriter("row_starts")
            ),
            buffer_size
        ),
        UIntWriter(
            std::make_unique<BP128_FOR_UIntWriter>(
                wb.createUIntWriter("val_data"),
                wb.createUIntWriter("val_idx")
            ),
            buffer_size
        ),
        wb.createUIntWriter("col_ptr"),
        wb.createUIntWriter("row_count"),
        wb.createStringWriter("row_names"),
        wb.createStringWriter("col_names")
    );
}
void StoredMatrixWriter::write(MatrixIterator<uint32_t> &mat, void (*checkInterrupt)(void)) {
    uint32_t col = 0;
    uint32_t idx = 0; // Index of for col_ptr array

    uint32_t max_capacity = std::max(row.maxCapacity(), val.maxCapacity());

    col_ptr.write_one(idx);

    while (mat.nextCol()) {
        if(checkInterrupt != NULL) checkInterrupt();
        if (mat.currentCol() < col)
            throw std::runtime_error("StoredMatrixWriter encountered out-of-order columns");
        while (col < mat.currentCol()) {
            col_ptr.write_one(idx);
            col += 1;
        }
        while(mat.load()) {
            uint32_t load_remaining = mat.capacity();
            while (load_remaining > 0) {
                uint32_t capacity = std::min(max_capacity, load_remaining);

                row.ensureCapacity(capacity);
                val.ensureCapacity(capacity);

                std::memmove(row.data(), mat.rowData(), capacity*sizeof(uint32_t));
                std::memmove(val.data(), mat.valData(), capacity*sizeof(uint32_t));

                row.advance(capacity);
                val.advance(capacity);
                idx += capacity;

                load_remaining -= capacity;
            }            

            if(checkInterrupt != NULL) checkInterrupt();
        }
    }
    row_count.write_one(mat.rows());
    col_ptr.write_one(idx);

    row.finalize();
    val.finalize();
    col_ptr.finalize();
    row_count.finalize();

    // Get cell and chromosome names. This probably incurs a few extra copies,
    // but it shouldn't matter since writing the actual fragments should dominate cost
    std::vector<std::string> col_names(0);
    for (int i = 0; ; i++) {
        const char* col_name = mat.colNames(i);
        if (col_name == NULL) break;
        col_names.push_back(std::string(col_name));
    }
    this->col_names->write(VecStringReader(col_names));
    
    std::vector<std::string> row_names(0);
    for (int i = 0; ; i++) {
        const char* row_name = mat.rowNames(i);
        if (row_name == NULL) break;
        row_names.push_back(std::string(row_name));
    }
    this->row_names->write(VecStringReader(row_names));
}

} // end namespace BPCells