## Purpose
The purpose of the arrayIO design is to allow for code sharing fragment/matrix
data across 3 different storage formats:

- In-memory arrays
- Binary files arranged in folders on disk
- Single hdf5 file

Because the data can mostly be represented by arrays of unsigned 32-bit integers, a lot of the logic for reading/writing is shared

## Data types

The data types are a little complicated right now since I'm erring on the side of templating rather than explicit inheritance with virtual function calls

#### High-level datatype representations

Fragments have the following main data fields:
  - `cell`, `start`, `end` - cell, start, and end values concatenated between all the chromosomes, ordered to be sorted by `start` within each chromosme
  - `chr_ptr` - indices into the `cell`,`start`, and `end` arrays that define the start/end of each chromosome. Chromosome `i` goes from `chr_ptr[2*i]` to `chr_ptr[2*i+1]`. Length is twice the number of chromosomes. `chr_ptr` is organized
  in this way because fragments.tsv.gz files can't always be seeked efficiently, 
  yet we also want to allow re-ordering chromosome IDs efficiently. 
  This means there are cases when fragments will seem to be loaded not in order by chromosome ID, so we need
  to be able to store data in the order it comes in, while being able to later load chromosomes in the correct order.
  This is not an issue for matrices, since all our input matrix formats are required to allow skipping arbitrarily 
  between columns.
  - `end_max` - running maximum of the end values within each chromosome, stored for every 128 values to enable seeking. `end_max[i]` is the maximum of `end[chr_start_idx:i*128]`. Note that because it's possble a chromosome ends at not a perfect multiple of 128, `end_max` will default to be maximum of `end` for the whole chunk of 128, rather than the maximum of `end` in the start of the new chromosome

Matrices have the following main data fields:
  - `val`, `row` - list of the value and row for each non-zero entry
  - `col_ptr` - Values in column `i` of the matrix go from index `col_ptr[i]` to `col_ptr[i+1]`.
    `col_ptr` is length number of columns + 1
  - `row_count` is a single integer listing the number of rows in the matrix
  

#### Data Loading classes
The main datatypes we want to represent are lists of integers, while also having support for
storing lists of strings (cell IDs, column/row names, etc.). Documentation for these classes is 
provided in the `array_interfaces.h` header file, but this is the high-level summary:

- `UIntBulkReader` and `UIntBulkWriter` are the low-level interfaces for loading data from disk, memory, etc.
- `UIntReader` and `UIntWriter` provide zero-copy, buffered interfaces for reading/writing data. This means
  that rather than copying data to/from a client array, they give clients a pointer to where data can be read/written
- `StringReader` and `StringWriter` provide (potentially inefficient) ways to persist lists of strings to disk.
  Most implementations will have an intermediate step where all IDs are stored in memory
- `ReaderBuilder` and `WriterBuilder` are interfaces to handle reading+writing a set of related UInt and String arrays.
  So e.g. for fragments stored in an hdf5 file, the ReaderBuilder and WriterBuilder would be responsible for opening
  the correct HDF5 DataSet for reading/writing.