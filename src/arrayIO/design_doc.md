## Purpose
The purpose of the arrayIO design is to allow for code sharing fragment/matrix
data across 3 different storage formats:

- In-memory arrays
- Binary files arranged in folders on disk
- Single hdf5 file

Because the data can mostly be represented by arrays of unsigned 32-bit integers, a lot of the logic for reading/writing is shared

NOTE: The interface API definitions are out-of-date right now since I just changed how things
work to simplify and make more efficient. Hopefully I'll come back to update this with
the new UIntReader and UIntWriter interfaces

## Data types

The data types are a little complicated right now since I'm erring on the side of templating rather than explicit inheritance with virtual function calls

#### Data representation

- `UnpackedFrags<T>` and `PackedFrags<T>`: Flexible struct representing fragment data for a single chromosome. `T` can be any kind of class -- `vector` `UIntWriter`, `UIntReader` depending on what is needed for the occasion.
  - `UnpackedFrags` includes start, end, cell_id arrays, as well as an end_max array where `end_max[i]` = `max(end[0:(i+1)*128])` [this is used to enable seeking]
  - `PackedFrags` is a bit more complicated due to the bitpacking
    - `start_data`, `end_data`, and `cell_data` store the main bitpacked data arrays. Starts are stored with the BP128 D1 encoding, ends are stored as BP128 of (end-start), and cells are a vanilla BP128 compression
    - `start_idx`, `end_idx`, and `cell_idx` contain offsets for each chunk of 128 integers within the `*_data` arrays. So entries `n*128 -> (n+1)*128` of `start` are read starting at `start_data[start_idx[n]]` and end at `start_data[start_idx[n+1]] `. Note that these arrays are length nchunks + 1
    - `start_starts` contains the beginning value of each chunk of 128 start values, needed for the difference encoding
    - `end_max` is defined the same as for `UnpackedFrags`
    - `count` stores the total number of fragments as a single integer, required so we can represent fragment lists that don't have an even multiple of 128 fragments

#### Headline classes

- `UnpackedFragmentsReader<Loader>` - Unpacked fragments reader object that relies on a given `FragmentsLoader` class responsible for interfacing with the underlying storage.
  - `PackedFragmentsReader<Loader>` - Same idea, but for packed fragments
- `UnpackedFragmentsWriter<Saver>` - Fragments writer object that relies on a given `FragmentsSaver` class responsible for interfacing with the underlying storage
  - `UnpackedFragmentsReader<Loader>` - Same idea, but for packed fragments

#### Helper classes

- `*UIntWriter` - writable stream of unsigned 32-bit integers
  - `void write(const uint32_t *buffer, uint32_t count)` - appends `count` integers read from `buffer` to the end of the stream.
  - `void finalize()` - Finalizes the stream (e.g. flushing any remaining data)
  - This class should support cheap move-assignment and move-constructors, but need-not support cheap copies
- `*UIntReader` - readable stream of unsigned 32-bit integers
  - `uint32_t read(uint32_t *buffer, uint32_t count)` - reads the up to `count` integers from the stream and writes to `buffer`. Returns how many are read, repeatedly returning 0 after stream is empty. 
  - `uint32_t size()` - Return the total number of integers in the array.
  - `void seek(const size_t pos)` - Change the stream to read the next integer from index `pos` (index 0 is first integer in the stream)
  - This class should support cheap move-assignment and move-constructors, but need-not support cheap copies
- `*FragmentsSaver` - Object that manages saving fragments data for each chromosome into the corresponding `*UIntWriter` objects. This can be split into separate objects for Packed and Unpacked fragments, or combined into one that implements all the methods
  - `UnpackedFrags<*UIntWriter> chrWriterUnpacked(uint32_t chr_id)` Initialize + return all the required arrays for writing unpacked fragments for a new chromosome
  - `PackedFrags<*UIntWriter> chrWriterPacked(uint32_t chr_id)` Initialize + return all the required arrays for writing packed fragments for a new chromosome
  - `void writeCellNames(std::vector<std::string> cell_names)`, and  `void writeChrNames(std::vector<std::string> chr_names)` - save the cell_names and chr_names, pretty self-explanatory
- `*FragmentsLoader` - Object that manages loading fragments data for each chromosome from `*UIntReader` objects.
  - `UnpackedFrags<*UIntWriter> chrReaderUnpacked(uint32_t chr_id)` Initialize + return all the required arrays for reading unpacked fragments for a new chromosome
  - `PackedFrags<*UIntWriter> chrReaderPacked(uint32_t chr_id)` Initialize + return all the required arrays for reading packed fragments for a new chromosome
  - `std::vector<std::string> readCellNames()`, and  `std::vector<std::string> readChrNames()` - read the cell_names and chr_names, pretty self-explanatory