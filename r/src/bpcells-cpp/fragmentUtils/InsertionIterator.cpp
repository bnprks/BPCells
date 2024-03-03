#include "InsertionIterator.h"

namespace BPCells {

// I made literally everything else inline, but here we are...
InsertionIterator::InsertionIterator(FragmentLoader &loader) : frags(loader) {}

} // end namespace BPCells
