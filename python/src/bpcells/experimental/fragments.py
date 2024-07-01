# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.

import bpcells
import bpcells.cpp

import json
import tempfile
import os.path

from typing import Dict, List, Optional, Tuple, Union
import sys
if sys.version_info >= (3, 9):
    from collections.abc import Sequence
else:
    from typing import Sequence

import numpy as np
import pandas as pd

def import_10x_fragments(input: str, output: str, shift_start: int = 0, shift_end: int = 0, keeper_cells: Optional[List[str]] = None):
    """Convert 10x fragment file to BPCells format

    Args:
        input (str): Path to 10x input file
        output (str): Path to BPCells output directory
        shift_start (int): Basepairs to add to start coordinates (generally positive number)
        shift_end (int): Basepairs to subtract from end coordinates (generally negative number)
        keeper_cells (list[str]): If not None, only save fragments from cells in the keeper_cells list
    """
    keeper_cells = np.asarray(keeper_cells) if keeper_cells is not None else keeper_cells
    bpcells.cpp.import_10x_fragments(input, output, shift_start, shift_end, keeper_cells)

def build_cell_groups(fragments: str, cell_ids: Sequence[str], group_ids: Sequence[str], group_order: Sequence[str]) -> np.ndarray:
    """Build cell_groups array for use in :func:`pseudobulk_insertion_counts()`

    Args:
        fragments (str): Path to BPCells fragments directory
        cell_ids (list[str]): List of cell IDs
        group_ids (list[str]): List of pseudobulk IDs for each cell (same length as ``cell_ids``)
        group_order (list[str]): Output order of pseudobulks (Contain the unique ``group_ids``)
    
    Returns:
        numpy.ndarray: 
        Numpy array suitable as input for ``cell_groups`` in :func:`pseudobulk_insertion_counts()`.
        Same length as total number of cells in the ``fragments`` input, specifying the output
        pseudobulk index for each cell (or -1 if the cell is excluded from consideration)
        
    See Also:
        :func:`pseudobulk_insertion_counts`
    """
    cell_index_lookup = {c: i for i, c in enumerate(bpcells.cpp.cell_names_fragments_dir(fragments))}
    group_index_lookup = {g: i for i, g in enumerate(group_order)}
    
    assert len(cell_ids) == len(group_ids)
    assert set(group_ids) <= set(group_index_lookup.keys())
    
    ret = np.full(len(cell_index_lookup), -1, np.int32)
    for cell_id, group_id in zip(cell_ids, group_ids):
        ret[cell_index_lookup[cell_id]] = group_index_lookup[group_id]
    
    return ret

def pseudobulk_insertion_counts(fragments: str, regions: pd.DataFrame, cell_groups: Sequence[int], bin_size: int = 1) -> np.ndarray:
    """Calculate a pseudobulk coverage matrix

    Coverage is calculated as the number of start/end coordinates falling into a given position bin.

    Args:
        fragments (str): Path to BPCells fragments directory
        regions (pandas.DataFrame): Pandas dataframe with columns (``chrom``, ``start``, ``end``) representing
          genomic ranges (0-based, end-exclusive like BED format). All regions must be the same size.
          ``chrom`` should be a string column; ``start``/``end`` should be numeric.
        cell_groups (list[int]): List of pseudbulk groupings as created by :func:`build_cell_groups()`
        bin_size (int): Size for bins within each region given in basepairs. If the region width is not
          an even multiple of ``resolution_bp``, then the last region may be truncated.
    
    Returns:
        numpy.ndarray: Numpy array with dimensions (region, psudobulks, position) and type numpy.int32
    
    See Also: 
        :func:`build_cell_groups`
    """
    chrs = bpcells.cpp.chr_names_fragments_dir(fragments)

    peak_order = sorted(
        range(regions.shape[0]),
        key = lambda i: (chrs.index(regions["chrom"].iloc[i]), regions["start"].iloc[i])
    )

    regions = regions.iloc[peak_order,]

    mat = bpcells.cpp.pseudobulk_coverage(
        fragments,
        np.asarray(regions["chrom"]),
        np.asarray(regions["start"]),
        np.asarray(regions["end"]),
        np.asarray(cell_groups),
        bin_size
    )
    return mat.reshape((mat.shape[0], -1, regions.shape[0]), order="F").transpose(2,0,1)


class PrecalculatedInsertionMatrix:
    """
    Disk-backed precalculated insertion matrix

    This reads per-base precalculated insertion matrices. The current implementation is EXPERIMENTAL, and will crash for matrices with more than
    2^32-1 non-zero entries.

    Args:
        dir (str): Path of the matrix directory
        
    See Also:
        :func:`precalculate_insertion_counts`
    """
    def __init__(self, path: str):
        self._dir = str(os.path.abspath(os.path.expanduser(path)))
        
        self._chrom_offsets = json.load(open(f"{self._dir}/chrom_offsets.json"))

    @property
    def shape(self) -> Tuple[int, int]:
        return tuple(np.fromfile(f"{self._dir}/shape", np.uint32, 2, offset=8))

    def __repr__(self):
        return f"<PrecalculatedInsertionMatrix with {self.shape[0]} pseudobulks and {len(self._chrom_offsets)} chromomsomes stored in \n\t{self._dir}"

    def get_counts(self, regions: pd.DataFrame):
        """Load pseudobulk insertion counts

        Args:
            regions (pandas.DataFrame): Pandas dataframe with columns (``chrom``, ``start``, ``end``) representing
                genomic ranges (0-based, end-exclusive like BED format). All regions must be the same size.
                ``chrom`` should be a string column; ``start``/``end`` should be numeric.
        
        Returns:
            numpy.ndarray: Numpy array of dimensions (region, psudobulks, position) and type numpy.int32
        """
        region_size = regions.end.iloc[0] - regions.start.iloc[0]
        assert (regions.end - regions.start == region_size).all()
        
        start_indices = [
            self._chrom_offsets[t.chrom] + t.start for t in regions.itertuples()
        ]
        
        return bpcells.cpp.query_precalculated_pseudobulk_coverage(
            self._dir,
            start_indices,
            region_size
        )\
            .reshape((region_size, regions.shape[0], -1), order="F")\
            .transpose((1,2,0))

def precalculate_insertion_counts(fragments: str, output_dir: str, cell_groups: Sequence[int], chrom_sizes: Union[str, Dict[str, int]], threads: int = 0):
    """Precalculate per-base insertion counts from fragment data

    The current implementation is EXPERIMENTAL, and will crash for matrices with more than
    2^32-1 non-zero entries.

    Args:
        fragments (str): Path to a BPCells fragments directory
        output_dir (str): Path to save the insertion counts in
        cell_groups (list[int]): List of pseudbulk groupings as created by :func:`build_cell_groups()`
        chrom_sizes (str | dict[str, int]): Path/URL of UCSC-style chrom.sizes file, or dictionary mapping chromosome names to sizes
        threads (int): Number of threads to use during matrix calculation (default = 1)
    
    Returns:
        A :class:`PrecalculatedInsertionMatrix` object

    See Also:
        :class:`PrecalculatedInsertionMatrix`
    """
    if isinstance(chrom_sizes, str):
        chrom_sizes = pd.read_csv(chrom_sizes, sep="\t", names=["chrom", "size"])
        chrom_sizes = {t.chrom: t.size for t in chrom_sizes.itertuples()}
    
    # Re-order chrom_sizes to match the fragment file chromosome order
    chrom_order = bpcells.cpp.chr_names_fragments_dir(fragments)
    chrom_sizes = dict(i for i in chrom_sizes.items() if i[0] in chrom_order)
    chrom_sizes = dict(sorted(chrom_sizes.items(), key = lambda x: chrom_order.index(x[0])))

    tmp = tempfile.TemporaryDirectory()
    bpcells.cpp.precalculate_pseudobulk_coverage(
        fragments,
        output_dir,
        tmp.name,
        list(chrom_sizes.keys()),
        list(chrom_sizes.values()),
        cell_groups,
        1,
        threads
    )
    
    chrom_offsets = dict(zip(chrom_sizes.keys(), [0] + np.cumsum(list(chrom_sizes.values()))[:-1].tolist()))
    json.dump(chrom_offsets, open(f"{output_dir}/chrom_offsets.json", "w"), indent=2)
    return PrecalculatedInsertionMatrix(output_dir)