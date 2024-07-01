# Copyright 2023 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


import bpcells.experimental

import pytest

import numpy as np
import scipy

@pytest.fixture
def mini_array(tmp_path):
    """Write a small matrix to disk, returning a tuple of (scipy.sparse, path)"""
    row = np.array([0, 2, 2, 0, 1, 2])
    col = np.array([0, 0, 1, 2, 2, 2])
    data = np.array([1, 2, 3, 4, 5, 6])
    mat = scipy.sparse.csc_matrix((data, (row, col)), shape=(3, 4))

    bpcells.experimental.DirMatrix.from_scipy_sparse(mat, str(tmp_path / "mini_array"))

    return (mat, tmp_path / "mini_array")

def sparse_equal(a, b):
    return (a != b).nnz == 0

def test_matrix_slice(mini_array):
    mat, path = mini_array
    bp_mat = bpcells.experimental.DirMatrix(path)

    assert sparse_equal(mat[:,[2,1]], bp_mat[:,[2,1]])
    assert sparse_equal(mat[[2,1],:], bp_mat[[2,1],:])
    assert sparse_equal(bp_mat[:,1:3], mat[:,1:3])
    assert sparse_equal(bp_mat[:,-1:1], mat[:,-1:1])
    assert sparse_equal(bp_mat[:,[True,False,True,False]], mat[:,[True,False,True,False]])
    
    bp_mat.threads = 2
    assert sparse_equal(mat[:,[2,1]], bp_mat[:,[2,1]])

def test_matrix_transpose_dir(mini_array):
    _, path = mini_array
    bp_mat = bpcells.experimental.DirMatrix(path)
    bp_mat_t = bp_mat.T

    assert bp_mat.shape == bp_mat_t.shape[::-1]
    assert sparse_equal(bp_mat[[1,2],:], bp_mat_t[:,[1,2]].T)

def test_matrix_transpose_mem(mini_array):
    _, path = mini_array
    bp_mat = bpcells.experimental.MemMatrix(path)
    bp_mat_t = bp_mat.T

    assert bp_mat.shape == bp_mat_t.shape[::-1]
    assert sparse_equal(bp_mat[[1,2],:], bp_mat_t[:,[1,2]].T)

def test_mem_matrix_slice(mini_array):
    mat, path = mini_array
    bp_mat = bpcells.experimental.MemMatrix(path)

    assert sparse_equal(mat[:,[2,1]], bp_mat[:,[2,1]])
    assert sparse_equal(mat[[2,1],:], bp_mat[[2,1],:])
    assert sparse_equal(bp_mat[:,1:3], mat[:,1:3])
    assert sparse_equal(bp_mat[:,-1:1], mat[:,-1:1])
    assert sparse_equal(bp_mat[:,[True,False,True,False]], mat[:,[True,False,True,False]])
    
    bp_mat.threads = 2
    assert sparse_equal(mat[:,[2,1]], bp_mat[:,[2,1]])

def test_matrix_dims(mini_array):
    mat, path = mini_array
    
    assert mat.shape == bpcells.experimental.DirMatrix(path).shape

def test_scipy_import(tmp_path):
    row = np.array([0, 2, 2, 0, 1, 2])
    col = np.array([0, 0, 1, 2, 2, 2])
    data = np.array([1, 2, 3, 4, 5, 6])
    mat = scipy.sparse.csc_array((data, (row, col)), shape=(3, 4))

    scipy_types = [
        scipy.sparse.bsr_array,
        scipy.sparse.coo_array,
        scipy.sparse.csc_array,
        scipy.sparse.csr_array,
        scipy.sparse.dia_array,
        scipy.sparse.dok_array,
        scipy.sparse.lil_array,
    ]

    for mat_type in scipy_types:
        r = bpcells.experimental.DirMatrix.from_scipy_sparse(mat_type(mat), str(tmp_path / mat_type.__name__))
        assert sparse_equal(r[:,:], mat)
        assert r._transpose == (mat_type == scipy.sparse.csr_array)

def test_hstack(tmp_path, mini_array):
    mat, path = mini_array

    bp_mat = bpcells.experimental.DirMatrix(path)
    
    vstack_1 = bpcells.experimental.DirMatrix.from_vstack([bp_mat]*5, tmp_path / "vstack_1")
    assert sparse_equal(scipy.sparse.vstack([mat]*5), vstack_1[:,:])

    vstack_2 = bpcells.experimental.DirMatrix.from_vstack([bp_mat.T]*5, tmp_path / "vstack_2")
    assert sparse_equal(scipy.sparse.vstack([mat.T]*5).tocsc(), vstack_2[:,:])

    hstack_1 = bpcells.experimental.DirMatrix.from_hstack([bp_mat]*5, tmp_path / "hstack_1")
    assert sparse_equal(scipy.sparse.hstack([mat]*5).tocsc(), hstack_1[:,:])

    hstack_2 = bpcells.experimental.DirMatrix.from_hstack([bp_mat.T]*5, tmp_path / "hstack_2")
    assert sparse_equal(scipy.sparse.hstack([mat.T]*5), hstack_2[:,:])
    

def test_h5ad(tmp_path, mini_array):
    mat, _ = mini_array
    anndata = pytest.importorskip("anndata")
    
    anndata.AnnData(mat).write(tmp_path / "mat1.h5ad")

    bp_mat = bpcells.experimental.DirMatrix.from_h5ad(tmp_path / "mat1.h5ad", tmp_path / "mat1_import")
    assert sparse_equal(bp_mat[:,:], mat)

    anndata.AnnData(mat.tocsr()).write(tmp_path / "mat2.h5ad")

    bp_mat = bpcells.experimental.DirMatrix.from_h5ad(tmp_path / "mat2.h5ad", tmp_path / "mat2_import")
    assert sparse_equal(bp_mat[:,:], mat)
