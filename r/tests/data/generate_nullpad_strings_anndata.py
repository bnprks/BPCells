# Writes an anndata file whose obs/_index and var/_index are variable-length
# string datasets encoded as STRPAD=NULLPAD/CSET=ASCII, which is what R's
# rhdf5 package writes for character vectors (and therefore what
# anndataR::write_h5ad() produces), instead of the STRPAD=NULLTERM/CSET=UTF8
# that anndata's own h5py-based writer always uses. Both are valid HDF5
# encodings of the same data. Written directly with h5py since anndata's
# writer has no option to control STRPAD/CSET.

import h5py
import numpy as np
from scipy import sparse


def write_vlen_strings(group, name, values, strpad, cset):
    tid = h5py.h5t.C_S1.copy()
    tid.set_size(h5py.h5t.VARIABLE)
    tid.set_strpad(strpad)
    tid.set_cset(cset)
    space = h5py.h5s.create_simple((len(values),))
    dset = h5py.Dataset(h5py.h5d.create(group.id, name.encode(), tid, space))
    dset[:] = np.array(values, dtype=object)


if __name__ == "__main__":
    x = np.array([[1., 2., 0., 3., 0.],
                  [0., 0., 2., 2., 1.],
                  [0., 0., 2., 0., 2.]])
    X = sparse.csr_matrix(x)
    n_obs, n_var = X.shape

    with h5py.File("mini_mat_nullpad_strings.anndata.h5ad", "w") as f:
        f.attrs["encoding-type"], f.attrs["encoding-version"] = "anndata", "0.1.0"

        xg = f.create_group("X")
        xg.attrs["encoding-type"], xg.attrs["encoding-version"] = "csr_matrix", "0.1.0"
        xg.attrs["shape"] = [n_obs, n_var]
        xg.create_dataset("data", data=X.data)
        xg.create_dataset("indices", data=X.indices.astype("int32"))
        xg.create_dataset("indptr", data=X.indptr.astype("int32"))

        for axis, n in [("obs", n_obs), ("var", n_var)]:
            g = f.create_group(axis)
            g.attrs["encoding-type"], g.attrs["encoding-version"] = "dataframe", "0.2.0"
            g.attrs["_index"], g.attrs["column-order"] = "_index", []
            write_vlen_strings(
                g,
                "_index",
                [str(i) for i in range(n)],
                h5py.h5t.STR_NULLPAD,
                h5py.h5t.CSET_ASCII,
            )
