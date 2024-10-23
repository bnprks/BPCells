
# The anndata version we want is from July 2019, which is when Python 3.7 was new
# AnnData versions 0.6.22 and 0.7.0 differ from each other and the current format

# Building python 3.7.17 (new at the time)
# wget https://www.python.org/ftp/python/3.7.17/Python-3.7.17.tgz
# tar -xzvf Python-3.7.17.tgz
# cd Python-3.7.17
# ./configure
# make -j 8
# cd ..

# Python-3.7.17/python -m venv venv-3.7
# source ./venv-3.7/bin/activate
# python -m pip install anndata==0.6.22 "pandas<0.25.0"

# For other anndata versions:
#    python -m pip install anndata==0.7.0 "pandas<0.25.0"


import anndata as ad
import numpy as np

import scipy.sparse

x = np.array([[1., 2., 0., 3., 0.],
        [0., 0., 2., 2., 1.],
        [0., 0., 2., 0., 2.]])

adata = ad.AnnData(
    X = scipy.sparse.csr_matrix(x),
    layers = {"transpose": scipy.sparse.csc_matrix(x), "dense": x},
    obsm = {"obs_mat": scipy.sparse.csc_matrix(x[:,:2])},
    varm = {"var_mat": x[:2,:].T},
)

adata.write_h5ad(f"mini_mat.anndata-v{ad.__version__}.h5ad")