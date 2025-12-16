import numpy as np
import anndata as ad
import pandas as pd
import scipy

if __name__ == "__main__":
    # Reproducibility
    rng = np.random.default_rng(42)

    # Dimensions
    n_cells = 5
    n_genes = 6

    # Create int64 matrix with negative values
    X_neg = rng.integers(
        low=-100,
        high=100,
        size=(n_cells, n_genes),
        dtype=np.int64
    )
    X = rng.integers(
        low = 0,
        high = 100,
        size = (n_cells, n_genes),
        dtype = np.int64
    )

    # Optional: cell and gene metadata
    obs = pd.DataFrame(
        {
            "cell_type": rng.choice(["A", "B", "C"], size=n_cells),
            "batch": rng.integers(1, 4, size=n_cells),
        },
        index=[f"cell_{i}" for i in range(n_cells)]
    )

    var = pd.DataFrame(
        {
            "gene_length": rng.integers(500, 3000, size=n_genes),
        },
        index=[f"gene_{j}" for j in range(n_genes)]
    )

    # Create AnnData object
    adata_neg = ad.AnnData(
        X=X_neg,
        layers = {"transpose": scipy.sparse.csc_matrix(X_neg), "dense": X_neg},
        obsm = {"obs_mat": scipy.sparse.csc_matrix(X_neg[:,:2])},
        varm = {"var_mat": X_neg[:2,:].T},
    )
    adata = ad.AnnData(
        X=X,
        layers = {"transpose": scipy.sparse.csc_matrix(X), "dense": X},
        obsm = {"obs_mat": scipy.sparse.csc_matrix(X[:,:2])},
        varm = {"var_mat": X[:2,:].T},
        )

    # Sanity checks
    assert adata.X.dtype == np.int64
    assert np.any(adata_neg.X < 0)

    adata_neg.write("mini_mat_int_64_neg.anndata-v0.12.6.h5ad")
    adata.write("mini_mat_int_64.anndata-v0.12.6.h5ad")