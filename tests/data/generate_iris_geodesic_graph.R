#Setup steps:
#python -m venv venv
#source ./venv/bin/activate
#python -m pip install umap-learn

# Make sure we're on the latest version of reticulate: https://github.com/dynverse/anndata/issues/25
# install.packages("reticulate") 


library(reticulate)
use_virtualenv("./venv")
umap <- import("umap")
pynndescent <- import("pynndescent")

data <- as.matrix(iris[,1:4])
pynn <- pynndescent$NNDescent(data)$neighbor_graph
knn <- list(idx = pynn[[1]][,1:5] + 1L, dist = pynn[[2]][,1:5])

x <- umap$umap_$fuzzy_simplicial_set(
    data,
    n_neighbors = 5,
    random_state = NULL,
    metric = NULL,
    knn_indices=np_array(knn$idx - 1L),
    knn_dists=np_array(knn$dist)
    # Don't include these already-set defaults
    # set_op_mix_ratio=1.0,
    # local_connectivity=1.0,
)
sigma <- x[[2]]
rho <- x[[3]]
x <- x[[1]]

saveRDS(list(
    knn=knn,
    graph=x
), "tests/data/iris_geodesic_graph.rds")