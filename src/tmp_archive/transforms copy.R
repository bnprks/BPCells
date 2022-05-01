

setClass("NormalizedMatrix",
    slots = c(
        matrix = "IterableMatrix",
        row_params = "matrix",
        col_params = "matrix",
        global_params = "numeric",
        is_fit = "logical"
    ),
    prototype = list(
        matrix = new("IterableMatrix"),
        row_params = matrix(0,0,0),
        col_params = matrix(0,0,0),
        global_params = numeric(0),
        is_fit = FALSE
    )
)

setGeneric("load_transform", function(x) standardGeneric("load_transform"))

setMethod("%*%", signature(x="NormalizedMatrix", y="matrix"), function(x, y) {
    ptr <- load_transform(x)
    return(transform_dense_multiply_right_cpp(ptr, y))
})

setMethod("%*%", signature(x="matrix", y="NormalizedMatrix"), function(x, y) {
    ptr <- load_transform(y)
    return(transform_dense_multiply_left_cpp(ptr, x))
})

setMethod("%*%", signature(x="NormalizedMatrix", y="numeric"), function(x, y) {
    ptr <- load_transform(x)
    return(transform_vec_multiply_right_cpp(ptr, y))
})

setMethod("%*%", signature(x="numeric", y="NormalizedMatrix"), function(x, y) {
    ptr <- load_transform(y)
    return(transform_vec_multiply_left_cpp(ptr, x))
})

assign_fit <- function(transform) {
    assert_is(transform, "NormalizedMatrix")
    ptr <- load_transform(transform)
    fit <- transform_get_fit_cpp(ptr)
    
    transform@row_params <- fit$row_params
    transform@col_params <- fit$col_params
    transform@global_params <- fit$global_params
    transform@is_fit <- TRUE
    return(transform)
}

setClass("TFIDF", contains="NormalizedMatrix",
    slots = c(
        scaleTo = "numeric"
    ),
    prototype = list(
        scaleTo = numeric(1)
    )
)
setMethod("load_transform", "TFIDF", function(x) {
    if (!x@is_fit)
        transform_tfidf_cpp(iterate_matrix(x@matrix), x@scaleTo, x@matrix@transpose)
    else
        transform_project_tfidf_cpp(iterate_matrix(x@matrix), x, 0, x@matrix@transpose)
})

#' TFIDF normalization
#' @param mat IterableMatrix to transform
#' @param scaleTo Value to scale to
#' @return NormalizedMatrix object
#' @details Transforms values of matrix according to the formula
#'     log(1 + scaleTo * x / (colSum * rowMean))
#' @export
normalize_TFIDF <- function(mat, scaleTo=1e4) {
    assert_is(mat, "IterableMatrix")
    assert_is(scaleTo, "numeric")

    ret <- new("TFIDF", matrix=mat, scaleTo=scaleTo)
    ret <- assign_fit(ret)

    return(ret)
}