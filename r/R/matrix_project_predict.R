# Copyright 2024 BPCells contributors
# 
# Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
# https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
# <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
# option. This file may not be copied, modified, or distributed
# except according to those terms.


#' Pipeline Base Class
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @name PipelineBase
#' @export
setClass("PipelineBase",
  contains = "VIRTUAL",
  slots = list(
    fitted = "logical"
  ),
  prototype = list(
    fitted = FALSE
  )
)

#' Fit the pipeline object to data
#' @param object (PipelineBase) The pipeline object to fit.
#' @param x (IterableMatrix) Input data to be fitted on.
#' @param y Optional output data to be fitted on.  Required if the final step is an Estimator, else ignored.
#' @return The fitted pipeline object.
#' @name fit(PipelineBase,IterableMatrix)
#' @export
setGeneric("fit", function(object, x, y = NULL, ...) standardGeneric("fit"))

#' @export
setMethod("fit", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, y = NULL, ...) {
  stop("fit() method not implemented for PipelineBase")
})

#' Transform the input data using a fitted pipeline
#' @param object (PipelineBase) The fitted pipeline object
#' @param x (IterableMatrix) Input data to be transformed
#' @return Data transformed by the pipeline
#' @name transform(PipelineBase,IterableMatrix)
#' @export
setGeneric("transform", function(object, x, ...) standardGeneric("transform"))

#' @export
setMethod("transform", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, ...) {
  stop("transform() method not implemented for PipelineBase")
})

#' Predict the output data using a fitted pipeline
#' @param object (PipelineBase) The fitted pipeline object
#' @param x (IterableMatrix) Input data to be predicted
#' @return Predicted output data
#' @name predict(PipelineBase,IterableMatrix)
#' @export
setGeneric("predict", function(object, x, ...) standardGeneric("predict"))

#' @export
setMethod("predict", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, ...) {
  stop("predict() method not implemented for PipelineBase")
})

#' Combine pipeline objects, to create a new pipeline object.
#' @param x (PipelineBase) The pipeline object to combine to.
#' @param ... (PipelineBase) The pipeline objects to combine from.
#' @name c(PipelineBase)
#' @export
setMethod("c", signature(x = "PipelineBase"), function(x, ...) {
  stop("c() method not implemented for PipelineBase")
})

setMethod("show", signature(object = "PipelineBase"), function(object) {
  stop("show() method not implemented for PipelineBase")
})

#' S4 class for combining multiple pipeline steps into a single pipeline
#' @slot steps (list) List of pipeline steps to execute in order
#' @name Pipeline
setClass(
  "Pipeline",
  contains = "PipelineBase",
  slots = list(
    steps = "list"
  ),
  prototype = list(
    steps = list()
  )
)

#' Return a new Pipeline object
#' @param steps A list of ordered steps to be executed in the pipeline.
#' @return A new Pipeline object.
#' @export
Pipeline <- function(steps = list()) {
  return(new("Pipeline", steps = steps))
}

#' Fit the pipeline object to data
#' @param object (Pipeline) The pipeline object to fit.
#' @param x (IterableMatrix) Input data to be fitted on.
#' @param y Optional output data to be fitted on.  Required if the final step is an Estimator, else ignored.
#' @return The fitted pipeline object.
#' @noRd
#' @export
setMethod("fit", signature(object = "Pipeline", x = "IterableMatrix"), function(object, x, y = NULL, ...) {
  steps <- object@steps
  # Check if all steps are transformers, with the final step being either an estimator or a transformer
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    # allow to fit with estimators as well
    if (i < length(steps)) {
      assert_is(step, "PipelineStep")
    } else {
      assert_is(step, c("PipelineStep"))
      if (!is.null(y)) {
        assert_is(step, "Estimator")
      }
    }
  }
  # Fit every step in the pipeline
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    # allow to fit with estimators as well
    if (i < length(steps) || is.null(y)) {
      step <- fit(step, x, ...)
      x <- transform(step, x)
      if (is(x, "dgCMatrix")) x <- as(x, "IterableMatrix")
    } else {
      step <- fit(step, x, y, ...)
    }
    steps[[i]] <- step
  }
  object@steps <- steps
  object@fitted <- TRUE
  return(object)
})


#' Transform the input data using a fitted pipeline
#' @param object (Pipeline) The fitted pipeline object
#' @param x (IterableMatrix) Input data to be transformed
#' @return Data transformed by the pipeline
#' @noRd
#' @export
setMethod("transform", signature(object = "Pipeline", x = "IterableMatrix"), function(object, x, ...) {
  if (!object@fitted) stop("Pipeline must be fitted before transforming")
  steps <- object@steps
  for (step in steps) {
    if (is(step, "Transformer")) x <- transform(step, x)
    # Some actions convert matrices to a different type, so we need to convert back to IterableMatrix
    # for following steps
    if (is(x, "dgCMatrix")) x <- as(x, "IterableMatrix")
  }
  return(x)
})


#' Predict the output data using a fitted pipeline
#' @param object (Pipeline) The fitted pipeline object
#' @param x (IterableMatrix) Input data to be predicted
#' @return Predicted output data
#' @noRd
#' @export
setMethod("predict", signature(object = "Pipeline", x = "IterableMatrix"), function(object, x, ...) {
  if (!object@fitted) stop("Pipeline must be fitted before predicting")
  steps <- object@steps
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    if (i < n_steps) {
      x <- transform(step, x)
    } else if (is(step, "Estimator")) {
        y_pred <- predict(step, x)
        return(y_pred)
    } else {
      stop("The final step must be an estimator with a predict method")
    }
  }
})

setMethod("short_description", "Pipeline", function(x) {
  character(0)
})

#' Show the pipeline steps, demonstrating how to recreate the pipeline with a function call.
#' @param object (Pipeline) The pipeline object to show
#' @noRd
#' @export
setMethod("show", signature(object = "Pipeline"), function(object) {
  fitted <- ifelse(object@fitted, "Fitted", "Unfitted")
  cat(fitted, " Pipeline with steps:\n")
  cat("Pipeline(\n")
  for (i in seq_along(object@steps)) {
    step <- object@steps[[i]]
    cat("\t", short_description(step))
    if (i < length(object@steps)) {
      cat(",")
    }
    cat("\n")
  }
  cat(")\n")
})

#' Add steps to a pipeline, where the first argument is the pipeline object and the rest are the steps to add.
#' Requires for every additional step to be a pipeline object
#' @param x (Pipeline) The PipelineBase object to add steps to
#' @param ... (PipelineBase) The steps to add to the pipeline
#' @noRd
#' @export
setMethod("c", signature(x = "Pipeline"), function(x, ...) {
  pipelines <- list(...)
  steps <- x@steps
  for (pipe in pipelines) {
    assert_is(pipe, "PipelineBase")
    # If the step is a pipeline, combine the steps.  Else, add the single step.
    steps <- ifelse(is(pipe, "PipelineStep"), c(steps, pipe), c(steps, pipe@steps))
  }

  # If all the steps are fitted, the pipeline overall is fitted.
  # We trust the user to have fitted the pipelines with the same data
  new_pipeline <- Pipeline(steps = steps)
  fitted <- TRUE
  for (step in steps) {
    if (!step@fitted) {
      fitted <- FALSE
    }
  }
  new_pipeline@fitted <- fitted
  return(new_pipeline)
})

#' S4 Class Representing a single transformer or predictor
#' @slot step_name (character) Name of the step
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @name PipelineStep
#' @export
setClass("PipelineStep",
  contains = "PipelineBase",
  slots = list(
    step_name = "character"
  ),
  prototype = list(
    step_name = ""
  )
)

#' Create a Pipeline out of pipeline steps.
#' @param x (PipelineStep) The initial pipeline step we want to add to the pipeline
#' @param ... (PipelineBase) The additional pipeline steps to add to the pipeline.  These can be either PipelineStep or Pipeline objects.
#' @return A new Pipeline object with the steps added.
#' @noRd 
#' @export
setMethod("c", signature(x = "PipelineStep"), function(x, ...) {
  pipelines <- list(...)
  steps <- list(x)
  for (pipe in pipelines) {
    assert_is(pipe, "PipelineBase")
    steps <- ifelse(is(pipe, "PipelineStep"), c(steps, pipe), c(steps, pipe@steps))
  }
  new_pipeline <- Pipeline(steps = steps)
  fitted <- TRUE
  for (step in steps) {
    if (!step@fitted) fitted <- FALSE
  }
})


setMethod("show", signature(object = "PipelineStep"), function(object) {
  cat(short_description(object))
  cat("\n")
})


#' S4 Class representing an operation that transforms data, and holds fitted parameters
#' @slot step_name (character) Name of the step
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @details Transformers represent single operations (derived from the PipelineStep class) that transform data.
#' They can be fit to data within an IterableMatrix object, which will be used to hold the fitted parameters.
#' Using the transform method on a fitted transformer will apply the transformation to the data and return
#' the transformed data as an IterableMatrix object.
#' 
#' These objects can be combined into a Pipeline object using the `c()` function, with other transformers, or estimators.
#' Transformers can also be combined with full pipelines, to create a new pipeline object.
#' Derived classes should implement the `fit()`, `transform()`, and `short_description()` methods.
#' @name Transformer
#' @export
setClass("Transformer",
  contains = "PipelineStep"
)


#' Perform latent semantic indexing (LSI) on a matrix.
#' @name LSITransformer
#' @export
setClass("LSITransformer",
  contains = "Transformer",
  slots = list(
    idf_ = "numeric",
    svd_attr_ = "list",
    z_score_norm = "logical",
    n_dimensions = "integer",
    scale_factor = "integer",
    threads = "integer"
  ),
  prototype = list(
    idf_ = numeric(0),
    svd_attr_ = list(),
    z_score_norm = FALSE,
    n_dimensions = 20L,
    scale_factor = 1e4L,
    threads = 1L
  )
)

#' Create a new LSITransformer object
#' @export
LSITransformer <- function(z_score_norm, n_dimensions, scale_factor, threads) {
  return(new(
    "LSITransformer", z_score_norm = z_score_norm, n_dimensions = n_dimensions, 
    scale_factor = scale_factor, threads = threads, step_name = "LSITransformer"))
}

setMethod("fit", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
  ret <- lsi(
    x, z_score_norm = object@z_score_norm, n_dimensions = object@n_dimensions, 
    scale_factor = object@scale_factor, threads = object@threads,
    save_lsi = TRUE
  )
  object@idf_ <- ret$idf
  object@svd_attr_ <- ret$svd_attr
  object@fitted <- TRUE
  return(object)
})

setMethod("transform", signature(object = "LSITransformer", x = "IterableMatrix"), function(object, x, ...) {
  # rudimentary implementation -- Works but is duplicate code.  
  assert_true(object@fitted)
  # Wait until LSI PR has been reviewed
  npeaks <- colSums(x) # Finding that sums are non-multithreaded and there's no interface to pass it in, but there is implementation in `ConcatenateMatrix.h`
  tf <- x %>% multiply_cols(1 / npeaks)
  mat_tfidf <- tf %>% multiply_rows(object@idf_)
  mat_log_tfidf <- log1p(object@scale_factor * mat_tfidf)
  mat_log_tfidf <- write_matrix_dir(mat_log_tfidf, tempfile("mat_log_tfidf"), compress = FALSE)
  if (object@z_score_norm) {
    cell_peak_stats <- matrix_stats(mat_log_tfidf, col_stats = "variance", threads = object@threads)$col_stats
    cell_means <- cell_peak_stats["mean",]
    cell_vars <- cell_peak_stats["variance",]
    mat_log_tfidf <- mat_log_tfidf %>%
      add_cols(-cell_means) %>%
      multiply_cols(1 / cell_vars)
  }
  pca_res <- t(object@svd_attr_$u) %*% mat_log_tfidf
  return(pca_res)
})

setMethod("short_description", "LSITransformer", function(x) {
  return(sprintf("LSITransformer(z_score_norm=%s, n_dimensions=%d, scale_factor=%d, threads=%d)",
    x@z_score_norm, x@n_dimensions, x@scale_factor, x@threads))
})

setClass("VarFeatSelectorTransformer",
  contains = "Transformer",
  slots = list(
    num_feats = "integer",
    n_bins = "integer"
  )
)


#' S4 Class representing an operation that predicts data, and holds fitted parameters.
#' @slot step_name (character) Name of the step
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @name Estimator
#' @export
setClass("Estimator",
  contains = "PipelineStep"
)