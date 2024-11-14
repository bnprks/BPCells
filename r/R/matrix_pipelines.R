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
#' @noRd
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

#' Fit a pipeline object to data
#' @param object (PipelineBase) The pipeline object to fit.
#' @param x (IterableMatrix) Input data to be fitted on.
#' @param y Optional output data to be fitted on.  Required if the final step is a supervised Estimator, else ignored.
#' @return The fitted pipeline object.
#' @details The `fit()` method is used to fit a pipeline object to data and a potential label output.  Within single estimators, the `fit()` method only
#' takes the input data to be fitted on. Within pipelines, the `fit()` method sequentially fits the transformers on each non-terminal step of the pipeline. More specifically,
#' The input data is transformed by each transformer, and used to fit the next transformer in the pipeline.  If the final step is an estimator, the input IterableMatrix
#' and label (if supervised) are used to fit the estimator.
#'
#' The fitted pipeline object is returned, allowing for projection of new data.
#' @name fit(PipelineBase,IterableMatrix)
#' @export
setGeneric("fit", function(object, x, y = NULL, ...) standardGeneric("fit"))
#' @export
setMethod("fit", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, y = NULL, ...) {
  stop("fit() method not implemented for PipelineBase")
})

#' Project input data using a fitted pipeline
#' @param object (PipelineBase) A fitted pipeline object
#' @param x (IterableMatrix) Input data to be transformed
#' @return Data projected by the pipeline
#' @name project(PipelineBase,IterableMatrix)
#' @export
setGeneric("project", function(object, x, ...) standardGeneric("project"))
#' @export
setMethod("project", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, ...) {
  stop("project() method not implemented for PipelineBase")
})

#' Estimate predictions on the output data using a fitted pipeline
#' @param object (PipelineBase) The fitted pipeline object. Either the final step is an Estimator, or the pipeline is a single Estimator.
#' @param x (IterableMatrix) Input data to be estimated on
#' @return Predicted output labels
#' @name estimate(PipelineBase,IterableMatrix)
#' @export
setGeneric("estimate", function(object, x, ...) standardGeneric("estimate"))
#' @export
setMethod("estimate", signature(object = "PipelineBase", x = "IterableMatrix"), function(object, x, ...) {
  stop("estimate() method not implemented for PipelineBase")
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

#' Return a new Pipeline object.
#' @param steps A list of ordered steps to be executed in the pipeline.
#' @return A new Pipeline object.
#' @details Creating a pipeline object can be done by passing a list of pipeline steps to the constructor.  
#' Creation only expects that all steps make logical sense.  i.e., the final step can be either an Estimator or a Transformer, 
#' but each intermediate step cannot be an Estimator.
#' @export
Pipeline <- function(steps = list()) {
  # Check if all steps are transformers, with the final step being either an estimator or a transformer
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    # allow to fit with estimators as well
    if (i < length(steps)) {
      assert_is(step, "Transformer")
    } else {
      assert_is(step, "PipelineStep")
    }
  }
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
      x <- project(step, x)
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


#' Project input data using a fitted pipeline
#' @param object (Pipeline) The fitted pipeline object
#' @param x (IterableMatrix) Input data to be used by the pipeline to project new data.
#' @return Data projected by the pipeline
#' @noRd
#' @export
setMethod("project", signature(object = "Pipeline", x = "IterableMatrix"), function(object, x, ...) {
  if (!object@fitted) stop("Pipeline must be fitted before projecting")
  steps <- object@steps
  for (step in steps) {
    if (is(step, "Transformer")) x <- project(step, x)
    # Some actions convert matrices to a different type, so we need to convert back to IterableMatrix
    # for following steps
    if (is(x, "dgCMatrix")) x <- as(x, "IterableMatrix")
  }
  return(x)
})


# #' Estimate predictions on the output data using a fitted pipeline
# #' @param object (Pipeline) The fitted pipeline object
# #' @param x (IterableMatrix) Input data to be estimated on
# #' @noRd
# #' @export
# setMethod("estimate", signature(object = "Pipeline", x = "IterableMatrix"), function(object, x, ...) {
#   if (!object@fitted) stop("Pipeline must be fitted before estimating")
#   steps <- object@steps
#   for (i in seq_along(steps)) {
#     step <- steps[[i]]
#     if (i < n_steps) {
#       x <- project(step, x)
#     } else if (is(step, "Estimator")) {
#         y_pred <- estimate(step, x)
#         return(y_pred)
#     } else {
#       stop("The final step must be an estimator with a estimate method")
#     }
#   }
# })

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
    # If the step is a pipeline step, add the single step.  Else, the step is a full pipeline and we want to move all the steps over.
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

#' PipelineBase representing a single step within a pipeline
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
    step_name = character(0)
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
  new_pipeline@fitted <- fitted
  return(new_pipeline)
})


setMethod("show", signature(object = "PipelineStep"), function(object) {
  cat(short_description(object))
  cat("\n")
})


#' PipelineStep representing an operation that transforms data, and holds fitted parameters
#' @slot step_name (character) Name of the step
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @details Transformers represent single operations (derived from the PipelineStep class) that project data from an IterableMatrix to another IterableMatrix/dgCMatrix.
#' They can be fit to data within an IterableMatrix object, which will be used to hold the fitted parameters.
#' Using the `project()`` method on a fitted transformer will apply the transformation to the data and return
#' the projected data as an IterableMatrix object, or a dgCMatrix.
#'
#' These objects can be combined into a Pipeline object using the `c()` function, with other transformers, or estimators.
#' Transformers can also be combined with full pipelines, to create a new pipeline object.
#' Derived classes should implement the `fit()`, `project()`, and `short_description()` methods.
#' @name Transformer
#' @export
setClass("Transformer",
  contains = "PipelineStep"
)


#' PipelineStep representing an operation that estimates data, and holds fitted parameters.
#' @slot step_name (character) Name of the step
#' @slot fitted (logical) Whether the pipeline has been fitted
#' @details Estimators represent single operations (derived from the PipelineStep class) that make predictions based on data given by an
#' IterableMatrix.  Additionally, supervised estimators will require a target numeric, or character array to be provided during a `fit()` call.
#' Unsupervised estimators, on the other hand, do not require a target array.  Following a `fit()` call, the estimator will hold the fitted parameters
#' such that data can be labeled using the `estimate()` method.
#'
#' Estimators can be combined into a Pipeline object using the `c()` function, with other transformers.  Estimators are required to be the terminal step
#' within a pipeline.  Estimators can also be combined with full pipelines, to create a new pipeline object.
#' Derived classes should implement the `fit()`, `estimate()`, and `short_description()` methods.
#' @name Estimator
#' @export
setClass("Estimator",
  contains = "PipelineStep"
)