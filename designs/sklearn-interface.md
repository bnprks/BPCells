# scikit-learn like S3 interface for fitting models and transforms to Iterable matrices

Linkage of data sets typically require a transform function, particularly a means to convert one data type into a latent shared space.  However, accomplishing this within scATAC analysis is difficult without tools that can allow for a transform/model to be fit within one dataset, and project other datasets using the same transform.  Building an interface that can handle transform/model weights and parameters is crucial to do tasks such as projection of multiple datasets using a singular latent semantic indexing tool, or for classification tasks such as logistic gene transfer.  

As we are currently also looking to develop a iterative latent semantic indexing (Iterative LSI) tool, we look towards creating it with a proof-of-concept scikit-learn like interface that allows for the following:

- Each transform should be able to either be standalone, or linked together with another transform, where each transform takes in and returns a two dimensional matrix.
- An estimator should be able to perform standalone, or be the last step following transforms of matrices.  The estimator should take in a two-dimensional matrix, and return a one dimensional label output for each sample.
- Transforms and matrices should be able to be linked together by inheriting from the same base `Pipeline` class, and concatenating each object together.
- The user should be able to call `fit()` on an iterable matrix, and a transform or estimator.  This would allow the transformers and estimators to change their internal weights and parameters to fit against one specific dataset.  
- Calling `fit()` should fit every transform/estimator within a pipeline.
- Calling `fit()` should take in labels if an estimator exists within a pipeline
- The user should be able to be able to call `transform()`, which would allow any dataset with the same features to be transformed using the same transformer, after being fit.
- The user should be able to call `predict()` after performing a `fit()`, which would perform all transformations included in a pipeline, and return a one dimensional vector of labels.

## Base class
`Pipeline`
* Represents a sequence of steps (transformers or estimators) to be applied to data.
* Fields
    * `steps`: An ordered list of pipeline components.
* Methods
    * `fit(pipeline, X, y = NULL, ...)`: Fits each component in the pipeline sequentially to X.  y is optional, and is only required if the last component is an estimator.
    * `transform(pipeline, X)`: Transforms X by applying each transformer in a sequence.  Returns a transformed matrix.
    * `predict(pipeline, X)`: Transforms X by applying each transformer in a sequence, and then applies the final estimator.  Returns a one dimensional vector of labels.
    * `show(pipeline)`: Prints out the steps in the pipeline, and the parameters of each step.  Print outa s a function call for recreating the pipeline.
* Operators
    * `c(c1, c2)`: Concatenates two pipelines together, where the output of `c1` is the input of `c2`.

## TF-IDF transformer
`TfidfTransformer`
* Inherits from Pipeline
* Transforms a matrix of counts into a matrix of term frequency-inverse document frequency (TF-IDF) values.
* Fields
    * `idf`: A vector of inverse document frequencies for each term.
    * `tf`: Term frequency normalization factors
* Methods
    * `fit(X, y = None)`: Calculates the inverse document frequencies for each term in X.
    * `transform(X)`: Transforms X into a matrix of TF-IDF values.


## Z score transformer
`ZScoreScaler`
* Inherits from Pipeline
* Standardizes each column of a matrix to have a mean of 0 and a standard deviation of 1.
* Fields
    * `means`: A vector of means for each column.
    * `stds`: A vector of standard deviations for each column.
* Methods
    * `fit(X, y = None)`: Calculates the means and standard deviations for each column of X.
    * `transform(X)`: Standardizes each column of X to have a mean of 0 and a standard deviation of 1.

## SVD transformer
`SVD`
* Inherits from Pipeline
* Performs singular value decomposition on a matrix, and returns the first `n` projected components.
* Fields
    * `n_components`: The number of components to return
    * `components`: Principal components computed during fitting
* Methods
    * `fit(X, y = None)`: Performs singular value decomposition on X, and stores the first `n` components.
    * `transform(X, y = None)`: Projects X onto the first `n` components.

## Feature selector transformer
`FeatureSelector`
* Inherits from Pipeline
* Selects a subset of features from a matrix, using scanpy default approach with binning.    Use their implementation of highly variable genes.
* Fields
    * `features`: A vector of indices of features to select
* Methods
    * `fit(X, y = None)`: Stores the indices of features to select.
    * `transform(X)`: Selects the features from X.


# Workflow
1. Create an abstract base class `Pipeline` that represents a sequence of steps to be applied to data.
2. Create a class `TfidfTransformer` that inherits from `Pipeline` and transforms a matrix of counts into a matrix of TF-IDF values.
3. Create a class `ZScoreScaler` that inherits from `Pipeline` and standardizes each column of a matrix to have a mean of 0 and a standard deviation of 1.  Ensure that tests exemplify the ability to pipeline this with the `TfidfTransformer`.
4. Create a class `SVD` that inherits from `Pipeline` and performs singular value decomposition on a matrix, and returns the first `n` projected components.
5. Create a class `FeatureSelector` that inherits from `Pipeline` and selects a subset of features from a matrix.
6. Create a pipeline wrapper around all aforementioned classes to perform LSI, and compare results against ArchR in a vignette.

# Considerations
* Each transform/model may need to take in additional metadata parameters, such as the number of features, or the number of samples. We should either consider finding a way to change these after object construction, or to pass them in as arguments to the `fit()` method.
    * As discussed previously, this might necessitate an automatic naming function of each step in a pipeline, similar to how objects are named in scikit-learn pipelines.
* All of the classes described here do not have a `predict()` method as they are not estimators/classifiers.  It might be good to consider adding a simple classifier (maybe a logistic regression) to test the pipeline functionality.


# Conclusion
By creating a scikit-learn like interface, we can use as common interface for future models and transforms within BPCells.  This will allow for the development of tools that can be used to project multiple datasets into a shared latent space, or to perform classification tasks on scATAC data.  Additionally, creating the first proof of concept functions, such as the `TfidfTransformer`, `ZScoreScaler`, `SVD`, and `FeatureSelector` will allow for bundling of these functions into a singular pipeline for usage of latent semantic indexing, allowing for projecting multiple datasets into a shared latent space.

