#' Partial Least Squares (PLS)
#'
#' This function implements the generic PLS algorithm
#'
#' @param X Matrix of descriptors (n x p)
#' @param Y Matrix of response variables (n x q)
#' @param num.components Number of PLS components (variates) to compute
#'
#' @return None
#'
#' @export

pls <- function(X, Y, num.components) {

  ### HELPER FUNCTIONS ###
  get.first.singular.vectors <- function(X, Y) {
    # Returns the 1st left and right singular vectors
    # of X.T.dot(Y)

    C = t(X) %*% Y
    # Compute only one component for each resulting matrix
    svd.res <- svd(C, nu = 1, nv = 1)
    U <- svd.res$u
    V <- svd.res$v

    weights <- list(X = U, Y = V)
    return(weights)
  }

  ### VARIABLES ###
  # Number of observations
  n <- dim(X)[1]
  # Number of explanatory variables or descriptors
  p <- dim(X)[2]
  # Number of response variables
  q <- dim(Y)[2]

  ### SCALING ###

  ### ALGORITHM ###
  for (iter in 1:num.components) {
    # Find 1st left and right singular vectors of the X.T.dot(Y)
    # cross-covariance matrix
    weights <- get.first.singular.vectors(X, Y)
    X.weights <- weights$X
    Y.weights <- weights$Y

    # Flip sign for consistency across solvers and archs

  }

}
