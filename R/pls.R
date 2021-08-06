#' Partial Least Squares (PLS)
#'
#' This function implements the generic PLS algorithm
#'
#' @param X Matrix of descriptors (n x p)
#' @param Y Matrix of response variables (n x q)
#' @param num.components Number of PLS components (variates) to compute
#'
#' @return Named list containing weights, scores, loadings, rotation matrices and coefficients
#'
#' @importFrom MASS ginv
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

  svd.flip.sign <- function(u, v) {
    # Function to flip signs of results SVD

    idx.largest.abs <- which.max(abs(u))
    sign.abs <- sign(u[idx.largest.abs])
    u <- u * sign.abs
    v <- v * sign.abs

    return(list(u, v))
  }

  ### VARIABLES ###
  # Number of observations
  n <- dim(X)[1]
  # Number of explanatory variables or descriptors
  p <- dim(X)[2]
  # Number of response variables
  q <- dim(Y)[2]

  # Matrices to store weights, scores and loadings of X and Y
  X.weights.comps  <- matrix(0, nrow = p, ncol = num.components)
  Y.weights.comps  <- matrix(0, nrow = q, ncol = num.components)
  X.scores.comps   <- matrix(0, nrow = n, ncol = num.components)
  Y.scores.comps   <- matrix(0, nrow = n, ncol = num.components)
  X.loadings.comps <- matrix(0, nrow = p, ncol = num.components)
  Y.loadings.comps <- matrix(0, nrow = q, ncol = num.components)

  ### SCALING ###
  Xk <- X
  Yk <- Y

  Y.sd <- sd(Y)
  Y.sd[Y.sd == 0] <- 1.0

  ### ALGORITHM ###
  for (k in 1:num.components) {
    # Find 1st left and right singular vectors of the X.T.dot(Y)
    # cross-covariance matrix
    weights <- get.first.singular.vectors(Xk, Yk)
    X.weights <- weights$X # (p x 1)
    Y.weights <- weights$Y # (q x 1)

    # Flip sign for consistency across solvers and archs
    res <- svd.flip.sign(X.weights, Y.weights)
    X.weights <- res[[1]]
    Y.weights <- res[[2]]

    # Compute scores (i.e., projections of X and Y)
    X.scores <- Xk %*% X.weights # (n x 1)
    Y.ss <- as.numeric(crossprod(Y.weights)) # Equivalent to `t(Y.weights) %*% Y.weights`
    Y.scores <- (Yk %*% Y.weights) / Y.ss # (n x 1)

    # Deflation (i.e., subtract rank-one approximation to obtain
    # X_(k+1) and Y_(k+1))
    X.loadings <- (t(X.scores) %*% Xk) / as.numeric(crossprod(X.scores)) # (1 x p)
    Xk <- Xk - (X.scores %*% X.loadings)

    Y.loadings <- (t(X.scores) %*% Yk) / as.numeric(crossprod(X.scores)) # (1 x q)
    Yk <- Yk - (X.scores %*% Y.loadings)

    # Store weights, scores and loadings for X and Y for current component k
    X.weights.comps[, k]  <- X.weights
    Y.weights.comps[, k]  <- Y.weights
    X.scores.comps[, k]   <- X.scores
    Y.scores.comps[, k]   <- Y.scores
    X.loadings.comps[, k] <- X.loadings
    Y.loadings.comps[, k] <- Y.loadings

  } # End loop for PLS

  # Compute transformation matrices
  X.rotations <- X.weights.comps %*% MASS::ginv(t(X.loadings.comps) %*% X.weights.comps)
  Y.rotations <- Y.weights.comps %*% MASS::ginv(t(Y.loadings.comps) %*% Y.weights.comps)

  coef <- X.rotations %*% t(Y.loadings)
  coef <- coef * Y.sd

  # Create named list to return
  ret <- list(
    X.weights   = X.weights.comps,
    Y.weights   = Y.weights.comps,
    X.scores    = X.scores.comps,
    Y.scores    = Y.scores.comps,
    X.loadings  = X.loadings.comps,
    Y.loadings  = Y.loadings.comps,
    X.rotations = X.rotations,
    Y.rotations = Y.rotations,
    coef = coef
  )
  return(ret)

}
