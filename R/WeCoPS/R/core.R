# Block preprocessing, block weightings and weighted consensus PCA (CPCA).

.as_blocks <- function(blocks) {
  if (!is.list(blocks) || length(blocks) < 2)
    stop("'blocks' must be a list of at least two matrices")
  out <- lapply(blocks, function(X) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X
  })
  n <- nrow(out[[1]])
  for (k in seq_along(out)) {
    if (nrow(out[[k]]) != n)
      stop("all blocks must have the same number of rows (individuals)")
    if (any(!is.finite(out[[k]])))
      stop(sprintf("block %d contains missing or infinite values", k))
  }
  out
}

snv <- function(X) {
  X <- as.matrix(X)
  m <- rowMeans(X)
  s <- apply(X, 1, stats::sd)
  (X - m) / s
}

prep_blocks <- function(blocks, center = TRUE, unit_inertia = TRUE) {
  lapply(.as_blocks(blocks), function(X) {
    if (center) X <- sweep(X, 2, colMeans(X))
    if (unit_inertia) {
      nrm <- sqrt(sum(X^2))
      if (nrm == 0) stop("a block has zero inertia")
      X <- X / nrm
    }
    X
  })
}

.statis <- function(blocks) {
  W <- lapply(blocks, function(X) X %*% t(X))
  nW <- vapply(W, function(Wk) sqrt(sum(Wk^2)), numeric(1))
  Wt <- Map(function(Wk, nk) Wk / nk, W, nW)
  K <- length(Wt)
  S <- matrix(0, K, K)
  for (a in seq_len(K)) for (b in seq_len(K)) S[a, b] <- sum(Wt[[a]] * Wt[[b]])
  beta <- abs(eigen(S, symmetric = TRUE)$vectors[, 1])
  list(weights = beta / nW, beta = beta / sum(beta), rv = S)
}

block_weights <- function(blocks, method = "default") {
  blocks <- .as_blocks(blocks)
  K <- length(blocks)
  if (!is.character(method)) {
    w <- as.numeric(method)
    if (length(w) != K || any(w <= 0)) stop("weights must be K positive numbers")
    return(w)
  }
  switch(tolower(method),
    default = rep(1, K),
    mfa = vapply(blocks, function(X) 1 / svd(X, nu = 0, nv = 0)$d[1]^2, numeric(1)),
    rv = vapply(blocks, function(X) 1 / sqrt(sum((X %*% t(X))^2)), numeric(1)),
    statis = .statis(blocks)$weights,
    nvar = vapply(blocks, function(X) 1 / ncol(X), numeric(1)),
    stop(sprintf("unknown weighting '%s'; use 'default', 'rv', 'mfa', 'statis', 'nvar' or a vector",
                 method))
  )
}

cpca <- function(blocks, weights = NULL, R = 2) {
  blocks <- .as_blocks(blocks)
  K <- length(blocks)
  w <- if (is.null(weights)) rep(1, K) else block_weights(blocks, weights)
  wn <- w / sum(w)
  n <- nrow(blocks[[1]])
  M <- matrix(0, n, n)
  for (k in seq_len(K)) M <- M + wn[k] * (blocks[[k]] %*% t(blocks[[k]]))
  e <- eigen(M, symmetric = TRUE)
  lam <- pmax(e$values, 0)
  if (R < 1 || R >= n) stop("R must be between 1 and n - 1")
  U <- e$vectors[, seq_len(R), drop = FALSE]
  list(U = U, eigenvalues = lam, explained = lam / sum(lam), weights = w)
}

saliences <- function(U, blocks) {
  U <- as.matrix(U)
  do.call(rbind, lapply(blocks, function(X) colSums((t(X) %*% U)^2)))
}
