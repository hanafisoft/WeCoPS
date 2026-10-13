# Common and proper components: concentration (INDORT criterion restricted to the retained
# subspace), data-driven target at threshold alpha, and rotation towards the target.
# All rotations act in the retained subspace: V = U Q, s_kh = q_h' (U' W_k U) q_h.

.projected <- function(U, blocks) {
  lapply(blocks, function(X) {
    L <- t(X) %*% U
    t(L) %*% L
  })
}

.sal <- function(Q, Wt) {
  do.call(rbind, lapply(Wt, function(Wk) colSums(Q * (Wk %*% Q))))
}

concentration <- function(Wt, nstart = 10, maxit = 1000, tol = 1e-12, seed = 0) {
  R <- nrow(Wt[[1]])
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  best_f <- -Inf
  best_Q <- NULL
  for (st in seq_len(nstart)) {
    Q <- if (st == 1) diag(R) else qr.Q(qr(matrix(stats::rnorm(R * R), R, R)))
    for (it in seq_len(maxit)) {
      S <- .sal(Q, Wt)
      G <- matrix(0, R, R)
      for (k in seq_along(Wt)) G <- G + sweep(Wt[[k]] %*% Q, 2, S[k, ], `*`)
      s <- svd(G)
      Qn <- s$u %*% t(s$v)
      if (max(abs(Qn - Q)) < tol) {
        Q <- Qn
        break
      }
      Q <- Qn
    }
    f <- sum(.sal(Q, Wt)^2)
    if (f > best_f + 1e-12) {
      best_f <- f
      best_Q <- Q
    }
  }
  best_Q
}

target_rotation <- function(Q, Wt, T, sweeps = 100, tol = 1e-12) {
  R <- ncol(Q)
  if (R < 2) return(Q)
  for (sw in seq_len(sweeps)) {
    moved <- 0
    for (i in seq_len(R - 1)) {
      for (j in (i + 1):R) {
        qi <- Q[, i]
        qj <- Q[, j]
        B <- 0
        C <- 0
        for (k in seq_along(Wt)) {
          zi <- T[k, i] == 0
          zj <- T[k, j] == 0
          if (zi == zj) next  # both or none excluded: the criterion does not depend on the angle
          a <- sum(qi * (Wt[[k]] %*% qi)); b <- sum(qj * (Wt[[k]] %*% qj)); d <- sum(qi * (Wt[[k]] %*% qj))
          sg <- if (zi) 1 else -1
          B <- B + sg * (a - b) / 2
          C <- C + sg * d
        }
        if (B * B + C * C < tol * tol) next
        t0 <- 0.5 * atan2(-C, -B)
        Q[, i] <- cos(t0) * qi + sin(t0) * qj
        Q[, j] <- -sin(t0) * qi + cos(t0) * qj
        moved <- max(moved, abs(t0))
      }
    }
    if (moved < 1e-10) break
  }
  Q
}

block_labels <- function(S, alpha = 0.05, previous = NULL, band = NULL) {
  S <- as.matrix(S)
  if (is.null(previous) || is.null(band)) {
    T <- (S >= alpha) * 1L
  } else {
    T <- previous
    T[S >= band[2]] <- 1L
    T[S < band[1]] <- 0L
  }
  storage.mode(T) <- "integer"
  T
}

label_names <- function(T, block_names = NULL) {
  K <- nrow(T)
  nm <- if (is.null(block_names)) as.character(seq_len(K)) else block_names
  vapply(seq_len(ncol(T)), function(h) {
    b <- nm[T[, h] == 1]
    if (length(b) > 1) paste("common", paste(b, collapse = "-"))
    else if (length(b) == 1) paste("proper", b)
    else "empty"
  }, character(1))
}

motif <- function(T, block_names = NULL) {
  sort(label_names(T, block_names), method = "radix")
}

label_groups <- function(T) {
  key <- apply(T, 2, paste, collapse = "")
  lapply(unique(key), function(L) which(key == L))
}

orient_groups <- function(Q, Wt, T, U, weights = NULL) {
  K <- length(Wt)
  w <- if (is.null(weights)) rep(1, K) else as.numeric(weights)
  for (idx in label_groups(T)) {
    if (length(idx) < 2) next
    ks <- which(T[, idx[1]] == 1)
    if (length(ks) == 0) ks <- seq_len(K)
    QL <- Q[, idx, drop = FALSE]
    M <- matrix(0, length(idx), length(idx))
    for (k in ks) M <- M + w[k] * (t(QL) %*% Wt[[k]] %*% QL)
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    Q[, idx] <- QL %*% e$vectors
  }
  V <- U %*% Q
  sgn <- sign(V[cbind(apply(abs(V), 2, which.max), seq_len(ncol(V)))])
  sgn[sgn == 0] <- 1
  sweep(Q, 2, sgn, `*`)
}

ccp <- function(blocks, U, alpha = 0.05, nstart = 10, seed = 0, maxit = 20, sweeps = 100,
                weights = NULL) {
  Wt <- .projected(U, blocks)
  Q <- concentration(Wt, nstart = nstart, seed = seed)
  T <- block_labels(.sal(Q, Wt), alpha)
  for (it in seq_len(maxit)) {
    Q <- target_rotation(Q, Wt, T, sweeps = sweeps)
    T2 <- block_labels(.sal(Q, Wt), alpha)
    if (identical(T2, T)) break
    T <- T2
  }
  Q <- orient_groups(Q, Wt, T, U, weights)
  S <- .sal(Q, Wt)
  T <- block_labels(S, alpha)  # labels of the final (oriented) components
  list(V = U %*% Q, Q = Q, T = T, S = S)
}

.ccp_from <- function(blocks, U, Q0, T_prev, band, maxit = 20, sweeps = 100, weights = NULL) {
  Wt <- .projected(U, blocks)
  Q <- Q0
  T <- block_labels(.sal(Q, Wt), previous = T_prev, band = band)
  for (it in seq_len(maxit)) {
    Q <- target_rotation(Q, Wt, T, sweeps = sweeps)
    T2 <- block_labels(.sal(Q, Wt), previous = T_prev, band = band)
    if (identical(T2, T)) break
    T <- T2
  }
  Q <- orient_groups(Q, Wt, T, U, weights)
  S <- .sal(Q, Wt)
  T <- block_labels(S, previous = T_prev, band = band)  # labels of the final (oriented) components
  list(V = U %*% Q, Q = Q, T = T, S = S)
}
