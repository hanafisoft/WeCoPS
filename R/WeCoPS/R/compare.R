# Mixing of principal components, presence, four-part decomposition, canonical factorisation
# and the equivalence rule between a weighting and the default weighting.

mixing <- function(U, V, T) {
  G <- t(V) %*% U
  G2 <- G^2
  mu <- 1 - apply(G2, 2, max)
  key <- apply(T, 2, paste, collapse = "")
  pooled <- do.call(rbind, lapply(unique(key), function(L) colSums(G2[key == L, , drop = FALSE])))
  mu_star <- 1 - apply(pooled, 2, max)
  list(gamma = G, mu = pmax(mu, 0), mu_star = pmax(mu_star, 0))
}

presence <- function(U_alt, V0) {
  colSums((t(U_alt) %*% V0)^2)
}

group_presence <- function(U_alt, V0, T0, block_names = NULL) {
  P <- t(U_alt) %*% V0
  lapply(label_groups(T0), function(idx) {
    G <- t(P[, idx, drop = FALSE]) %*% P[, idx, drop = FALSE]
    list(label = label_names(T0[, idx[1], drop = FALSE], block_names), components = idx,
         kappa = min(eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values))
  })
}

four_parts <- function(S, T) {
  K <- nrow(S)
  nb <- colSums(T)
  P <- matrix(0, K, 4, dimnames = list(NULL, c("common", "own_proper", "other_proper", "residual")))
  for (k in seq_len(K)) {
    P[k, 1] <- sum(S[k, nb >= 2])
    P[k, 2] <- sum(S[k, nb == 1 & T[k, ] == 1])
    P[k, 3] <- sum(S[k, nb < 2 & !(nb == 1 & T[k, ] == 1)])
    P[k, 4] <- 1 - sum(S[k, ])
  }
  P
}

canonical_factorization <- function(blocks, B) {
  lapply(blocks, function(X) {
    L <- t(X) %*% B
    d <- sqrt(colSums(L^2))
    Vk <- sweep(L, 2, ifelse(d > 0, d, 1), "/")
    Vk[, d == 0] <- 0
    list(D = diag(d, nrow = length(d)), V = Vk)
  })
}

reweighting_bound <- function(blocks, R, weights, V0 = NULL) {
  blocks <- .as_blocks(blocks)
  w <- block_weights(blocks, weights)
  n <- nrow(blocks[[1]])
  M1 <- matrix(0, n, n)
  for (X in blocks) M1 <- M1 + X %*% t(X)
  lam <- pmax(eigen(M1, symmetric = TRUE, only.values = TRUE)$values, 0)
  ratio <- if (lam[R + 1] > 0) lam[R] / lam[R + 1] else Inf
  A <- max(w) / min(w)
  out <- list(A = A, ratio = ratio, guaranteed = A < ratio)
  if (!is.null(V0)) {
    rho <- colSums(V0 * (M1 %*% V0))
    out$kappa_lower_bound <- if (out$guaranteed)
      (rho / A - lam[R + 1]) / (lam[1] - lam[R + 1]) else NULL
  }
  out
}

# Hungarian algorithm (minimum-cost assignment of a square matrix); returns, for each row,
# the assigned column.
.hungarian <- function(a) {
  n <- nrow(a)
  u <- numeric(n + 1); v <- numeric(n + 1)
  p <- integer(n + 1); way <- integer(n + 1)
  for (i in seq_len(n)) {
    p[1] <- i
    j0 <- 0L
    minv <- rep(Inf, n + 1)
    used <- rep(FALSE, n + 1)
    repeat {
      used[j0 + 1] <- TRUE
      i0 <- p[j0 + 1]
      delta <- Inf
      j1 <- 0L
      for (j in seq_len(n)) {
        if (!used[j + 1]) {
          cur <- a[i0, j] - u[i0 + 1] - v[j + 1]
          if (cur < minv[j + 1]) {
            minv[j + 1] <- cur
            way[j + 1] <- j0
          }
          if (minv[j + 1] < delta) {
            delta <- minv[j + 1]
            j1 <- j
          }
        }
      }
      for (j in 0:n) {
        if (used[j + 1]) {
          u[p[j + 1] + 1] <- u[p[j + 1] + 1] + delta
          v[j + 1] <- v[j + 1] - delta
        } else {
          minv[j + 1] <- minv[j + 1] - delta
        }
      }
      j0 <- j1
      if (p[j0 + 1] == 0) break
    }
    repeat {
      j1 <- way[j0 + 1]
      p[j0 + 1] <- p[j1 + 1]
      j0 <- j1
      if (j0 == 0) break
    }
  }
  ans <- integer(n)
  for (j in seq_len(n)) ans[p[j + 1]] <- j
  ans
}

.match <- function(U0, U) {
  C <- abs(t(U0) %*% U)
  col <- .hungarian(-C)
  Um <- U[, col, drop = FALSE]
  Um <- sweep(Um, 2, sign(colSums(U0 * Um)), `*`)
  list(U = Um, cos2 = C[cbind(seq_len(nrow(C)), col)]^2)
}

.procrustes <- function(A) {
  s <- svd(A)
  s$u %*% t(s$v)
}

compare_weighting <- function(blocks, R, weights, alpha = 0.05, kappa_min = 0.9, band = 0.02,
                              nstart = 10, seed = 0, block_names = NULL) {
  blocks <- .as_blocks(blocks)
  K <- length(blocks)
  w <- block_weights(blocks, weights)
  name <- if (is.character(weights)) weights else "user"
  bnd <- if (band > 0) c(max(alpha - band, alpha / 2), alpha + band) else c(alpha, alpha)
  if (is.null(block_names)) block_names_out <- paste("block", seq_len(K)) else block_names_out <- block_names

  a0 <- cpca(blocks, NULL, R)
  U0 <- a0$U
  c0 <- ccp(blocks, U0, alpha, nstart = nstart, seed = seed)
  V0 <- c0$V
  T0 <- c0$T

  ag <- cpca(blocks, w, R)
  Ug <- ag$U
  cg <- .ccp_from(blocks, Ug, .procrustes(t(Ug) %*% V0), T0, bnd, weights = w)
  Vg <- cg$V
  Tg <- cg$T

  kap_c <- presence(Ug, V0)
  groups <- group_presence(Ug, V0, T0, block_names)
  kap <- vapply(groups, function(g) g$kappa, numeric(1))
  m0 <- motif(T0, block_names)
  mg <- motif(Tg, block_names)
  same <- identical(m0, mg)
  equivalent <- same && min(kap) >= kappa_min

  mt <- .match(U0, Ug)
  S_pc0 <- saliences(U0, blocks)
  S_pcg <- saliences(mt$U, blocks)
  pc_same <- min(mt$cos2) >= kappa_min && identical(S_pc0 >= alpha, S_pcg >= alpha)

  S0 <- saliences(V0, blocks)
  Sg <- saliences(Vg, blocks)
  p0 <- four_parts(S0, T0)
  pg <- four_parts(Sg, Tg)
  rownames(p0) <- rownames(pg) <- block_names_out

  structure(list(
    weighting_name = name, weights = w, R = R, alpha = alpha, kappa_min = kappa_min, band = bnd,
    equivalent = equivalent, same_motif = same, kappa = kap, groups = groups, kappa_components = kap_c,
    motif_default = m0, motif_alt = mg,
    labels_default = label_names(T0, block_names), labels_alt = label_names(Tg, block_names),
    T_default = T0, T_alt = Tg, V_default = V0, V_alt = Vg, U_default = U0, U_alt = Ug,
    saliences_default = S0, saliences_alt = Sg,
    parts_default = p0, parts_alt = pg, parts_change = pg - p0,
    eigenvalues_default = a0$explained, eigenvalues_alt = ag$explained,
    mixing_default = mixing(U0, V0, T0), mixing_alt = mixing(Ug, Vg, Tg),
    pc_same = pc_same, pc_cos2 = mt$cos2,
    bound = reweighting_bound(blocks, R, w, V0),
    factorization = list(
      default_ccp = canonical_factorization(blocks, V0), default_pc = canonical_factorization(blocks, U0),
      alt_ccp = canonical_factorization(blocks, Vg), alt_pc = canonical_factorization(blocks, Ug)),
    block_names = block_names_out
  ), class = "wecops_comparison")
}

summary.wecops_comparison <- function(object, ...) {
  x <- object
  lines <- c(
    sprintf("WeCoPS comparison: weighting '%s' versus default (R = %d, alpha = %g)",
            x$weighting_name, as.integer(x$R), x$alpha),
    sprintf("  verdict: %s (same motif: %s, min presence: %.3f, threshold %g)",
            if (x$equivalent) "EQUIVALENT" else "INTERPRETATION CHANGED", x$same_motif,
            min(x$kappa), x$kappa_min),
    paste0("  default motif:     ", paste(x$motif_default, collapse = ", ")),
    paste0("  alternative motif: ", paste(x$motif_alt, collapse = ", ")),
    paste0("  presence of each labelled space of the default: ",
           paste(vapply(x$groups, function(g) sprintf("%s%s: %.2f", g$label,
             if (length(g$components) > 1) sprintf(" (%d comp.)", length(g$components)) else "",
             g$kappa), character(1)), collapse = ", ")),
    sprintf("  weight ratio A = %.3g; lambda_R/lambda_R+1 = %.3g%s", x$bound$A, x$bound$ratio,
            if (x$bound$guaranteed) "  -> below the reweighting threshold: no component can leave the subspace entirely" else ""),
    sprintf("  principal-component reading: %s (max mixing mu* default %.2f, alternative %.2f)",
            if (x$pc_same) "same" else "changed", max(x$mixing_default$mu_star), max(x$mixing_alt$mu_star))
  )
  if (x$equivalent && !x$pc_same)
    lines <- c(lines, "  note: reading the principal components one by one would wrongly signal a change")
  structure(lines, class = "wecops_summary")
}

print.wecops_comparison <- function(x, ...) {
  cat(summary.wecops_comparison(x), sep = "\n")
  invisible(x)
}

print.wecops_summary <- function(x, ...) {
  cat(unclass(x), sep = "\n")
  invisible(x)
}
