# One weighted CPCA read in two bases: principal components and common and proper components.

wecops <- function(blocks, R, weights = "default", alpha = 0.05, nstart = 10, seed = 0,
                   block_names = NULL) {
  blocks <- .as_blocks(blocks)
  w <- block_weights(blocks, weights)
  a <- cpca(blocks, w, R)
  U <- a$U
  cc <- ccp(blocks, U, alpha, nstart = nstart, seed = seed, weights = w)
  V <- cc$V
  T <- cc$T
  S_pc <- saliences(U, blocks)
  T_pc <- block_labels(S_pc, alpha)
  S <- saliences(V, blocks)
  bn <- if (is.null(block_names)) paste("block", seq_along(blocks)) else block_names
  parts <- four_parts(S, T)
  rownames(parts) <- bn
  structure(list(
    weights = w, R = R, alpha = alpha, eigenvalues = a$explained,
    U = U, saliences_pc = S_pc, labels_pc = label_names(T_pc, block_names), motif_pc = motif(T_pc, block_names),
    V = V, saliences = S, T = T, labels = label_names(T, block_names), motif = motif(T, block_names),
    mixing = mixing(U, V, T), parts = parts,
    factorization = list(ccp = canonical_factorization(blocks, V), pc = canonical_factorization(blocks, U)),
    block_names = bn
  ), class = "wecops")
}

print.wecops <- function(x, ...) {
  cat(sprintf("WeCoPS analysis (R = %d, alpha = %g)\n", as.integer(x$R), x$alpha))
  cat("  principal components:          ", paste(x$labels_pc, collapse = ", "), "\n")
  cat("  common and proper components:  ", paste(x$labels, collapse = ", "), "\n")
  cat(sprintf("  max mixing of principal components: mu = %.2f, mu* = %.2f\n",
              max(x$mixing$mu), max(x$mixing$mu_star)))
  invisible(x)
}
