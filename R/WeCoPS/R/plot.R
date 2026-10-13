# Diagnostic sheet (fiche) of a comparison between a weighting and the default weighting.

plot_fiche <- function(res, file = NULL, title = NULL, width = 12, height = 5) {
  if (!is.null(file)) {
    grDevices::png(file, width = width, height = height, units = "in", res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  kap <- res$kappa
  lab <- vapply(res$groups, function(g) paste0(g$label,
    if (length(g$components) > 1) sprintf(" (%d comp.)", length(g$components)) else ""), character(1))
  tryCatch({
    .draw_fiche(res, kap, lab, title)
    TRUE
  }, error = function(e) {
    if (grepl("margins too large", conditionMessage(e)))
      stop("the plot window is too small for the diagnostic sheet: enlarge it, or use ",
           "plot_fiche(res, file = \"fiche.png\")", call. = FALSE)
    stop(e)
  })
  invisible(res)
}

.draw_fiche <- function(res, kap, lab, title) {
  graphics::layout(matrix(1:2, 1, 2), widths = c(1, 1.15))
  graphics::par(oma = c(0, 0, 2, 0), cex = 0.8)
  graphics::par(mar = c(4, min(12, 1 + 0.45 * max(nchar(lab))), 2.5, 1))
  cols <- ifelse(kap >= res$kappa_min, "#2a78d6", "#eb6834")
  bp <- graphics::barplot(rev(kap), horiz = TRUE, names.arg = rev(lab), las = 1, col = rev(cols),
                          border = NA, xlim = c(0, 1.1), xlab = "presence kappa under the alternative weighting")
  graphics::abline(v = res$kappa_min, lty = 2)
  graphics::text(pmin(rev(kap), 1) + 0.01, bp, sprintf("%.2f", rev(kap)), pos = 4, cex = 0.8, offset = 0)
  graphics::title(main = sprintf("(a) %s; motif %s",
                                 if (res$equivalent) "equivalent" else "interpretation changed",
                                 if (res$same_motif) "identical" else "modified"),
                  adj = 0, cex.main = 0.95, font.main = 1)
  D <- res$parts_change
  K <- nrow(D)
  graphics::par(mar = c(4, min(8, 1 + 0.6 * max(nchar(res$block_names))), 2.5, 1))
  pal <- grDevices::colorRampPalette(c("#c2185b", "#f7f6f3", "#4e9a2a"))(101)
  Dc <- pmin(pmax(D, -0.3), 0.3)
  graphics::image(x = 1:4, y = 1:K, z = t(Dc[K:1, , drop = FALSE]), col = pal, zlim = c(-0.3, 0.3),
                  axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(K)) for (j in 1:4) {
    v <- D[K - i + 1, j]
    if (abs(v) < 0.005) v <- 0
    graphics::text(j, i, if (v == 0) "0.00" else sprintf("%+.2f", v), cex = 0.85,
                   font = if (abs(v) >= 0.05) 2 else 1)
  }
  graphics::axis(1, at = 1:4, labels = c("common", "own proper", "other proper", "residual"),
                 tick = FALSE, cex.axis = 0.8)
  graphics::axis(2, at = 1:K, labels = rev(res$block_names), las = 1, tick = FALSE, cex.axis = 0.85)
  graphics::title(main = "(b) change of the inertia of each block (alternative - default)",
                  adj = 0, cex.main = 0.95, font.main = 1)
  graphics::mtext(if (is.null(title)) sprintf("Weighting '%s' versus default, R = %d",
                                              res$weighting_name, as.integer(res$R)) else title,
                  side = 3, outer = TRUE, line = 0.5, adj = 0.01, font = 2)
  invisible(NULL)
}

plot.wecops_comparison <- function(x, ...) plot_fiche(x, ...)
