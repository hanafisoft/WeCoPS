# Vérification sous R des résultats de l'article 2a (WeCoPS 0.2.1).
# Valeurs attendues : calcul Python de référence (WeCoPS 0.2.1). Chaque ligne affiche OK ou PROBLEM.
# Durée : quelques minutes (partie 3 : 196 analyses ; partie 4 : 224 analyses).
library(WeCoPS)
n_ok <- 0; n_pb <- 0
check <- function(label, ok) {
  if (isTRUE(ok)) { n_ok <<- n_ok + 1; cat("OK      ", label, "\n") }
  else { n_pb <<- n_pb + 1; cat("PROBLEM ", label, "\n") }
}
near <- function(a, b, tol = 2e-3) isTRUE(all(abs(a - b) <= tol))
W <- c("rv", "mfa", "statis", "nvar")
D <- list(S1 = load_simulated(1)$blocks, S2 = load_simulated(2)$blocks, S3 = load_simulated(3)$blocks,
          S4 = load_simulated(4)$blocks, vins = load_wines()$blocks, beurres = load_butters()$blocks)
cases <- list(c("S1", 5), c("S2", 10), c("S3", 5), c("S4", 5), c("vins", 3), c("beurres", 4), c("beurres", 3))

cat("\n=== 1. Tableau 2 : max mu* (et nombre de composantes principales mélangées) ===\n")
exp_t2 <- list(
  list("S1", 5, c(0.1540, 0.0102, 0.0085, 0.0102, 0.1540), c(2, 0, 0, 0, 2)),
  list("S2", 10, c(0.5565, 0.4181, 0.5185, 0.4853, 0.5162), c(4, 4, 5, 4, 7)),
  list("S3", 5, c(0.2258, 0.1982, 0.1896, 0.2429, 0.2258), c(2, 2, 2, 4, 2)),
  list("S4", 5, c(0.0250, 0.0814, 0.1049, 0.0416, 0.0250), c(0, 0, 1, 0, 0)),
  list("vins", 3, c(0.0136, 0.0203, 0.0365, 0.0079, 0.0285), c(0, 0, 0, 0, 0)),
  list("beurres", 4, c(0.2858, 0.1760, 0.0943, 0.2648, 0.4324), c(3, 1, 0, 3, 3)),
  list("beurres", 3, c(0.2990, 0.3859, 0.0381, 0.3208, 0.1478), c(2, 2, 0, 2, 2)))
for (e in exp_t2) {
  b <- D[[e[[1]]]]; R <- e[[2]]; mus <- numeric(5); nm <- integer(5)
  for (i in seq_along(W)) {
    r <- compare_weighting(b, R, W[i])
    if (i == 1) { mus[1] <- max(r$mixing_default$mu_star); nm[1] <- sum(r$mixing_default$mu_star > 0.1) }
    mus[i + 1] <- max(r$mixing_alt$mu_star); nm[i + 1] <- sum(r$mixing_alt$mu_star > 0.1)
  }
  check(sprintf("%-8s R=%-2d mu* = %s (attendu %s) ; mélangées = %s", e[[1]], R,
                paste(sprintf("%.3f", mus), collapse = " "), paste(sprintf("%.3f", e[[3]]), collapse = " "),
                paste(nm, collapse = " ")), near(mus, e[[3]], 5e-3) && all(nm == e[[4]]))
}

cat("\n=== 2. Tableau des verdicts (annexe) : verdict, min kappa_L, lecture des CP, min cos2 ===\n")
exp_v <- list(
  list("S1", 5, "rv", TRUE, 1.0000, FALSE, 0.8885),
  list("S1", 5, "mfa", TRUE, 0.9999, FALSE, 0.8796),
  list("S1", 5, "statis", TRUE, 0.9487, FALSE, 0.8151),
  list("S1", 5, "nvar", TRUE, 1.0000, TRUE, 1.0000),
  list("S2", 10, "rv", TRUE, 0.9995, FALSE, 0.8055),
  list("S2", 10, "mfa", TRUE, 0.9980, FALSE, 0.4810),
  list("S2", 10, "statis", TRUE, 0.9998, FALSE, 0.9341),
  list("S2", 10, "nvar", TRUE, 0.9965, FALSE, 0.5635),
  list("S3", 5, "rv", FALSE, 0.0132, FALSE, 0.0043),
  list("S3", 5, "mfa", FALSE, 0.0076, FALSE, 0.0026),
  list("S3", 5, "statis", TRUE, 0.9995, TRUE, 0.9266),
  list("S3", 5, "nvar", TRUE, 1.0000, TRUE, 1.0000),
  list("S4", 5, "rv", FALSE, 0.7134, FALSE, 0.7069),
  list("S4", 5, "mfa", FALSE, 0.6112, FALSE, 0.5906),
  list("S4", 5, "statis", FALSE, 0.7611, FALSE, 0.7443),
  list("S4", 5, "nvar", TRUE, 1.0000, TRUE, 1.0000),
  list("vins", 3, "rv", TRUE, 0.9998, TRUE, 0.9987),
  list("vins", 3, "mfa", TRUE, 0.9992, TRUE, 0.9907),
  list("vins", 3, "statis", TRUE, 0.9996, FALSE, 0.9980),
  list("vins", 3, "nvar", TRUE, 0.9988, TRUE, 0.9951),
  list("beurres", 4, "rv", TRUE, 0.9970, FALSE, 0.5605),
  list("beurres", 4, "mfa", TRUE, 0.9848, FALSE, 0.7331),
  list("beurres", 4, "statis", TRUE, 0.9991, TRUE, 0.9936),
  list("beurres", 4, "nvar", FALSE, 0.3714, FALSE, 0.3900),
  list("beurres", 3, "rv", FALSE, 0.5562, FALSE, 0.4289),
  list("beurres", 3, "mfa", FALSE, 0.3975, FALSE, 0.2412),
  list("beurres", 3, "statis", TRUE, 0.9952, TRUE, 0.9936),
  list("beurres", 3, "nvar", FALSE, 0.2630, FALSE, 0.1179))
n_eq <- 0; n_fa <- 0
for (e in exp_v) {
  r <- compare_weighting(D[[e[[1]]]], e[[2]], e[[3]])
  n_eq <- n_eq + r$equivalent; n_fa <- n_fa + (r$equivalent && !r$pc_same)
  check(sprintf("%-8s R=%-2d %-7s %s kappa %.3f | CP %s cos2 %.3f", e[[1]], e[[2]], e[[3]],
                if (r$equivalent) "equiv." else "change", min(r$kappa), if (r$pc_same) "idem " else "change", min(r$pc_cos2)),
        identical(r$equivalent, e[[4]]) && near(min(r$kappa), e[[5]]) && identical(r$pc_same, e[[6]]) && near(min(r$pc_cos2), e[[7]]))
}
check(sprintf("28 cas : %d équivalents (attendu 19), %d fausses alertes (attendu 10)", n_eq, n_fa), n_eq == 19 && n_fa == 10)

cat("\n=== 3. Choix de R : 196 cas (tableau 4, figure choix de R) ===\n")
Rmax <- c(S1 = 8, S2 = 12, S3 = 8, S4 = 8, vins = 5, beurres = 8)
tot <- c(n = 0, eq = 0, fa = 0, below = 0, eq_below = 0, eq_above = 0, missed = 0)
for (nm in names(Rmax)) for (R in seq_len(Rmax[[nm]])) for (w in W) {
  r <- compare_weighting(D[[nm]], R, w)
  g <- r$bound$guaranteed
  tot <- tot + c(1, r$equivalent, r$equivalent && !r$pc_same, g, g && r$equivalent, !g && r$equivalent, !r$equivalent && r$pc_same)
}
print(tot)
check("196 cas, 152 équivalents, 62 fausses alertes, 0 changement manqué", all(tot[c("n", "eq", "fa", "missed")] == c(196, 152, 62, 0)))
check("sous le seuil : 92 cas dont 91 équivalents ; au-dessus : 61 équivalents sur 104", all(tot[c("below", "eq_below", "eq_above")] == c(92, 91, 61)))

cat("\n=== 4. Sensibilité au seuil alpha (figure seuils a) ===\n")
attendus <- c("S3|5|mfa", "S3|5|rv", "S4|5|mfa", "S4|5|rv", "S4|5|statis", "beurres|3|mfa", "beurres|3|nvar", "beurres|3|rv", "beurres|4|nvar")
for (a in c(0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20)) {
  ch <- character(0); fa <- 0
  for (e in cases) for (w in W) {
    r <- compare_weighting(D[[e[1]]], as.integer(e[2]), w, alpha = a)
    if (!r$equivalent) ch <- c(ch, paste(e[1], e[2], w, sep = "|"))
    fa <- fa + (r$equivalent && !r$pc_same)
  }
  check(sprintf("alpha = %4.2f : %d changements (attendu 9, mêmes cas), %d fausses alertes (attendu 9 ou 10)", a, length(ch), fa),
        setequal(ch, attendus) && fa %in% c(9, 10))
}

cat("\n=== 5. Sensibilité au seuil tau (figure seuils b) ===\n")
res <- lapply(cases, function(e) lapply(W, function(w) compare_weighting(D[[e[1]]], as.integer(e[2]), w)))
res <- unlist(res, recursive = FALSE)
nch <- sapply(seq(0.5, 0.99, by = 0.01), function(t) sum(!sapply(res, function(r) r$same_motif && min(r$kappa) >= t - 1e-12)))
names(nch) <- sprintf("%.2f", seq(0.5, 0.99, by = 0.01))
print(nch)
plateau <- names(nch)[nch == 9]
check(sprintf("plateau à 9 changements : tau de %s à %s (attendu 0.77 à 0.94)", min(plateau), max(plateau)),
      min(plateau) == "0.77" && max(plateau) == "0.94")
cat(sprintf("\n===== %d contrôles OK, %d problèmes =====\n", n_ok, n_pb))
if (n_pb > 0) cat("Merci d'envoyer les lignes PROBLEM (ou toute la sortie).\n")
