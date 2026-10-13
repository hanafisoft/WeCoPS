# Reproduce the 28 verdicts of the article (four weightings, seven cases).
library(WeCoPS)

cases <- list(
  "Simulated 1"    = list(load_simulated(1)$blocks, 5),
  "Simulated 2"    = list(load_simulated(2)$blocks, 10),
  "Simulated 3"    = list(load_simulated(3)$blocks, 5),
  "Simulated 4"    = list(load_simulated(4)$blocks, 5),
  "Wines"          = list(load_wines()$blocks, 3),
  "Butters, R = 4" = list(load_butters()$blocks, 4),
  "Butters, R = 3" = list(load_butters()$blocks, 3)
)
weightings <- c("rv", "mfa", "statis", "nvar")

n_eq <- 0; n_fa <- 0
tab <- matrix("", length(cases), length(weightings), dimnames = list(names(cases), weightings))
for (nm in names(cases)) {
  for (w in weightings) {
    r <- compare_weighting(cases[[nm]][[1]], cases[[nm]][[2]], w)
    fa <- r$equivalent && !r$pc_same
    n_eq <- n_eq + r$equivalent
    n_fa <- n_fa + fa
    tab[nm, w] <- sprintf("%s %.2f%s", if (r$equivalent) "=" else "#", min(r$kappa), if (fa) " c" else "")
  }
}
print(noquote(tab))
cat(sprintf("\nequivalent: %d / 28; false alarms of the principal-component reading: %d / %d\n",
            n_eq, n_fa, n_eq))
cat("(= equivalent, # interpretation changed, value: min presence, c: false alarm)\n")
# Expected (article): 19 / 28 equivalent, 10 false alarms.
