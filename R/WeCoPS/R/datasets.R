# Datasets of the article: wines (sensory), butters and margarines (spectral), simulated data.

.read_ext <- function(name, header = TRUE) {
  f <- system.file("extdata", name, package = "WeCoPS")
  if (f == "") stop("data file not found: ", name)
  if (header) {
    d <- utils::read.csv(f, check.names = FALSE, row.names = 1)
    as.matrix(d)
  } else {
    unname(as.matrix(utils::read.csv(f, header = FALSE)))
  }
}

load_wines <- function(prep = TRUE) {
  blocks <- lapply(1:4, function(k) .read_ext(sprintf("wines_judge%d.csv", k)))
  list(blocks = if (prep) prep_blocks(blocks) else blocks,
       variables = lapply(blocks, colnames), individuals = rownames(blocks[[1]]),
       block_names = paste("judge", 1:4))
}

load_butters <- function(prep = TRUE) {
  nm <- c("VIS", "NIR", "MIR")
  raw <- lapply(nm, function(b) .read_ext(sprintf("butters_%s.csv", b)))
  blocks <- if (prep) prep_blocks(lapply(raw, snv)) else raw
  list(blocks = blocks, variables = lapply(raw, colnames), individuals = rownames(raw[[1]]),
       block_names = nm)
}

load_simulated <- function(case) {
  nb <- c(3, 6, 3, 3)
  Rs <- c(5, 10, 5, 5)
  if (!case %in% 1:4) stop("case must be 1, 2, 3 or 4")
  blocks <- lapply(seq_len(nb[case]), function(k) .read_ext(sprintf("sim%d_block%d.csv", case, k), header = FALSE))
  list(blocks = blocks, R = Rs[case], block_names = paste("block", seq_len(nb[case])))
}
