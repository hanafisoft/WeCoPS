\name{prep_blocks}
\alias{prep_blocks}
\alias{snv}
\title{Preprocessing of blocks}
\description{
\code{prep_blocks} centres each block by column and scales it to unit total inertia, which is
the default weighting. \code{snv} applies the standard normal variate transform to spectra
(each row centred and scaled to unit standard deviation).
}
\usage{
prep_blocks(blocks, center = TRUE, unit_inertia = TRUE)
snv(X)
}
\arguments{
\item{blocks}{a list of at least two numeric matrices (or data frames) with the same number of
rows; rows are the individuals, in the same order in every block.}
\item{center}{logical; centre each variable.}
\item{unit_inertia}{logical; divide each block by its Frobenius norm so that
\eqn{\mathrm{trace}(X_k X_k') = 1}.}
\item{X}{a numeric matrix of spectra (one spectrum per row).}
}
\value{\code{prep_blocks}: a list of matrices. \code{snv}: a matrix.}
\examples{
d <- load_wines(prep = FALSE)
blocks <- prep_blocks(d$blocks)
sapply(blocks, function(X) sum(X^2))
}
