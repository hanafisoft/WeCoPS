\name{load_wines}
\alias{load_wines}
\alias{load_butters}
\alias{load_simulated}
\title{Datasets of the article}
\description{
\code{load_wines}: free-choice profiling of 8 wines by 4 judges, with 4, 3, 4 and 3 personal
descriptors scored from 0 to 10 (Williams and Langron, 1984, as reproduced in Hanafi and Qannari,
2008); the article uses \eqn{R = 3}.

\code{load_butters}: visible (VIS, 200 variables), near infrared (NIR, 850) and mid infrared
(MIR, 1816) spectra of 21 fat products, 8 butters and 13 margarines (Jaillais, Morrin and Downey,
2007). With \code{prep = TRUE}, each spectrum is SNV-transformed, then each block is centred and
scaled to unit inertia, as in the article (\eqn{R = 4}, and \eqn{R = 3}).

\code{load_simulated}: the four simulated datasets of the article (blocks already centred and of
unit inertia): 1, three blocks with one-dimensional structures (\eqn{R = 5}); 2, six blocks of 10
to 40 variables with one-dimensional structures (\eqn{R = 10}); 3, one block with surplus proper
structures (\eqn{R = 5}); 4, two blocks with surplus proper structures (\eqn{R = 5}).
}
\usage{
load_wines(prep = TRUE)
load_butters(prep = TRUE)
load_simulated(case)
}
\arguments{
\item{prep}{logical; preprocess the blocks as in the article.}
\item{case}{1, 2, 3 or 4.}
}
\value{A list with \code{blocks}, \code{block_names} and, for the real datasets,
\code{variables} and \code{individuals}; \code{load_simulated} also returns the number \code{R}
of components used in the article.}
\references{
Williams A.A., Langron S.P. (1984). The use of free-choice profiling for the evaluation of
commercial ports. Journal of the Science of Food and Agriculture 35, 558--568.

Hanafi M., Qannari E.M. (2008). Nouvelles propriétés de l'analyse en composantes communes et
poids spécifiques. Journal de la Société Française de Statistique 149(2), 75--97.

Jaillais B., Morrin V., Downey G. (2007). Image processing of outer-product matrices -- A new way
to classify samples: examples using visible/NIR/MIR spectral data. Chemometrics and Intelligent
Laboratory Systems 86, 179--188.
}
\examples{
w <- load_wines()
sapply(w$blocks, dim)
}
