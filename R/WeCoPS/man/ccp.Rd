\name{ccp}
\alias{ccp}
\alias{concentration}
\alias{target_rotation}
\alias{block_labels}
\alias{label_names}
\alias{motif}
\alias{label_groups}
\alias{orient_groups}
\title{Common and proper components of a subspace}
\description{
\code{ccp} rotates the basis \code{U} of the retained subspace into common and proper components
\eqn{V = UQ}: (1) concentration, which maximises the INDORT criterion
\eqn{\sum_k\sum_h s_{kh}^2} restricted to the subspace (\code{concentration}, polar algorithm with
several starts); (2) target \eqn{T_{kh} = 1} if \eqn{s_{kh} \ge \alpha}; (3) rotation towards the
target (\code{target_rotation}, sweeps of plane rotations with the exact optimal angle,
minimising the inertia of each block on the components it should not carry); steps 2 and 3 are
repeated until the target is stable. This determines the space spanned by each group of
components sharing a label, but not their orientation inside it. (4) Orientation inside each
labelled space (\code{orient_groups}): the components of a group are the axes of the multiple
co-inertia analysis (Chessel and Hanafi, 1996) of the blocks of the label restricted to that space,
which makes the basis unique up to signs when the eigenvalues are distinct. The labels returned
are those of the oriented components (threshold \code{alpha} on their saliences).
\code{label_groups} returns the groups of components sharing a label.

\code{block_labels} computes binary labels from saliences, optionally with a tolerance band
around \eqn{\alpha}; \code{label_names} names each component (\code{"common 1-2"},
\code{"proper 3"}, \code{"empty"}); \code{motif} is the sorted list of labels.
}
\usage{
ccp(blocks, U, alpha = 0.05, nstart = 10, seed = 0, maxit = 20, sweeps = 100,
    weights = NULL)
concentration(Wt, nstart = 10, maxit = 1000, tol = 1e-12, seed = 0)
target_rotation(Q, Wt, T, sweeps = 100, tol = 1e-12)
orient_groups(Q, Wt, T, U, weights = NULL)
label_groups(T)
block_labels(S, alpha = 0.05, previous = NULL, band = NULL)
label_names(T, block_names = NULL)
motif(T, block_names = NULL)
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{U}{an \eqn{n \times R} orthonormal basis of the retained subspace.}
\item{alpha}{salience threshold.}
\item{nstart}{number of starts of the concentration (the first is the identity).}
\item{seed}{seed of the random starts; the random state of the session is restored.}
\item{maxit}{maximum number of iterations.}
\item{sweeps}{number of sweeps of plane rotations.}
\item{Wt}{list of projected blocks \eqn{U'X_kX_k'U}.}
\item{tol}{convergence tolerance.}
\item{Q}{an \eqn{R \times R} orthogonal matrix.}
\item{T}{a \eqn{K \times R} binary matrix of labels.}
\item{weights}{block weights of the analysis (\code{NULL}: default weighting).}
\item{S}{a \eqn{K \times R} matrix of saliences.}
\item{previous}{previous labels, used with \code{band}.}
\item{band}{\code{c(lo, hi)}: a membership changes only if the salience leaves the band.}
\item{block_names}{optional names of the blocks.}
}
\value{
\code{ccp}: a list with \code{V} (common and proper components), \code{Q}, \code{T} (labels)
and \code{S} (saliences). The other functions return a matrix or a character vector.
}
\references{
Chessel D., Hanafi M. (1996). Analyses de la co-inertie de K nuages de points. Revue de
Statistique Appliquée 44(2), 35--60.

Ten Berge J.M.F. (1984). A joint treatment of varimax rotation and the problem of diagonalizing
symmetric matrices simultaneously in the least-squares sense. Psychometrika 49, 347--358.
}
\examples{
b <- load_simulated(1)$blocks
U <- cpca(b, R = 5)$U
cc <- ccp(b, U)
label_names(cc$T)
motif(cc$T)
}
