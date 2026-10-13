\name{WeCoPS-package}
\alias{WeCoPS-package}
\alias{WeCoPS}
\docType{package}
\title{Weighting CPCA for Common and Proper Structures}
\description{
In exploratory multiblock analysis, blocks are usually scaled to unit total inertia (the default
weighting), but other weightings are common: multiple factor analysis (MFA), RV normalisation,
STATIS compromise, number of variables. WeCoPS answers the user's question: \emph{I have my
results with the default weighting; would they have changed with another weighting?}

For a fixed number \eqn{R} of components, the retained subspace of the weighted consensus PCA
(CPCA) is read in two bases: the principal components and the common and proper components
(\code{\link{ccp}}). The package measures how principal components mix common and proper
components (\code{\link{mixing}}), applies the equivalence rule (\code{\link{compare_weighting}}:
same motif of common and proper components, and presence of every default component above a
threshold), returns the canonical factorisation \eqn{X_k \approx B D_k V_k'} of each block on
both bases (\code{\link{canonical_factorization}}), and draws a diagnostic sheet
(\code{\link{plot_fiche}}).
}
\details{
Typical use: preprocess the blocks with \code{\link{prep_blocks}}, then call
\code{\link{compare_weighting}} for each alternative weighting of interest, and
\code{\link{wecops}} to read one analysis in both bases.
}
\author{Mohamed Hanafi \email{mohammed44hanafi@gmail.com}}
\references{
Hanafi M. When do two weightings give the same interpretation in common and proper structures
between several blocks of data? (submitted).
}
\keyword{package}
\keyword{multivariate}
