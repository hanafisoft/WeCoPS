\name{compare_weighting}
\alias{compare_weighting}
\alias{print.wecops_comparison}
\alias{summary.wecops_comparison}
\alias{print.wecops_summary}
\title{Does a weighting give the same interpretation as the default weighting?}
\description{
Applies the equivalence rule. The default analysis is read in common and proper components;
the alternative analysis is aligned on them (Procrustes rotation) and read with the same target,
labels being compared through a tolerance band around \eqn{\alpha}. The weighting is
\emph{equivalent} to the default if both analyses have the same motif and every labelled space of
the default (space spanned by the components sharing a label) has a presence at least equal to
\code{kappa_min} in the alternative subspace (see \code{\link{group_presence}}). The usual reading of the principal components one by one is also compared.
}
\usage{
compare_weighting(blocks, R, weights, alpha = 0.05, kappa_min = 0.9, band = 0.02,
                  nstart = 10, seed = 0, block_names = NULL)
\method{print}{wecops_comparison}(x, ...)
\method{summary}{wecops_comparison}(object, ...)
\method{print}{wecops_summary}(x, ...)
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{R}{number of retained components (fixed).}
\item{weights}{alternative weighting: \code{"rv"}, \code{"mfa"}, \code{"statis"}, \code{"nvar"}
or a vector of \eqn{K} positive weights relative to the default.}
\item{alpha}{salience threshold of the labels.}
\item{kappa_min}{presence threshold of the rule.}
\item{band}{half-width of the tolerance band around \code{alpha} (0: strict comparison).}
\item{nstart, seed}{starts of the concentration; see \code{\link{ccp}}.}
\item{block_names}{optional names of the blocks.}
\item{x, object}{a \code{wecops_comparison} object, or its summary.}
\item{...}{unused.}
}
\value{
An object of class \code{wecops_comparison}: a list with the verdict (\code{equivalent},
\code{same_motif}), the presence \code{kappa} of each labelled space (with \code{groups}) and of
each component (\code{kappa_components}), the motifs and labels of both analyses, the
components (\code{U_default}, \code{V_default}, \code{U_alt}, \code{V_alt}), their saliences, the
four-part decompositions and their change (\code{parts_change}), the mixing of the principal
components of both analyses, the reading of the principal components (\code{pc_same},
\code{pc_cos2}), the reweighting threshold (\code{bound}), and the canonical factorisations
(\code{factorization}: \code{default_ccp}, \code{default_pc}, \code{alt_ccp}, \code{alt_pc}).
}
\examples{
d <- load_butters()
res <- compare_weighting(d$blocks, R = 4, weights = "nvar", block_names = d$block_names)
res
}
