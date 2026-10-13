\name{presence}
\alias{presence}
\alias{group_presence}
\alias{four_parts}
\alias{reweighting_bound}
\title{Presence, four-part decomposition and reweighting threshold}
\description{
\code{presence} returns \eqn{\kappa_j = \|U_{alt}'v_j\|^2}, the presence of each component
\eqn{v_j} in the subspace spanned by \code{U_alt}.

\code{group_presence} returns the presence of each labelled space \eqn{E_L} of \code{V0}
(the space spanned by the components sharing a label): the squared cosine of the largest
principal angle between \eqn{E_L} and the subspace of \code{U_alt}. It equals \eqn{\kappa_j} for
a label carried by one component, and is the quantity used by the equivalence rule.

\code{four_parts} splits the inertia of each block (unit inertia) into the parts carried by the
common components, by its own proper components, by the proper components of other blocks (and
empty components), and the residual outside the retained subspace.

\code{reweighting_bound} compares \eqn{A = \max_k w_k / \min_k w_k} with
\eqn{\lambda_R/\lambda_{R+1}}, the ratio of eigenvalues of the default analysis. If
\eqn{A < \lambda_R/\lambda_{R+1}}, every default component keeps a presence at least equal to
\eqn{(\rho_j/A - \lambda_{R+1})/(\lambda_1 - \lambda_{R+1})}, with \eqn{\rho_j = v_j'M(1)v_j}.
}
\usage{
presence(U_alt, V0)
group_presence(U_alt, V0, T0, block_names = NULL)
four_parts(S, T)
reweighting_bound(blocks, R, weights, V0 = NULL)
}
\arguments{
\item{U_alt}{orthonormal basis of the subspace of the alternative weighting.}
\item{V0}{common and proper components of the default weighting.}
\item{S}{saliences (\eqn{K \times R}).}
\item{T}{labels (\eqn{K \times R}).}
\item{T0}{labels of \code{V0}.}
\item{block_names}{optional names of the blocks.}
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{R}{number of retained components.}
\item{weights}{a weighting name or a vector of weights; see \code{\link{block_weights}}.}
}
\value{
\code{presence}: a vector. \code{group_presence}: a list of groups with \code{label},
\code{components} and \code{kappa}. \code{four_parts}: a \eqn{K \times 4} matrix.
\code{reweighting_bound}: a list with \code{A}, \code{ratio}, \code{guaranteed} and, if
\code{V0} is given, \code{kappa_lower_bound}.
}
\examples{
b <- load_simulated(1)$blocks
reweighting_bound(b, 5, "mfa")
}
