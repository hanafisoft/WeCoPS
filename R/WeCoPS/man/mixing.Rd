\name{mixing}
\alias{mixing}
\title{Mixing of principal components}
\description{
Writes each principal component as a combination of common and proper components,
\eqn{U = V\Gamma} with \eqn{\Gamma = V'U}; \eqn{\gamma_{jh}^2} is the share of \eqn{u_h} carried
by \eqn{v_j}. The mixing index is \eqn{\mu_h = 1 - \max_j \gamma_{jh}^2}; \eqn{\mu^\star_h} is the
same index after pooling the components that share a label, so that it measures the mixing of
components with different labels only.
}
\usage{
mixing(U, V, T)
}
\arguments{
\item{U}{principal components (\eqn{n \times R}).}
\item{V}{common and proper components of the same subspace (\eqn{n \times R}).}
\item{T}{labels of \code{V} (\eqn{K \times R}).}
}
\value{A list with \code{gamma}, \code{mu} and \code{mu_star}.}
\examples{
a <- wecops(load_simulated(1)$blocks, R = 5)
round(a$mixing$mu_star, 2)
}
