\name{plot_fiche}
\alias{plot_fiche}
\alias{plot.wecops_comparison}
\title{Diagnostic sheet of a comparison}
\description{
Two-panel diagnostic sheet of a \code{\link{compare_weighting}} result: (a) presence of each
common or proper component of the default under the alternative weighting, with the threshold;
(b) change (alternative minus default) of the four parts of the inertia of each block.
}
\usage{
plot_fiche(res, file = NULL, title = NULL, width = 12, height = 5)
\method{plot}{wecops_comparison}(x, ...)
}
\arguments{
\item{res, x}{a \code{wecops_comparison} object.}
\item{file}{optional PNG file.}
\item{title}{optional title.}
\item{width, height}{size of the PNG file, in inches.}
\item{...}{passed to \code{plot_fiche}.}
}
\value{The comparison, invisibly.}
\examples{
d <- load_butters()
res <- compare_weighting(d$blocks, R = 4, weights = "nvar", block_names = d$block_names)
plot(res)
}
