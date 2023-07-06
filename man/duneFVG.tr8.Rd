\name{duneFVG.tr8}
\alias{duneFVG.tr8}
\docType{data}
\title{
Functional Traits of Vascular Plants Data in Northern Adriatic Dunes
}
\description{
This dataset is a list which contains the functional traits of 71 vascular plant species sampled in 128 16 m^{2} squared plots in the Northernn Adriatic Dunes. The species are the ones contained in the dataset \code{\link{duneFVG}}. See Tordoni et al. (2019) for details.
}
\usage{
data("duneFVG.tr8")
}
\format{
  A list which contains 3 objects:
  \describe{
    \item{\code{traits.tot}}{a dataframe which contains 15 functional traits in columns and 71 species as rows. Species correspond to those contained in \code{\link{duneFVG}}\code{$total}.}
    \item{\code{traits.nat}}{a dataframe which contains 15 functional traits in columns and 62 species as rows. Species correspond to those contained in \code{\link{duneFVG}}\code{$native}.}
    \item{\code{traits.ali}}{a dataframe which contains 15 functional traits in columns and 9 species as rows. Species correspond to those contained in \code{\link{duneFVG}}\code{$alien}.}
  }
  There are two linked lists to \code{duneFVG.tr8}, one with the abundance data of the species (\code{\link{duneFVG}}) and one with the spatial coordinates of the plots (\code{\link{duneFVG.xy}}).
}
\details{
The following leaf functional and mechanistic traits were measured: specific leaf area (SLA, mm2 mg-1), major and minor vein length per unit area (VLAmaj VLAmin, mm mm^{-2}), osmotic potential at full turgor (\eqn{\pi_0}, -MPa), water potential at turgor loss point (\eqn{\Psi}tlp, -MPa), N and C content (N \% and C \%), C and N stable isotope composition (\eqn{\delta}13C per thousand, and \eqn{\delta}15N per thousand) and C to N ratio (C:N).

The following categorical traits associated with phenology and reproductive strategies derived from literature are also included: life form, clonality, photosynthetic pathway, sprouting and flowering time-span.
}
\source{
Tordoni, E., Petruzzellis, F., Nardini, A., Savi, T., Bacaro, G. (2019) Make it simpler: alien species decrease functional diversity of coastal plant communities. \emph{Journal of Vegetation Science}, \bold{30}, 498--509.
}
\examples{
data(duneFVG.tr8)
}
\keyword{datasets}
