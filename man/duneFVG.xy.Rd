\name{duneFVG.xy}
\alias{duneFVG.xy}
\docType{data}
\title{
Geographical Coordinates of Plots, in Northen Adriatic Dunes
}
\description{
This data set is a list that contains the geographic coordinates of 128 16 m^{2} squared plots in the Northern Adriatic Dunes where the alien and native vascular plants of the data set \code{\link{duneFVG}} were collected. See Tordoni et al. (2019) for details.}
\usage{
data("duneFVG.xy")
}
\format{
  A list that contains 3 objects:
  \describe{
    \item{\code{tot.xy}}{a numeric dataframe that contains the coordinates of 128 plots in the Northern Adriatic Dunes where all the species of \code{\link{duneFVG}}\code{$total} are collected.}
    \item{\code{nat.xy}}{a numeric dataframe that contains the coordinates of 127 plots in the Northern Adriatic Dunes where the native species of \code{\link{duneFVG}}\code{$native} are collected.}
    \item{\code{ali.xy}}{a numeric dataframe that contains the coordinates of 106 plots in the Northern Adriatic Dunes where the alien species of \code{\link{duneFVG}}\code{$alien} are collected.}
  }
  There are two linked lists to \code{duneFVG.xy}, one with the abundance data of the species (\code{\link{duneFVG}}) and one with the functional traits of the species (\code{\link{duneFVG.tr8}}).

}

\source{
Tordoni, E., Petruzzellis, F., Nardini, A., Savi, T., Bacaro, G. (2019) Make it simpler: alien species decrease functional diversity of coastal plant communities. \emph{Journal of Vegetation Science}, \bold{30}, 498--509.
}

\examples{
data(duneFVG.xy)
}
\keyword{datasets}
