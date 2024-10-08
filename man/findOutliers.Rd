% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findOutliers.r
\name{findOutliers}
\alias{findOutliers}
\title{Detect Spatial and Environmental Outliers in Occurrence Data}
\usage{
findOutliers(
  pres,
  spatial = TRUE,
  environmental = TRUE,
  method = "grubbs",
  distEnvMethod = "euclidean",
  scaleData = TRUE,
  pval = 1e-05,
  checkPairs = FALSE,
  kRosner = NULL,
  verbose = TRUE
)
}
\arguments{
\item{pres}{An `sf` object with POINT geometry representing species 
occurrence data.}

\item{spatial}{Logical; if `TRUE`, performs spatial outlier detection.}

\item{environmental}{Logical; if `TRUE`, performs environmental outlier 
detection.}

\item{method}{Character; method for detecting outliers. Options are 'grubbs',
'iqr', 'dixon', or 'rosner'. Only a single method can be used at a time.}

\item{distEnvMethod}{Character; method for calculating environmental 
distances. Options are 'euclidian', 'manhattan', 'cosine, or 'mahalanobis'.}

\item{scaleData}{Logical; if `TRUE`, scales environmental data before 
calculating distances.}

\item{pval}{Numeric; user-specified p-value for assessing the significance 
of the test statistic.}

\item{checkPairs}{Logical; if `TRUE`, checks for a pair of outliers using 
the Grubbs test (only applicable for sample sizes < 30).}

\item{kRosner}{Integer; number of outliers to test for with the Rosner test. 
Applicable only when `method = 'rosner'`.}

\item{verbose}{Logical; if `TRUE`, prints messages about the number of d
detected outliers.}
}
\description{
This function identifies spatial and environmental outliers in 
species occurrence data. It can apply different statistical methods to 
detect geographic outliers and environmental anomalies.
}
\details{
The function takes as input an `sf` object with POINT geometry and 
applies spatial and/or environmental outlier detection methods. Users can 
choose between multiple outlier detection methods, including Grubbs, Dixon, 
IQR, or Rosner. For environmental outlier detection, the function calculates 
the distance between occurrences based on either Euclidian, Manhattan, 
Cosine, or Mahalanobis distance.
}
\author{
Cory Merow <cory.merow@gmail.com>, Gonzalo E. Pinilla-Buitrago
}
