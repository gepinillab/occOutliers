#' @title Omit Outlying Points Based on Percentile Distance from Centroid
#' @description
#' This function removes outliers from a spatial data set based on the distance
#' of each point to the centroid. Points that are farther than the specified 
#' percentile are considered outliers and are excluded.
#' @param xy An `sf` object containing the spatial points.
#' @param percent A numeric value between 0 and 100 indicating the percentile
#'   threshold for outlier removal. If NULL, no outliers are removed. Default is NULL.
#' @return A list with two elements: `xy`, the cleaned `sf` object without
#'   outliers, and `out_quantile`, an `sf` object containing the excluded outliers.
#' @export
.presPercentile <- function(xy, percent = NULL) {
  
  # Input validation
  if (!inherits(xy, "sf")) {
    stop("The 'xy' parameter must be an 'sf' object.")
  }
  
  if (!is.null(percent)) {
    if (!is.numeric(percent) || length(percent) != 1) {
      stop("The 'percent' parameter must be a single numeric value.")
    }
    
    if (percent > 100) {
      warning("Percent value is greater than 100. Using all points.")
      percent <- 100
    }
    
    if (percent < 0) {
      stop("The 'percent' parameter must be between 0 and 100.")
    }
  }
  
  # Check if CRS is assigned
  if (is.null(sf::st_crs(xy))) {
    stop(paste0("The 'xy' object does not have a CRS assigned. Please assign ",
                "a CRS before using this function."))
  }
  
  # Compute centroid
  cent <- sf::st_combine(xy) |> sf::st_centroid()
  
  # Calculate distances to centroid
  dist_cent <- sf::st_distance(xy, cent)
  xy$dist_cent <- dist_cent
  xy$out_quantile <- NULL
  # Remove outliers based on percentile
  if (!is.null(percent) && percent < 100) {
    dist_quan <- stats::quantile(dist_cent, percent / 100)
    out_indices <- (dist_cent > dist_quan)
    xy$out_quantile <- out_indices
  }
  return(xy)
}

#' find records mor than 1.5 times the interquartile range beyond the upper quartile 
#' @param dists a numeric vector
.iqrOutlier <- function(dists) {
  q3 <- stats::quantile(dists, .25)
  iqr <- stats::IQR(dists)
  return(dists > (q3 + 1.5 * iqr))
}

#' Vectorized version of grep
#' @description vectorized version of grep
#' @param pattern character string containing a regular expression (or character string for ‘fixed = TRUE’) to be matched in the given character vector.  Coerced by ‘as.character’ to a character string if possible.  If a character vector of length 2 or more is supplied, the first element is used with  a warning.
#' @param x a character vector where matches are sought, or an object which can be coerced by ‘as.character’ to a character vector.
# @export
.vgrep=function(pattern,x){mapply(function(y){grep(y,pattern)},x)}
