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
.presPercentile2 <- function(xy, percent = NULL) {
  
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

#' omit outlying pres
#' @param xy data.frame with 2 columns
#' @param percent numeric on [0, 100]
#' @export
.presPercentile <- function (xy, percent = 95) {
  if (ncol(sp::coordinates(xy)) > 2)
    stop("xy should be defined in two dimensions")
  pfs <- sp::proj4string(xy)
  if (!is.null(percent)) {
    if (length(percent) > 1)
      stop("only one value is required for percent")
    if (percent > 100) {
      warning("Using all relocations (percent > 100)")
      percent <- 100
    }
  }
  if (inherits(xy, "SpatialPointsDataFrame")) {
    if (ncol(xy) != 1) {
      id <- factor(rep("a", nrow(as.data.frame(xy))))
    } else {
      id <- xy[[1]]
    }
  } else {
    id <- factor(rep("a", nrow(as.data.frame(xy))))
  }
  
  if (min(table(id)) < 4) stop("must have 4 records to proceed")
  id <- factor(id)
  xy <- as.data.frame(sp::coordinates(xy))
  r <- split(xy, id)
  print(r)
  est.cdg <- function(xy) apply(xy, 2, mean)
  cdg <- lapply(r, est.cdg)
  print(cdg)
  levid <- levels(id)
  print(levid)
  res <- lapply(1:length(r), function(i) {
    k <- levid[i]
    df.t <- r[[levid[i]]]
    cdg.t <- cdg[[levid[i]]]
    dist.cdg <- function(xyt) {
      d <- sqrt(((xyt[1] - cdg.t[1])^2) + ((xyt[2] - cdg.t[2])^2))
      return(d)
    }
    di <- apply(df.t, 1, dist.cdg)
    key <- c(1:length(di))
    if (!is.null(percent)) {
      acons <- key[di <= stats::quantile(di, percent / 100)]
    } else { acons = key }
    xy.t <- df.t[acons, ]
    sp::coordinates(xy.t) <- c(1,2)
    return(list(xy.t = xy.t, dist.from.centroid = di))
  })
  return(res)
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
