#' @title Find outlying occurrence data in environmental space
#' @description Identifies environmental  outliers in species occurrence data.
#' @param pres An sf object (POINTS) or data.frame describing the environmental data
#'  of species records.
#' @param pvalSet Numeric; user-specified p-value for assessing the significance of 
#'   outlier test statistics. Default is 1e-5.
#' @param method Character; options are 'iqr', 'grubbs', 'dixon', 'rosner'. Default is 'grubbs'.
#' @param distEnvMethod Character; options are 'euclidian', 'manhattan', 'cosine,  or
#' 'mahalanobis'. Default is 'euclidean'. 
#' @param scaleData Boolean; Scale data using base::scale(). Default is TRUE.
#' @param checkPairs Logical; check for a single pair of outliers using the 
#'   Grubbs test. This can only be performed for sample sizes <30. Only a single 
#'   test is used because repeating it tends to throw out more points than seem 
#'   reasonable, by eye. The value has no effect unless `method = 'grubbs'`. Default is FALSE.
#' @param kRosner Integer between 1 and 10. Determines the number of outliers 
#'   suspected with a Rosner test. The value has no effect unless `method = 'rosner'`. 
#'   Default is NULL.
#'
#' @return Returns the input sf object or data.frame with an additional column 'out_env' 
#'   indicating whether each point is considered a environmental outlier (TRUE) or not (FALSE).
#'
#' @export
#'
#' @examples
#' library(sf)
#' library(terra)

#' # Read data and convert to sf object
#' myPres <- read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                                package='occOutliers'))
#' myPres <- myPres[complete.cases(myPres),]
#' myPres <- st_as_sf(myPres, coords = c("Longitude", "Latitude"))

#' # Load environmental data
#' env <- terra::rast(system.file('extdata/envs.tif', package = 'occOutliers'))
#' envData <- terra::extract(env, myPres, ID = FALSE)
#' myPres <- cbind(myPres, envData)
#' # Find environmental outliers
#' presOut <- envOutliers(pres = myPres, pvalSet = 1e-5)
#' 
#' @author Cory Merow <cory.merow@@gmail.com>, Gonzalo E. Pinilla-Buitrago
envOutliers <- function(pres,
                        pvalSet = 1e-5,
                        method = 'grubbs',
                        distEnvMethod = 'euclidean',
                        scaleData = TRUE,
                        checkPairs = FALSE,
                        kRosner = NULL) {
  # Validate input parameters 
  if (!(inherits(pres, "sf") || is.data.frame(pres))) {
    stop("pres must be an sf object (POINT) or a data.frame")
  }
  
  if (inherits(pres, "sf") && sf::st_geometry_type(pres, by_geometry = FALSE) != "POINT") {
    stop("If pres is an sf object, it must contain POINT geometry")
  }
  
  if (!method %in% c('grubbs', 'iqr', 'dixon', 'rosner')) {
    stop("Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
  }
  
  if (!distEnvMethod %in% c('euclidean', 'manhattan', 'cosine', 'mahalanobis')) {
    stop("Invalid distance. Choose from 'euclidean', 'manhattan', 'cosine', or 'mahalanobis'")
  }
  
  if (method == "grubbs" && nrow(pres) >= 30 && checkPairs) {
    warning("Grubbs test is only appropriate for sample sizes < 30. Ignoring checkPairs.")
    checkPairs <- FALSE
  }
  
  # if (distEnvMethod == "jaccard" && !all(apply(p.env, 2, function(x) all(x %in% c(0, 1))))) {
  #   stop("Jaccard distance can only be used with binary data (0 and 1 values).")
  # }
  
  # Create column spatial outlier if it doesn't exist
  if (!"outlier" %in% colnames(pres)) {
    pres$outlier <- FALSE
  }
  
  # Get environmental data
  if (inherits(pres, "sf")) {
    p.env <- sf::st_drop_geometry(pres)
  } else {
    p.env <- pres
  }
  
  # Scale variables if needed
  if (scaleData) {
    p.env <- base::scale(p.env)
  }
  
  # Remove variables with no variation
  same.val <- which(apply(p.env, 2, function(x) !all(is.nan(x))))
  p.env <- p.env[, same.val]
  
  # Calculate centroid
  cent <- colMeans(p.env)
  
  # Calculate the covariance matrix (for Mahalanobis)
  if (distEnvMethod == "mahalanobis") {
    cov_matrix <- cov(p.env)
  }
  
  
  # Calculate distances
  dists <- switch(
    distEnvMethod,
    euclidean = apply(p.env, 1, euclidean_dist, cent),
    manhattan = apply(p.env, 1, manhattan_dist, cent),
    cosine = apply(p.env, 1, cosine_similarity, cent),
    # jaccard = apply(p.env, 1, jaccard_dist, cent),
    # gower = apply(p.env, 1, gower_dist, cent),
    mahalanobis = apply(p.env, 1, mahalanobis_dist, cent, cov_matrix)
  )
  
  # Apply the chosen method
  pres$out_env <- switch(
    method,
    grubbs = detect_grubbs_outliers(pres, dists, pvalSet, checkPairs),
    iqr = .iqrOutlier(dists),
    dixon = detect_dixon_outliers(pres, dists, pvalSet),
    rosner = detect_rosner_outliers(pres, dists, pvalSet, kRosner)
  )
  
  # Remove outlier column
  pres <- pres |> dplyr::select(-outlier)
  
  return(pres)
}

euclidean_dist <- function(x, cent) {
  sqrt(sum((x - cent)^2))
}

manhattan_dist <- function(x, cent) {
  sum(abs(x - cent))
}

cosine_similarity <- function(x, cent) {
  sum(x * cent) / (sqrt(sum(x^2)) * sqrt(sum(cent^2)))
}

mahalanobis_dist <- function(x, cent, cov_matrix) {
  mahalanobis(x, cent, cov_matrix)
}

# jaccard_dist <- function(x, cent) {
#   intersection <- sum(x & cent)
#   union <- sum(x | cent)
#   1 - (intersection / union)
# }
# 
# gower_dist <- function(x, cent) {
#   as.matrix(cluster::daisy(as.data.frame(rbind(x, cent)), 
#                            metric = "gower"))[1, 2]
# }
