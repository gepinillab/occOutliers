#' @title Find outlying occurrence data in geographic space
#' @description Identifies spatial outliers in species occurrence data.
#' @param pres An sf object (POINTS) describing the locations of species records. 
#'   If environmental outliers are to be detected (not implemented in this function), 
#'   the sf object should also contain the values of environmental variables to be used.
#' @param pvalSet Numeric; user-specified p-value for assessing the significance of 
#'   outlier test statistics. Default is 1e-5.
#' @param method Character; options are 'iqr', 'grubbs', 'dixon', 'rosner'. Default is 'grubbs'.
#' @param checkPairs Logical; check for a single pair of outliers using the 
#'   Grubbs test. This can only be performed for sample sizes <30. Only a single 
#'   test is used because repeating it tends to throw out more points than seem 
#'   reasonable, by eye. The value has no effect unless `method = 'grubbs'`. Default is FALSE.
#' @param kRosner Integer between 1 and 10. Determines the number of outliers 
#'   suspected with a Rosner test. The value has no effect unless `method='rosner'`. 
#'   Default is NULL.
#'
#' @return Returns the input sf object with an additional column 'out_spatial' 
#'   indicating whether each point is considered a spatial outlier (TRUE) or not (FALSE).
#'
#' @export
#'
#' @examples
#' library(sf)
#' 
#' # Read data and convert to sf object
#' myPres <- read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                                package='occOutliers'))
#' myPres <- myPres[complete.cases(myPres),]
#' myPres <- st_as_sf(myPres, coords = c("longitude", "latitude"), crs = 4326)
#' 
#' # Find spatial outliers
#' presOut <- spatialOutliers(pres = myPres, pvalSet = 1e-5)
#' 
#' @author Cory Merow <cory.merow@@gmail.com>, Gonzalo E. Pinilla-Buitrago
spatialOutliers <- function(pres,
                            pvalSet = 1e-5,
                            method = 'grubbs',
                            checkPairs = FALSE,
                            kRosner = NULL) {
  # Validate input parameters
  if (!inherits(pres, "sf") || sf::st_geometry_type(pres, by_geometry = FALSE) != "POINT") {
    stop("pres must be an sf object with POINT geometry")
  }
  
  if (!method %in% c('grubbs', 'iqr', 'dixon', 'rosner')) {
    stop("Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
  }
  
  # Create column spatial outlier if it doesn't exist
  if (!"outlier" %in% colnames(pres)) {
    pres$outlier <- FALSE
  }
  
  # Calculate distances once
  dists <- .presPercentile(pres)$dist_cent |> as.numeric()
  
  # Apply the chosen method
  pres$out_spatial <- switch(
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

# Helper functions for each method
detect_grubbs_outliers <- function(pres, dists, pvalSet, checkPairs) {
  outlier <- pres$outlier
  while (sum(!outlier) > 3) {
    gt <- outliers::grubbs.test(dists[!outlier])
    if (is.na(gt$p.value) || gt$p.value >= pvalSet) break
    outlier[!outlier][which.max(dists[!outlier])] <- TRUE
  }
  
  if (checkPairs && nrow(pres) > 3 && nrow(pres) < 31) {
    gt_pair <- outliers::grubbs.test(dists, type = 20)
    if (!is.na(gt_pair$p.value) && gt_pair$p.value < pvalSet) {
      toss <- utils::tail(order(dists), 2)
      outlier[toss] <- TRUE
    }
  }
  
  return(outlier)
}

detect_dixon_outliers <- function(pres, dists, pvalSet) {
  if (nrow(pres) < 3 || nrow(pres) > 30) {
    warning("Dixon test only applies to sample sizes between [3, 30]. Skipping this analysis.")
    return(pres$outlier)
  }
  
  if (length(unique(dists)) == 1) {
    warning("All records are the same distance from the centroid. Skipping this analysis.")
    return(pres$outlier)
  }
  
  dt <- outliers::dixon.test(dists, type = 0, two.sided = FALSE)
  outlier <- pres$outlier
  if (dt$p.value < pvalSet) {
    outlier[which.max(dists)] <- TRUE
  }
  return(outlier)
}

detect_rosner_outliers <- function(pres, dists, pvalSet, kRosner) {
  if (is.null(kRosner) || kRosner >= length(dists)) {
    warning("Invalid kRosner value. Skipping Rosner test.")
    return(pres$outlier)
  }
  
  rt <- EnvStats::rosnerTest(dists, kRosner, alpha = pvalSet)
  outlier <- pres$outlier
  if (any(rt$all.stats$Outlier)) {
    outlier[utils::tail(order(dists), kRosner)] <- TRUE
  }
  return(outlier)
}