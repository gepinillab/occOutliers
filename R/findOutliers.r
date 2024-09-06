#' @title Detect Spatial and Environmental Outliers in Occurrence Data
#' @description This function identifies spatial and environmental outliers in 
#' species occurrence data. It can apply different statistical methods to 
#' detect geographic outliers and environmental anomalies.
#' @details The function takes as input an `sf` object with POINT geometry and 
#' applies spatial and/or environmental outlier detection methods. Users can 
#' choose between multiple outlier detection methods, including Grubbs, Dixon, 
#' IQR, or Rosner. For environmental outlier detection, the function calculates 
#' the distance between occurrences based on either Euclidean or Gower distance.
#' @param pres An `sf` object with POINT geometry representing species 
#' occurrence data.
#' @param spatial Logical; if `TRUE`, performs spatial outlier detection.
#' @param environmental Logical; if `TRUE`, performs environmental outlier 
#' detection.
#' @param method Character; method for detecting outliers. Options are 'grubbs',
#'  'iqr', 'dixon', or 'rosner'. Only a single method can be used at a time.
#' @param distEnvMethod Character; method for calculating environmental 
#' distances. Options are 'euclidian' or 'gower'.
#' @param scaleData Logical; if `TRUE`, scales environmental data before 
#' calculating distances.
#' @param pval Numeric; user-specified p-value for assessing the significance 
#' of the test statistic.
#' @param checkPairs Logical; if `TRUE`, checks for a pair of outliers using 
#' the Grubbs test (only applicable for sample sizes < 30).
#' @param kRosner Integer; number of outliers to test for with the Rosner test. 
#' Applicable only when `method = 'rosner'`.
#' @param verbose Logical; if `TRUE`, prints messages about the number of d
#' detected outliers.
#' @author Cory Merow <cory.merow@@gmail.com>, Gonzalo E. Pinilla-Buitrago
#' @export
findOutliers <- function(pres,
                         spatial = TRUE,
                         environmental = TRUE,
                         method = 'grubbs',
                         distEnvMethod = 'euclidean',
                         scaleData = TRUE,
                         pval = 1e-5,
                         checkPairs = FALSE,
                         kRosner = NULL,
                         verbose = TRUE) {
  
  # Detect spatial outliers
  if (spatial) {
    pres <- spatialOutliers(pres = pres, 
                            pvalSet = pval, 
                            method = method, 
                            checkPairs = checkPairs, 
                            kRosner = kRosner)
    
    if (verbose) {
      sp_count <- sum(pres$out_spatial, na.rm = TRUE)
      print(paste0(sp_count, " geographic outlier(s) found with method ", method))
    }
  }
  
  # Detect environmental outliers
  if (environmental) {
    pres <- envOutliers(pres = pres, 
                        pvalSet = pval, 
                        method = method, 
                        distEnvMethod = distEnvMethod,
                        scaleData = scaleData,
                        checkPairs = checkPairs, 
                        kRosner = kRosner)
    
    if (verbose) {
      env_count <- sum(pres$out_env, na.rm = TRUE)
      print(paste0(env_count, " environmental outlier(s) found with methods ", 
                   method, " and " , distEnvMethod, " (distance)."))
      
      if (nrow(pres) - env_count < 2) {
        warning(paste0("Almost all presences were flagged as environmental ",
                       "outliers. This often happens when there are two clear ",
                       "outliers and all other records have the same exact ",
                       "environmental data."))
      }
    }
  }
  return(pres)
}
