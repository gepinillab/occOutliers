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
#'   reasonable, by eye. The value has no effect unless `method = 'grubbs'`. Default is TRUE.
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
#' presOut <- findSpatialOutliers2(pres = myPres, pvalSet = 1e-5)
#' 
#' @author Cory Merow <cory.merow@@gmail.com>, Gonzalo E. Pinilla-Buitrago
spatialOutliers <- function(pres,
                            pvalSet = 1e-5,
                            method = 'grubbs',
                            checkPairs = TRUE,
                            kRosner = NULL) {
  # Validate input parameters
  if (!inherits(pres, "sf") || sf::st_geometry_type(pres, by_geometry = FALSE) != "POINT") {
    stop("pres must be an sf object with POINT geometry")
  }
  
  if (!method %in% c('grubbs', 'iqr', 'dixon', 'rosner')) {
    stop("Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
  }
  
  # Create column spatial outlier if it doesn't exist
  if (!"out_spatial" %in% colnames(pres)) {
    pres$out_spatial <- FALSE
  }
  
  # Calculate distances once
  dists <- .presPercentile2(pres)$dist_cent |> as.numeric()
  
  # Apply the chosen method
  pres$out_spatial <- switch(
    method,
    grubbs = detect_grubbs_outliers(pres, dists, pvalSet, checkPairs),
    iqr = .iqrOutlier(dists),
    dixon = detect_dixon_outliers(pres, dists, pvalSet),
    rosner = detect_rosner_outliers(pres, dists, pvalSet, kRosner)
  )
  
  return(pres)
}

# Helper functions for each method
detect_grubbs_outliers <- function(pres, dists, pvalSet, checkPairs) {
  out_spatial <- pres$out_spatial
  while (sum(!out_spatial) > 3) {
    gt <- outliers::grubbs.test(dists[!out_spatial])
    if (is.na(gt$p.value) || gt$p.value >= pvalSet) break
    out_spatial[!out_spatial][which.max(dists[!out_spatial])] <- TRUE
  }
  
  if (checkPairs && nrow(pres) > 3 && nrow(pres) < 31) {
    gt_pair <- outliers::grubbs.test(dists, type = 20)
    if (!is.na(gt_pair$p.value) && gt_pair$p.value < pvalSet) {
      toss <- utils::tail(order(dists), 2)
      out_spatial[toss] <- TRUE
    }
  }
  
  return(out_spatial)
}

detect_dixon_outliers <- function(pres, dists, pvalSet) {
  if (nrow(pres) < 3 || nrow(pres) > 30) {
    warning("Dixon test only applies to sample sizes between [3, 30]. Skipping this analysis.")
    return(pres$out_spatial)
  }
  
  if (length(unique(dists)) == 1) {
    warning("All records are the same distance from the spatial centroid. Skipping this analysis.")
    return(pres$out_spatial)
  }
  
  dt <- outliers::dixon.test(dists, type = 0, two.sided = FALSE)
  out_spatial <- pres$out_spatial
  if (dt$p.value < pvalSet) {
    out_spatial[which.max(dists)] <- TRUE
  }
  return(out_spatial)
}

detect_rosner_outliers <- function(pres, dists, pvalSet, kRosner) {
  if (is.null(kRosner) || kRosner >= length(dists)) {
    warning("Invalid kRosner value. Skipping Rosner test.")
    return(pres$out_spatial)
  }
  
  rt <- EnvStats::rosnerTest(dists, kRosner, alpha = pvalSet)
  out_spatial <- pres$out_spatial
  if (any(rt$all.stats$Outlier)) {
    out_spatial[utils::tail(order(dists), kRosner)] <- TRUE
  }
  return(out_spatial)
}

#' @title Find outlying occurrence data in geographic space
#' @description Spatial outliers
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing 
#' the locations of species records. A `SpatialPointsDataFrame` containing the 
#' values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of 
#' Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the 
#' Grubbs test. This can only be performed for sample sizes <30. Only a single 
#' test is used because repeating it tends to throw out more points than seem 
#' reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers 
#' suspected with a Rosner test. The value has no effect unless `method='rosner'`.
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' presOut=findSpatialOutliers(pres=myPres, pvalSet=1e-5)
#' 
#' @return Returns the indices of spatial outliers
#' @author Cory Merow <cory.merow@@gmail.com>
findSpatialOutliers2 <- function(pres,
                                 pvalSet = 1e-5,
                                 method = 'grubbs',
                                 checkPairs = TRUE,
                                 kRosner = NULL) {
  # Create column spatial outlier
  if (!"out_spatial" %in% colnames(pres)) {
    pres$out_spatial <- FALSE
  }
  
  if (any(method == 'grubbs')) {
    pval <-  0
    while (pval < pvalSet & sum(!pres$out_spatial) > 3) {
      sp_pres <- pres[pres$out_spatial == FALSE, ]
      # i want to recompute the distances once an outlier is removed because the 
      # outlier biased the centroid of the group, which influences distances
      dists <- .presPercentile2(sp_pres)$dist_cent |> as.numeric()
      gt <- outliers::grubbs.test(dists)
      pval <- gt$p.value
      # conservative way to toss outliers. this checks whether the single largest
      # distance is an outlier. this is repeated until no more outliers are found
      if (is.na(pval)) {
        warning(paste0("The p-value for grubbs test was NA. Sample size is too ",
                       "small to implement the test"))
        break
      }
      if (gt$p.value < pvalSet) {
        toss <- which.max(dists)
        pres[toss, "out_spatial"] <- TRUE
      }
    }
    if (length(dists) < 3) {
      warning(paste0('All but two records were deemed outliers. The Grubbs test ',
                     'may not be appropriate for these data.'))
    }
    # toss pairs
    if (checkPairs) {
      if (nrow(pres) < 31 & nrow(pres) > 3) {
        pval <- 0
        dists <- .presPercentile2(pres)$dist_cent |> as.numeric()
        
        # By turning off this loop, I'm ensuring that you can only toss 1 pair 
        # of outliers with the loop, it tends to find lots of supposed outliers 
        # very confidently, but by eye, it tends to omit clusters
        gt <- outliers::grubbs.test(dists, type = 20)
        
        pval <- gt$p.value
        # conservative way to toss outliers. this checks whether the single 
        # largest distance is an outlier. this is repeated until no more outliers are found
        if (is.na(pval)) {
          warning(paste0(
            'The p-value for grubbs test checking for pairs of outliers was NA. ',
            'Sample size is too small to implement the test'))
          break
        }
        if (gt$p.value < pvalSet){
          toss <- utils::tail(order(dists), 2)
          # IDs in the original data frame
          pres[toss, "out_spatial"] <- TRUE
        }
      }	
    }
  } # end grubbs
  
  if (any(method == 'iqr')) {
    dists <- .presPercentile2(pres)$dist_cent |> as.numeric()
    pres$out_spatial <- .iqrOutlier(dists)
  }
  
  if (any(method == 'dixon')) {
    if (nrow(pres) < 3 | nrow(pres) > 30) {
      warning(paste0('Dixon test only applies to sample sizes between [3, 30]. ',
                     'Skipping this analysis.'))
      return(sp.toss.id)
    }
    dists <- .presPercentile2(pres)$dist_cent |> as.numeric()
    if (length(unique(dists)) == 1) {
      warning(paste0('All records are the same distance from the spatial ',
                     'centroid so outliers cannot be detected. maybe your ',
                     'records come from gridded data. Skipping this analysis.'))
      return(sp.toss.id)
    }
    dt <- outliers::dixon.test(dists, type = 0, two.sided = FALSE)
    if (dt$p.value < pvalSet) {
      pres[which.max(dists), "out_spatial"] <- TRUE
    } 
  }
  
  if (any(method == 'rosner')) {
    dists <- .presPercentile2(pres)$dist_cent |> as.numeric()
    if (kRosner >= length(dists)) {
      warning(paste0('kRosner must be an integer less than the number of presence',
                     ' records, skipping this taxon.'))
    } else {
      rt <- EnvStats::rosnerTest(dists, kRosner, alpha = pvalSet)
      if (any(rt$all.stats$Outlier)) {
        sp.toss.id <- utils::tail(order(dists), kRosner)
        pres[which.max(dists), "out_spatial"] <- TRUE
      }
    }
  }
  return(pres)
}
#======================================================================
#======================================================================
#' @title Find outlying occurrence data in geographic space
#'
#' @description Spatial outliers
#'
#' @details
#' See Examples.
#' @param pres a `SpatialPoints` or `SpatialPointsDataFrame` object describing the locations of species records. A `SpatialPointsDataFrame` containing the values of environmental variables to be used must be supplied if `envOutliers=TRUE`
#' @param method character; options are 'iqr', 'grubbs', 'dixon', 'rosner'
#' @param pvalSet user-specified p-value for assessing the significance of Grubbs test statistic.
#' @param checkPairs logical; check for a single pair of outliers using the Grubbs test. This can only be performed for sample sizes <30. Only a single test is used because repeating it tends to throw out more points than seem reasonable, by eye. The value has no effect unless `method='grubbs'`.
#' @param kRosner integer between 1 and 10. Determines the number of outliers suspected with a Rosner test. The value has no effect unless `method='rosner'`.
# @keywords
#' @export
#'
#' @examples
#' myPres=read.csv(system.file('extdata/SpeciesCSVs/Camissonia_tanacetifolia.csv',
#'                             package='occOutliers'))
#' myPres=myPres[complete.cases(myPres),]
#' sp::coordinates(myPres)=c(1,2)
#' presOut=findSpatialOutliers(pres=myPres, pvalSet=1e-5)
#' 
#' @return Returns the indices of spatial outliers
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findSpatialOutliers=function(pres,
                             pvalSet=1e-5,
                             method='grubbs',
                             checkPairs=TRUE,
                             kRosner=NULL){
  #  for testing
  #   pvalSet=1e-5; checkPairs=T
  
  pres.inliers=pres
  sp.toss.coord=NULL
  sp.toss.id=NULL
  
  if(any(method=='grubbs')){
    pval=0
    #tmp.dists=dists
    #toss singles
    while(pval<pvalSet & length(pres.inliers)>3){
    	# i want to recompute the distances once an outlier is removed because the outlier biased the centroid of the group, which influences distances
      dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
      gt=outliers::grubbs.test(dists)
      #gt=dixon.test(dists)
    	#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
      (pval=gt$p.value)
      # conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
      if(is.na(pval)) {
      	warning('p-value for grubbs test was NA. sample size is too small to implement the test')
      	break
      }
      if(gt$p.value<pvalSet){
        toss=which.max(dists)
        # IDs in the original data frame
        sp.toss.coord=rbind(sp.toss.coord,sp::coordinates(pres.inliers)[toss,])
        pres.inliers=pres.inliers[-toss,]
        #tmp.dists=tmp.dists[-toss]
      }
    }
    if(length(dists)<3) warning('All but two records were deemed outliers. The Grubbs test may not be appropriate for these data.')
    # toss pairs
    if(checkPairs){
    	if(length(pres.inliers)<31 & length(pres.inliers)>3){
  			pval=0
  			#tmp.dists=dists
      	dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid

  			# By turning off this loop, I'm ensuring that you can only toss 1 pair of outliers. with the loop, it tends to find lots of supposed outliers very confidently, but by eye, it tends to omit clusters
  			#while(pval<pvalSet){
  				gt=outliers::grubbs.test(dists,type=20)
  
  				#gt=dixon.test(dists)
  				#cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
  				(pval=gt$p.value)
  				# conservative way to toss outliers. this checks whether the single largest distance is an outlier. this is repeated until no more outliers are found
  				if(is.na(pval)) {
      			warning('p-value for grubbs test checking for pairs of outliers was NA. sample size is too small to implement the test')
      			break
      	}
  				if(gt$p.value<pvalSet){
  					toss=utils::tail(order(dists),2)
  					# IDs in the original data frame
  					sp.toss.coord=rbind(sp.toss.coord, sp::coordinates(pres.inliers)[toss,])
  					pres.inliers=pres.inliers[-toss,]
  				}
  			#}
  		}	
    }
    if(!is.null(sp.toss.coord)){
      coor=sp::coordinates(pres)
      sp.toss.id= apply(sp.toss.coord,1,function(x) which(x[1]==coor[,1] & x[2]==coor[,2]))
    } 
  } # end grubbs
  
  if(any(method=='iqr')) {
  	dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
  	sp.toss.id=.iqrOutlier(dists)
  }
  
  if(any(method=='dixon')){
  	if(length(pres.inliers)<3 | length(pres.inliers)>30 ){
  		warning('dixon test only applies to sample sizes between [3,30]. skipping this taxon')
    	return(sp.toss.id)
  	}
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(length(unique(dists))==1) {
    	warning('all records are the same distance from the spatial centroid so outliers cannot be detected. maybe your records come from gridded data. skipping this taxon')
    	return(sp.toss.id)
    }
    dt=outliers::dixon.test(dists,type=0,two.sided = FALSE)
    #cot=chisq.out.test(dists,variance = var(dists),opposite = FALSE)
    if(dt$p.value<pvalSet) sp.toss.id=which.max(dists)
  }
  
  if(any(method=='rosner')){
    dists=.presPercentile(pres.inliers,percent=NULL)[[1]]$dist.from.centroid
    if(kRosner <= length(dists)) {
    	warning('kRosner must be an integer less than the number of presence records, skipping this taxon')
    	return(sp.toss.id)
    }
    rt=EnvStats::rosnerTest(dists,kRosner,alpha=pvalSet)
    if(any(rt$all.stats$Outlier)) sp.toss.id=utils::tail(order(dists),kRosner)
  }
  
  sp.toss.id
}