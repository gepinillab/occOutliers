# Helper function to create test data
create_test_data_find <- function(n = 20, type = "sf") {
  set.seed(123)
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n) + 1,
    z = rnorm(n) + 2,
    zz = rnorm(n) + 3
  )
  if (type == "sf") {
    if(requireNamespace("sf", quietly = TRUE)){
      sf_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = 4326)
      return(sf_data)
    } else {
      stop("The package sf is required to generate sf data")
    }
    
  } else {
    return(data)
  }
}

create_test_data_env_outliers <- function(n = 10, type = "sf", all_equal = FALSE) {
  set.seed(123)
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n) + 1
  )
  if (all_equal) {
    data$x[1:8] <- 1
    data$y[1:8] <- 1
    data$x[9:10] <- c(100, 101)
    data$y[9:10] <- c(100, 101)
  }
  
  if (type == "sf") {
    if(requireNamespace("sf", quietly = TRUE)){
      sf_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = 4326)
      return(sf_data)
    } else {
      stop("The package sf is required to generate sf data")
    }
    
  } else {
    return(data)
  }
}


test_that("findOutliers input validation", {
  # Invalid 'pres'
  expect_error(findOutliers(pres = "not_a_sf"),
               "pres must be an sf object with POINT geometry")
  
  # Invalid 'pres' sf geometry
  if(requireNamespace("sf", quietly = TRUE)){
    non_point_sf <- sf::st_sf(sf::st_sfc(sf::st_linestring(matrix(c(0,0,1,1),ncol = 2, byrow = TRUE))))
    expect_error(findOutliers(pres = non_point_sf),
                 "pres must be an sf object with POINT geometry")
  }
  
  # Invalid method
  sf_data <- create_test_data_find()
  expect_error(findOutliers(pres = sf_data, method = "invalid"),
               "Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
  
  # Invalid distEnvMethod
  expect_error(findOutliers(pres = sf_data, distEnvMethod = "invalid"),
               "Invalid distance. Choose from 'euclidean', 'manhattan', 'cosine', or 'mahalanobis'")
})


test_that("findOutliers basic functionality", {
  sf_data <- create_test_data_find()
  
  # Test without spatial or environmental outliers
  result_no_outliers <- findOutliers(pres = sf_data, spatial = FALSE, environmental = FALSE, verbose = FALSE)
  expect_false("out_spatial" %in% names(result_no_outliers))
  expect_false("out_env" %in% names(result_no_outliers))
  expect_s3_class(result_no_outliers, "sf")
  
  # Test with only spatial outliers
  result_spatial <- findOutliers(pres = sf_data, spatial = TRUE, environmental = FALSE, verbose = FALSE)
  expect_true("out_spatial" %in% names(result_spatial))
  expect_false("out_env" %in% names(result_spatial))
  expect_s3_class(result_spatial, "sf")
  
  # Test with only environmental outliers
  result_env <- findOutliers(pres = sf_data, spatial = FALSE, environmental = TRUE, verbose = FALSE)
  
  # Remove columns with zero variance
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data_no_geom <- sf::st_drop_geometry(sf_data)
    same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
    sf_data <- sf_data[, names(sf_data_no_geom)[same.val]]
  }
  
  if(ncol(sf_data) > 1){
    result_env <- findOutliers(pres = sf_data, spatial = FALSE, environmental = TRUE, verbose = FALSE)
    expect_false("out_spatial" %in% names(result_env))
    expect_true("out_env" %in% names(result_env))
    expect_s3_class(result_env, "sf")
  }
  
  
  # Test with both spatial and environmental outliers
  # Remove columns with zero variance
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data_no_geom <- sf::st_drop_geometry(sf_data)
    same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
    sf_data <- sf_data[, names(sf_data_no_geom)[same.val]]
  }
  
  if(ncol(sf_data) > 1){
    result_both <- findOutliers(pres = sf_data, spatial = TRUE, environmental = TRUE, verbose = FALSE)
    expect_true("out_spatial" %in% names(result_both))
    expect_true("out_env" %in% names(result_both))
    expect_s3_class(result_both, "sf")
  }
  
  # Test with different methods
  methods <- c("grubbs", "iqr", "dixon", "rosner")
  for (method in methods) {
    if (method == "rosner"){
      result_method_rosner <- findOutliers(pres = sf_data, spatial = TRUE, environmental = TRUE, method = method, verbose = FALSE, kRosner = 2)
      expect_true("out_spatial" %in% names(result_method_rosner))
      expect_true("out_env" %in% names(result_method_rosner))
      expect_s3_class(result_method_rosner, "sf")
    } else {
      # Remove columns with zero variance
      if(requireNamespace("sf", quietly = TRUE)){
        sf_data_no_geom <- sf::st_drop_geometry(sf_data)
        same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
        sf_data <- sf_data[, names(sf_data_no_geom)[same.val]]
      }
      if(ncol(sf_data) > 1){
        result_method <- findOutliers(pres = sf_data, spatial = TRUE, environmental = TRUE, method = method, verbose = FALSE)
        expect_true("out_spatial" %in% names(result_method))
        expect_true("out_env" %in% names(result_method))
        expect_s3_class(result_method, "sf")
      }
    }
  }
  
  
  # Test different distEnvMethods
  distances <- c("euclidean", "manhattan", "cosine", "mahalanobis")
  for (dist in distances){
    # Remove columns with zero variance
    if(requireNamespace("sf", quietly = TRUE)){
      sf_data_no_geom <- sf::st_drop_geometry(sf_data)
      same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
      sf_data <- sf_data[, names(sf_data_no_geom)[same.val]]
    }
    if(ncol(sf_data) > 1){
      result_dist <- findOutliers(pres = sf_data, spatial = FALSE, environmental = TRUE, distEnvMethod = dist, verbose = FALSE)
      expect_true("out_env" %in% names(result_dist))
      expect_s3_class(result_dist, "sf")
    }
  }
  
  # Test verbose = TRUE
  expect_output(findOutliers(pres = sf_data, spatial = TRUE, environmental = FALSE, verbose = TRUE),
                "\\d+ geographic outlier\\(s\\) found with method grubbs")
  expect_output(findOutliers(pres = sf_data, spatial = FALSE, environmental = TRUE, verbose = TRUE),
                "\\d+ environmental outlier\\(s\\) found with methods grubbs and euclidean \\(distance\\)\\.")
  
  # Test the warning when almost all presences are flagged as environmental outliers
  sf_data_equal <- create_test_data_env_outliers(all_equal = TRUE)
  expect_warning(findOutliers(pres = sf_data_equal, spatial = FALSE, environmental = TRUE, verbose = TRUE),
                 "Almost all presences were flagged as environmental outliers. This often happens when there are two clear outliers and all other records have the same exact environmental data.")
  
})

test_that("findOutliers output structure", {
  sf_data <- create_test_data_find()
  # Remove columns with zero variance
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data_no_geom <- sf::st_drop_geometry(sf_data)
    same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
    sf_data <- sf_data[, names(sf_data_no_geom)[same.val]]
  }
  if(ncol(sf_data) > 1){
    result <- findOutliers(pres = sf_data, spatial = TRUE, environmental = TRUE)
    expect_s3_class(result, "sf")
    expect_true("out_spatial" %in% names(result))
    expect_true("out_env" %in% names(result))
  }
})