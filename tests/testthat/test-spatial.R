# Helper function to create test data
create_test_data_spatial <- function(n = 20, type = "sf") {
  set.seed(123) # Make tests reproducible
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
  if (type == "sf") {
    if(requireNamespace("sf", quietly = TRUE)){
      sf_data <- sf::st_as_sf(data, coords = c("x", "y"))
      return(sf_data)
    } else {
      stop("The package sf is required to generate sf data")
    }
    
  } else {
    return(data)
  }
}

# Helper function to create test data with a clear outlier for dixon test
create_test_data_dixon <- function(n = 10, outlier = TRUE) {
  set.seed(123)
  data <- data.frame(
    x = rep(1, 10),
    y = rep(2, 10)
  )
  if(outlier)
  {
    data[1, ] <- c(10,10)
  }
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data <- sf::st_as_sf(data, coords = c("x", "y"))
    return(sf_data)
  } else {
    stop("The package sf is required to generate sf data")
  }
}

# Helper function to create test data with a clear outlier for rosner test
create_test_data_rosner <- function(n = 10, k = 2, outlier = TRUE) {
  set.seed(123)
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
  if(outlier)
  {
    data[1:k, ] <- data[1:k,] + 50
  }
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data <- sf::st_as_sf(data, coords = c("x", "y"))
    return(sf_data)
  } else {
    stop("The package sf is required to generate sf data")
  }
}


test_that("Input validation checks work", {
  # Invalid 'pres'
  expect_error(spatialOutliers(pres = "not_a_sf"),
               "pres must be an sf object with POINT geometry")
  
  # Invalid 'pres' sf geometry
  if(requireNamespace("sf", quietly = TRUE)){
    non_point_sf <- sf::st_sf(sf::st_sfc(sf::st_linestring(matrix(c(0,0,1,1),ncol = 2, byrow = TRUE))))
    expect_error(spatialOutliers(pres = non_point_sf),
                 "pres must be an sf object with POINT geometry")
  }
  
  # Invalid method
  sf_data <- create_test_data_spatial()
  expect_error(spatialOutliers(pres = sf_data, method = "invalid"),
               "Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
})

test_that("Grubbs test Functionality", {
  sf_data <- create_test_data_spatial(n = 10)
  result_no_pairs <- spatialOutliers(pres = sf_data, method = "grubbs", checkPairs = FALSE)
  expect_true("out_spatial" %in% names(result_no_pairs))
  expect_s3_class(result_no_pairs, "sf")
  
  result_pairs <- spatialOutliers(pres = sf_data, method = "grubbs", checkPairs = TRUE)
  expect_true("out_spatial" %in% names(result_pairs))
  expect_s3_class(result_pairs, "sf")
  
})


test_that("Dixon test Functionality", {
  # Test that the function skips the analysis when the sample size is out of range
  sf_data_small <- create_test_data_spatial(n = 2)
  expect_warning(spatialOutliers(pres = sf_data_small, method = "dixon"))
  
  sf_data_large <- create_test_data_spatial(n = 31)
  expect_warning(spatialOutliers(pres = sf_data_large, method = "dixon"))
  
  # Test if it skips if all the values are equal
  equal_data <- create_test_data_dixon(n = 10, outlier = FALSE)
  expect_warning(spatialOutliers(pres = equal_data, method = "dixon"),
                 "All records are the same distance from the centroid. Skipping this analysis.")
  
  # Test if it identifies an outlier properly
  outlier_data <- create_test_data_dixon(n=10)
  result_dixon <- spatialOutliers(pres = outlier_data, method = "dixon")
  expect_true(any(result_dixon$out_spatial))
  expect_s3_class(result_dixon, "sf")
  
})

test_that("Rosner test Functionality", {
  # Test that the function skips when kRosner is NULL
  sf_data <- create_test_data_spatial(n = 10)
  expect_warning(spatialOutliers(pres = sf_data, method = "rosner", kRosner = NULL),
                 "Invalid kRosner value. Skipping Rosner test.")
  
  # Test that the function skips when kRosner is equal or greater than the length of the distances
  expect_warning(spatialOutliers(pres = sf_data, method = "rosner", kRosner = 10),
                 "Invalid kRosner value. Skipping Rosner test.")
  expect_warning(spatialOutliers(pres = sf_data, method = "rosner", kRosner = 11),
                 "Invalid kRosner value. Skipping Rosner test.")
  
  # Test if it identifies the outliers correctly
  outlier_data <- create_test_data_rosner(n=10, k = 2, outlier = TRUE)
  result_rosner <- suppressWarnings(spatialOutliers(pres = outlier_data, method = "rosner", kRosner = 2))
  expect_true(sum(result_rosner$out_spatial) == 2)
  expect_s3_class(result_rosner, "sf")
  
})


test_that("Basic functionality with different methods", {
  methods <- c("grubbs", "iqr", "dixon", "rosner")
  
  # Test with sf object
  for (method in methods) {
    sf_data <- create_test_data_spatial()
    if (method == "rosner"){
      result <- spatialOutliers(pres = sf_data, method = method, kRosner = 2)
      expect_true("out_spatial" %in% names(result))
      expect_s3_class(result, "sf")
    } else {
      result <- spatialOutliers(pres = sf_data, method = method)
      expect_true("out_spatial" %in% names(result))
      expect_s3_class(result, "sf")
    }
    
  }
})


test_that("Output structure is correct", {
  sf_data <- create_test_data_spatial()
  result_sf <- spatialOutliers(pres = sf_data)
  expect_s3_class(result_sf, "sf")
  expect_true("out_spatial" %in% names(result_sf))
})