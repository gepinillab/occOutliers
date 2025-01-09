# Helper function to create test data
create_test_data_percentile <- function(n = 20, crs = 4326) {
  set.seed(123)
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n)
  )
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data <- sf::st_as_sf(data, coords = c("x", "y"), crs = crs)
    return(sf_data)
  } else {
    stop("The package sf is required to generate sf data")
  }
}

test_that(".presPercentile input validation checks work", {
  # Invalid 'xy'
  expect_error(.presPercentile(xy = "not_a_sf"),
               "The 'xy' parameter must be an 'sf' object.")
  
  # Invalid 'percent'
  sf_data <- create_test_data_percentile()
  expect_error(.presPercentile(xy = sf_data, percent = "not_numeric"),
               "The 'percent' parameter must be a single numeric value.")
  expect_error(.presPercentile(xy = sf_data, percent = c(1,2)),
               "The 'percent' parameter must be a single numeric value.")
  
  # Invalid percent value
  expect_error(.presPercentile(xy = sf_data, percent = -10),
               "The 'percent' parameter must be between 0 and 100.")
  
  # No CRS assigned
  if(requireNamespace("sf", quietly = TRUE)){
    no_crs_sf <- create_test_data_percentile(crs = NA)
    expect_error(.presPercentile(xy = no_crs_sf),
                 "The 'xy' object does not have a CRS assigned. Please assign a CRS before using this function.")
  }
  
})

test_that(".presPercentile basic functionality", {
  # Test with a percent value, check that the function returns the correct values
  sf_data <- create_test_data_percentile()
  result_with_percent <- .presPercentile(xy = sf_data, percent = 25)
  expect_true("out_quantile" %in% names(result_with_percent))
  expect_true(all(sf::st_geometry_type(result_with_percent) == "POINT"))
  
  # Test with NULL percent
  result_no_percent <- .presPercentile(xy = sf_data, percent = NULL)
  expect_false("out_quantile" %in% names(result_no_percent))
  expect_true(all(sf::st_geometry_type(result_no_percent) == "POINT"))
  
  # Test with a percent > 100
  expect_warning(.presPercentile(xy = sf_data, percent = 120),
                 "Percent value is greater than 100. Using all points.")
  
})


test_that(".presPercentile output structure", {
  sf_data <- create_test_data_percentile()
  result <- .presPercentile(xy = sf_data)
  expect_s3_class(result, "sf")
  expect_true("dist_cent" %in% names(result))
})

test_that(".iqrOutlier basic functionality", {
  # Check that it correctly identifies outliers
  dists <- c(1, 2, 3, 4, 5, 10)
  outliers <- .iqrOutlier(dists)
  expect_equal(outliers, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
  
  # Test with data with no outliers
  dists_no_outliers <- c(1,2,3,4,5)
  outliers_no_outliers <- .iqrOutlier(dists_no_outliers)
  expect_equal(outliers_no_outliers, c(FALSE, FALSE, FALSE, FALSE, FALSE))
  
})

test_that(".iqrOutlier output structure", {
  dists <- c(1, 2, 3, 4, 5, 10)
  outliers <- .iqrOutlier(dists)
  expect_type(outliers, "logical")
  expect_length(outliers, length(dists))
})
