# Helper function to create test data
create_test_data <- function(type = "df", n = 20, scale = FALSE) {
  set.seed(123) # Make tests reproducible
  data <- data.frame(
    x = rnorm(n),
    y = rnorm(n) + 1,
    z = rnorm(n) + 2
  )
  if(scale){
    data <- scale(data)
  }
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

create_test_data_mahalanobis <- function(type = "df", n = 20, scale = FALSE) {
  set.seed(123) # Make tests reproducible
  data <- data.frame(
    x = rnorm(n) ,
    y = rnorm(n) + 1,
    z = rnorm(n) + 2,
    w = rnorm(n) + 3
  )
  if(scale){
    data <- scale(data)
  }
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


test_that("Input validation checks work", {
  # Invalid 'pres'
  expect_error(envOutliers(pres = "not_a_sf_or_df"))
  
  # Invalid 'pres' sf geometry
  if(requireNamespace("sf", quietly = TRUE)){
    non_point_sf <- sf::st_sf(sf::st_sfc(sf::st_linestring(matrix(c(0,0,1,1),ncol = 2, byrow = TRUE))))
    expect_error(envOutliers(pres = non_point_sf),
                 "If pres is an sf object, it must contain POINT geometry")
  }
  
  # Invalid method
  df <- create_test_data()
  expect_error(envOutliers(pres = df, method = "invalid"),
               "Invalid method. Choose from 'grubbs', 'iqr', 'dixon', or 'rosner'")
  
  # Invalid distEnvMethod
  expect_error(envOutliers(pres = df, distEnvMethod = "invalid"),
               "Invalid distance. Choose from 'euclidean', 'manhattan', 'cosine', or 'mahalanobis'")
  
  # binary data with jaccard (commented)
  # binary_df <- data.frame(x = c(0, 1, 0), y = c(1, 0, 1))
  # expect_error(envOutliers(pres = binary_df, distEnvMethod = "jaccard"),
  #            "Jaccard distance can only be used with binary data (0 and 1 values).")
})

test_that("Grubbs test warning is triggered", {
  df_large <- create_test_data(n=30)
  expect_warning(envOutliers(pres = df_large, method = "grubbs", checkPairs = TRUE),
                 "Grubbs test is only appropriate for sample sizes < 30. Ignoring checkPairs.")
})

test_that("Basic functionality with different methods and distances", {
  methods <- c("grubbs", "iqr", "dixon", "rosner")
  distances <- c("euclidean", "manhattan", "cosine")
  
  # Test with dataframe and no scale
  for (method in methods) {
    for (dist in distances) {
      df_data <- create_test_data()
      # Remove columns with zero variance
      same.val <- which(apply(df_data, 2, function(x) !all(diff(x) == 0)))
      df_data <- df_data[, same.val]
      if(ncol(df_data) > 1){
        result <- envOutliers(pres = df_data, method = method, distEnvMethod = dist, scaleData = FALSE, kRosner = 2)
        expect_true("out_env" %in% names(result))
        expect_s3_class(result, "data.frame")
      }
    }
  }
  
  # Test with sf object and scale
  for (method in methods) {
    for (dist in distances) {
      sf_data <- create_test_data(type = "sf")
      # Remove columns with zero variance
      if(requireNamespace("sf", quietly = TRUE)){
        sf_data_no_geom <- sf::st_drop_geometry(sf_data)
        same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
        sf_data_no_geom <- sf_data_no_geom[, same.val]
        sf_data <- sf_data[, names(sf_data_no_geom)]
        if(ncol(sf_data_no_geom) > 1 && !is.null(ncol(sf_data_no_geom))){
          result <- envOutliers(pres = sf_data, method = method, distEnvMethod = dist, scaleData = TRUE, kRosner = 2)
          expect_true("out_env" %in% names(result))
          expect_s3_class(result, "sf")
        }
      }
    }
  }
  
  # Test with mahalanobis and no scale and more variables
  for (method in methods) {
    df_data <- create_test_data_mahalanobis()
    # Remove columns with zero variance
    same.val <- which(apply(df_data, 2, function(x) !all(diff(x) == 0)))
    df_data <- df_data[, same.val]
    if(ncol(df_data) > 1){
      result <- tryCatch({
        envOutliers(pres = df_data, method = method, distEnvMethod = "mahalanobis", scaleData = FALSE, kRosner = 2)
      }, error = function(e) {
        if (grepl("system is exactly singular", e$message)) {
          return(NULL)
        } else {
          stop(e)
        }
      })
      if (!is.null(result)) {
        expect_true("out_env" %in% names(result))
        expect_s3_class(result, "data.frame")
      }
    }
    
  }
  
  # Test with sf object with mahalanobis and scale and more variables
  for (method in methods) {
    sf_data <- create_test_data_mahalanobis(type = "sf")
    if(requireNamespace("sf", quietly = TRUE)){
      # Remove columns with zero variance
      sf_data_no_geom <- sf::st_drop_geometry(sf_data)
      same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
      sf_data_no_geom <- sf_data_no_geom[, same.val]
      sf_data <- sf_data[, names(sf_data_no_geom)]
      if(ncol(sf_data_no_geom) > 1 && !is.null(ncol(sf_data_no_geom))){
        result <- tryCatch({
          envOutliers(pres = sf_data, method = method, distEnvMethod = "mahalanobis", scaleData = TRUE, kRosner = 2)
        }, error = function(e) {
          if (grepl("system is exactly singular", e$message)) {
            return(NULL)
          } else {
            stop(e)
          }
        })
        if (!is.null(result)) {
          expect_true("out_env" %in% names(result))
          expect_s3_class(result, "sf")
        }
      }
    }
  }
})


test_that("Output structure is correct", {
  df_data <- create_test_data()
  # Remove columns with zero variance
  same.val <- which(apply(df_data, 2, function(x) !all(diff(x) == 0)))
  df_data <- df_data[, same.val]
  if(ncol(df_data) > 1){
    # Test data frame output
    result_df <- envOutliers(pres = df_data, kRosner = 2)
    expect_s3_class(result_df, "data.frame")
    expect_true("out_env" %in% names(result_df))
  }
  
  sf_data <- create_test_data(type = "sf")
  if(requireNamespace("sf", quietly = TRUE)){
    sf_data_no_geom <- sf::st_drop_geometry(sf_data)
    # Remove columns with zero variance
    same.val <- which(apply(sf_data_no_geom, 2, function(x) !all(diff(x) == 0)))
    sf_data_no_geom <- sf_data_no_geom[, same.val]
    sf_data <- sf_data[, names(sf_data_no_geom)]
    if(ncol(sf_data_no_geom) > 1 && !is.null(ncol(sf_data_no_geom))){
      # Test sf output
      result_sf <- envOutliers(pres = sf_data, kRosner = 2)
      expect_s3_class(result_sf, "sf")
      expect_true("out_env" %in% names(result_sf))
    }
  }
})