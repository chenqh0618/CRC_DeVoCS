
calculate_spatial_stats_3d <- function(data,
                                       feature_col,
                                       coord_cols = c("x", "y", "z"),
                                       k = 6,
                                       style = "W",
                                       zero.policy = TRUE,
                                       randomisation = TRUE) {
  require(spdep)
  # --- Input Validation ---
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!all(coord_cols %in% names(data))) {
    stop("Coordinate columns specified in 'coord_cols' not found in data.")
  }
  if (length(coord_cols) != 3) {
    stop("'coord_cols' must specify exactly three column names for x, y, and z.")
  }
  if (!feature_col %in% names(data)) {
    stop("Feature column '", feature_col, "' not found in data.")
  }
  if (!is.numeric(data[[feature_col]])) {
    stop("Feature column '", feature_col, "' must be numeric.")
  }
  if (!is.numeric(k) || k <= 0 || k >= nrow(data) || floor(k) != k) {
    stop("'k' must be a positive integer smaller than the number of samples.")
  }
  
  # Extract coordinates and feature values
  coords <- dplyr::select(data, any_of(coord_cols))
  feature_values <- data[[feature_col]]
  
  # Check for NAs
  if (any(is.na(coords))) {
    stop("Coordinate data contains NA values.")
  }
  if (any(is.na(feature_values))) {
    warning("Feature column contains NA values. Removing corresponding samples for analysis.")
    complete_cases <- complete.cases(feature_values)
    coords <- coords[complete_cases, ]
    feature_values <- feature_values[complete_cases]
    if(nrow(coords) <= k) {
      stop("Not enough non-NA samples remaining to perform analysis with k = ", k)
    }
    if(nrow(coords) < 10) { # Arbitrary small number check
      warning("Number of samples after NA removal is small (", nrow(coords), "), results may be unreliable.")
    }
  }
  
  # Check for zero variance in the feature
  if (var(feature_values, na.rm = TRUE) == 0) {
    warning("The feature variable '", feature_col, "' has zero variance. Spatial autocorrelation tests are not meaningful.")
    return(list(
      Moran = list(statistic = NA, expected = NA, p.value = NA, message = "Zero variance in feature"),
      Geary = list(statistic = NA, expected = NA, p.value = NA, message = "Zero variance in feature")
    ))
  }
  
  # --- Spatial Analysis ---
  
  # 1. Define Neighbors using k-Nearest Neighbors (kNN) on 3D coordinates
  #    `knearneigh` finds the k nearest neighbors for each point based on Euclidean distance.
  knn <- knearneigh(coords, k = k)
  
  # 2. Convert kNN object to a neighbours list (`nb` object)
  #    An `nb` object lists the neighbors for each sample index.
  nb <- knn2nb(knn)
  
  # Check for samples with no neighbours (should not happen with kNN unless k=0 or duplicates, but good practice)
  if (any(card(nb) == 0)) {
    warning("Some samples have 0 neighbours identified. This might indicate duplicate coordinates or issues with 'k'. ",
            "Using zero.policy=TRUE, but results might be affected.")
  }
  
  # 3. Create Spatial Weights List (`listw` object)
  #    This assigns weights to the neighbor relationships. `style="W"` means row-standardized
  #    (weights for each sample's neighbors sum to 1).
  #    `zero.policy=TRUE` allows processing even if some samples have no neighbors (e.g., isolated).
  listw <- nb2listw(nb, style = style, zero.policy = zero.policy)
  
  # 4. Perform Moran's I test
  #    `moran.test` calculates Moran's I, its expected value under randomness, variance, and p-value.
  #    `randomisation=TRUE` uses a permutation approach for the p-value (recommended).
  moran_result <- moran.test(feature_values, listw, randomisation = randomisation, zero.policy = zero.policy)
  
  # 5. Perform Geary's C test
  #    `geary.test` calculates Geary's C, its expected value (typically 1), variance, and p-value.
  geary_result <- geary.test(feature_values, listw, randomisation = randomisation, zero.policy = zero.policy)
  
  # 6. Format and Return Results
  
  results <- data.frame(idx = c("Moran", "Geary"),
                        statistic = c(moran_result$statistic[[1]], geary_result$statistic[[1]]),
                        expected = c(moran_result$estimate[1], geary_result$estimate[1]),
                        p.value = c(moran_result$p.value, geary_result$p.value))
  rownames(results) <- NULL
  
  return(results)
}








