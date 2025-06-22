# R/04_evaluation_metrics.R

#' Calculate performance metrics for changepoint detection
#'
#' @param predicted_cps A numeric vector of predicted changepoint locations.
#' @param true_cp An integer, the location of the true changepoint.
#' @param tolerance The margin of error for a correct detection.
#' @return A named list of metrics: precision, MAE, ADD, is_correct, FP, FN.
evaluate_detection <- function(predicted_cps, true_cp, tolerance = 5) {
  
  # False Negative (Miss): No changepoint was detected
  if (length(predicted_cps) == 0) {
    return(list(precision = 0, MAE = NA, ADD = NA, is_correct = 0, FP = 0, FN = 1))
  }
  
  # False Positives (Extra): Detections outside the tolerance window of the true CP
  fp <- sum(!(predicted_cps >= true_cp - tolerance & predicted_cps <= true_cp + tolerance))
  
  # True Detections: Detections within the tolerance window
  true_detections <- predicted_cps[predicted_cps >= true_cp - tolerance & predicted_cps <= true_cp + tolerance]
  
  # Precision
  precision <- ifelse(length(predicted_cps) > 0, length(true_detections) / length(predicted_cps), 0)
  
  # Mean Absolute Error (MAE) for true detections
  mae <- ifelse(length(true_detections) > 0, mean(abs(true_detections - true_cp)), NA)
  
  # Average Detection Delay (ADD)
  delays <- true_detections[true_detections >= true_cp] - true_cp
  add <- ifelse(length(delays) > 0, mean(delays), NA)
  
  # Corrected Point (is_correct): Was the true CP found?
  is_correct <- as.integer(length(true_detections) > 0)
  
  return(list(
    precision = precision,
    MAE = mae,
    ADD = add,
    is_correct = is_correct,
    FP = fp,
    FN = 0
  ))
}

#' Summarize metrics over multiple simulation runs
#'
#' @param results_matrix A matrix where each column is the run-length sequence from one simulation.
#' @param true_cp The location of the true changepoint.
#' @return A data frame summarizing the performance metrics.
summarize_simulation_results <- function(results_matrix, true_cp) {
  num_reps <- ncol(results_matrix)
  all_metrics <- list()
  
  for (i in 1:num_reps) {
    run_lengths <- results_matrix[, i]
    # Identify changepoints where run length resets to 0
    predicted_cps <- which(diff(run_lengths) < 0) + 1
    
    # Your original code used which(run_lengths == 0), which is also a valid way
    # predicted_cps <- which(run_lengths == 0)
    
    all_metrics[[i]] <- evaluate_detection(predicted_cps, true_cp)
  }
  
  # Aggregate the results
  summary_df <- data.frame(
    Precision = mean(sapply(all_metrics, `[[`, "precision"), na.rm = TRUE),
    MAE = mean(sapply(all_metrics, `[[`, "MAE"), na.rm = TRUE),
    ADD = mean(sapply(all_metrics, `[[`, "ADD"), na.rm = TRUE),
    Corrected_Points = sum(sapply(all_metrics, `[[`, "is_correct")),
    Total_FP = sum(sapply(all_metrics, `[[`, "FP")),
    Total_FN = sum(sapply(all_metrics, `[[`, "FN"))
  )
  
  return(summary_df)
}
