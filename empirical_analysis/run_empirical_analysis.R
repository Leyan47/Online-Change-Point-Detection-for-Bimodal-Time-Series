# empirical_analysis/run_empirical_analysis.R
library(geckor)
library(zoo)
library(mixtools) # For normalmixEM to estimate initial parameters
source("R/01_bocd_core.R")

#' Analyze cryptocurrency data for changepoints
#'
#' @param coin_id The ID of the coin from geckor (e.g., "ethereum").
#' @param train_start Start date for the training period.
#' @param train_end End date for the training period.
#' @param test_start Start date for the testing period.
#' @param test_end End date for the testing period.
#' @param analysis_name A name for this analysis, used for saving files.
#' @param model_config A list specifying which parameter to estimate ('mu1', 'sigma1', 'p', etc.).
analyze_crypto <- function(coin_id, train_start, train_end, test_start, test_end, analysis_name, model_config) {
  
  cat(paste("\n--- Starting Analysis for:", analysis_name, "---\n"))
  
  # ---- 1. Data Fetching and Preparation ----
  train_data_raw <- coin_history(coin_id, vs_currency = "usd", from = train_start, to = train_end)
  test_data_raw <- coin_history(coin_id, vs_currency = "usd", from = test_start, to = test_end)
  
  train_log_return <- diff(log(train_data_raw$price)) * 100
  test_log_return <- diff(log(test_data_raw$price)) * 100
  test_dates <- as.Date(test_data_raw$timestamp)[-1]
  
  # ---- 2. Estimate Fixed Parameters from Training Data ----
  # Using normalmixEM to get initial estimates for the bimodal model
  cat("Estimating initial parameters from training data...\n")
  mix_em_fit <- normalmixEM(train_log_return, k = 2)
  
  # Your thesis estimates Kendall's Tau for alpha, let's replicate that
  tau <- cor(train_log_return[-length(train_log_return)], train_log_return[-1], method = "kendall")
  estimated_alpha <- 2 * tau / (1 - tau) # A common transformation
  
  # Setup the fixed parameters for the BOCD run
  fixed_params <- list(
    p = mix_em_fit$lambda[1],
    mu1 = mix_em_fit$mu[1],
    mu2 = mix_em_fit$mu[2],
    sigma1 = mix_em_fit$sigma[1],
    sigma2 = mix_em_fit$sigma[2]
  )
  # Remove the parameter that we want to estimate
  fixed_params[[model_config$param_to_estimate]] <- NULL
  
  # ---- 3. Run BOCD ----
  # Define model and MCMC parameters
  bocd_model_params <- list(
    fixed_params = fixed_params,
    prior_params = list(mu = list(mean=0, sd=10), sigma = list(shape=1, rate=0.01)),
    mcmc_params = list(n_mcmc = 6500, burn_in = 1300, proposal_sd = list(p=0.1, mu=10, sigma=10)),
    copula_alpha = estimated_alpha
  )
  
  cat("Running BOCD on test data...\n")
  bocd_results <- run_bocd(
    data = test_log_return,
    hazard_rate = 1 / 150, # From your code
    model_params = bocd_model_params,
    predictive_fn = calculate_posterior_predictive_bimodal # Assuming bimodal model
  )
  
  # ---- 4. Plot and Save Results ----
  results_dir <- "results/empirical"
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  # Identify changepoints
  predicted_cps <- which(diff(bocd_results$max_rl) < 0) + 1
  cp_dates <- test_dates[predicted_cps]
  
  # Save plots
  png(file.path(results_dir, paste0(analysis_name, ".png")), width = 800, height = 1000)
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  
  # a) Log return plot
  plot(test_dates, test_log_return, type = 'l', xlab = "Date", ylab = "Log Return (%)", main = paste("Log Return of", toupper(coin_id)))
  abline(v = cp_dates, col = "red", lty = 2)
  
  # b) Run length plot
  plot(test_dates, bocd_results$max_rl, type = 'l', xlab = "Date", ylab = "Run Length", main = "BOCD Run Length Plot")
  abline(v = cp_dates, col = "red", lty = 2)
  
  # c) Price plot
  plot(as.Date(test_data_raw$timestamp), test_data_raw$price, type = 'l', xlab = "Date", ylab = "Price (USD)", main = paste("Price of", toupper(coin_id)))
  abline(v = cp_dates, col = "red", lty = 2)
  
  dev.off()
  
  cat(paste("Analysis complete. Results saved to:", results_dir, "\n"))
}

# ---- Execute Analyses ----
# Example: ETH ICO Boom, estimating the mean of the first component
analyze_crypto(
  coin_id = "ethereum",
  train_start = "2016-08-01", train_end = "2017-01-01",
  test_start = "2017-01-01", test_end = "2017-04-01",
  analysis_name = "ETH_ICO_Boom_est_mu1",
  model_config = list(param_to_estimate = "mu1")
)

# Example: BTC ICO Bust, estimating the mixture ratio 'p'
analyze_crypto(
  coin_id = "bitcoin",
  train_start = "2017-06-10", train_end = "2017-11-10",
  test_start = "2017-11-10", test_end = "2018-02-10",
  analysis_name = "BTC_ICO_Bust_est_p",
  model_config = list(param_to_estimate = "p")
)
