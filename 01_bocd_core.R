# R/01_bocd_core.R
source("R/02_posterior_predictive.R")

#' Run the Bayesian Online Changepoint Detection Algorithm
#'
#' @param data A numeric vector of the time series data.
#' @param hazard_fn A function that calculates the hazard rate. For constant hazard,
#'        use `function(r) 1/lambda`.
#' @param model_params A list containing parameters for the predictive model, including
#'        `fixed_params`, `prior_params`, `mcmc_params`, and `copula_alpha`.
#' @param predictive_fn The function to calculate posterior predictive probability
#'        (e.g., `calculate_posterior_predictive_bimodal`).
#' @return A list containing the run length probability matrix `R` and the most
#'         probable run lengths at each time step `max_rl`.
run_bocd <- function(data, hazard_rate, model_params, predictive_fn) {
  
  T <- length(data)
  # R[t, r] is the probability that the run length is r-1 at time t
  R <- matrix(0, T + 1, T + 1)
  R[1, 1] <- 1 # At t=0, run length is 0 with probability 1
  
  max_rl <- numeric(T)
  
  # Loop through each new data point
  for (t in 1:T) {
    # 1. Calculate predictive probabilities for each possible run length
    pred_probs <- numeric(t)
    for (r in 1:t) {
      # The data segment for run length r-1 at time t is data[t-r+1 ... t]
      # The previous observation is at t-1
      data_segment <- data[(t - r + 1):t]
      
      # The first posterior P(x_t | r_{t-1}=0) is a special case (integral over prior)
      if (r == 1) {
        # This part requires a separate function to integrate over the prior,
        # which was 'mcc' in your code. For simplicity, we'll approximate.
        # A proper implementation needs the mcc function.
        pred_probs[r] <- 1e-5 # Placeholder for the integral over prior
      } else {
        # Note: Your original code passed data[1:t] to the f(t) function, which
        # seems incorrect. The posterior should only depend on data since the
        # last *potential* changepoint. This is a key change.
        pred_probs[r] <- predictive_fn(
          data_segment = data_segment[-length(data_segment)], # data up to t-1
          prev_obs = data[t-1],
          new_obs = data[t],
          fixed_params = model_params$fixed_params,
          prior_params = model_params$prior_params,
          mcmc_params = model_params$mcmc_params,
          copula_alpha = model_params$copula_alpha
        )
      }
    }
    
    # 2. Calculate growth probabilities
    # P(r_t = r_{t-1} + 1 | x_{1:t})
    growth_probs <- R[t, 1:t] * pred_probs * (1 - hazard_rate)
    
    # 3. Calculate changepoint probability
    # P(r_t = 0 | x_{1:t})
    cp_prob <- sum(R[t, 1:t] * pred_probs * hazard_rate)
    
    # 4. Update the run length distribution for time t
    R[t + 1, 2:(t + 1)] <- growth_probs
    R[t + 1, 1] <- cp_prob
    
    # 5. Normalize the distribution
    R[t + 1, ] <- R[t + 1, ] / sum(R[t + 1, ], na.rm = TRUE)
    
    # 6. Store the most likely run length
    max_rl[t] <- which.max(R[t + 1, ]) - 1
    
    # Progress update
    cat(sprintf("\rProcessing t = %d / %d", t, T))
  }
  
  return(list(R = R, max_rl = max_rl))
}