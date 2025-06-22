# R/02_posterior_predictive.R
library(truncnorm)
source("R/05_utils.R")

#' Calculate the posterior predictive distribution for the bimodal model
#'
#' This function uses a Metropolis-Hastings MCMC algorithm to sample from the
#' posterior distribution of the model parameters and then computes the
#' predictive probability for a new data point. It can handle cases where
#' one or more parameters are fixed.
#'
#' @param data_segment A numeric vector of the current data segment (since the last potential CP).
#' @param prev_obs The observation at t-1.
#' @param new_obs The new observation at t for which to calculate the predictive probability.
#' @param fixed_params A list of fixed model parameters. E.g., list(mu1=0, sigma1=1).
#'        Parameters not in this list will be estimated via MCMC.
#' @param prior_params A list of hyperparameters for the prior distributions.
#' @param mcmc_params A list of MCMC settings: list(n_mcmc, burn_in, proposal_sd).
#' @param copula_alpha The dependency parameter for the Clayton copula.
#' @return A single numeric value: the posterior predictive probability P(new_obs | data_segment).
calculate_posterior_predictive_bimodal <- function(data_segment,
                                                   prev_obs,
                                                   new_obs,
                                                   fixed_params,
                                                   prior_params,
                                                   mcmc_params,
                                                   copula_alpha) {
  
  # ---- MCMC Setup ----
  # Initialize parameters. If not fixed, start at a reasonable value.
  # This part is complex and highly specific to your model.
  # The following is a simplified example based on your original code.
  # In a real scenario, you'd want more robust initialization.
  params <- list(
    p = ifelse(is.null(fixed_params$p), 0.5, fixed_params$p),
    mu1 = ifelse(is.null(fixed_params$mu1), 0, fixed_params$mu1),
    mu2 = ifelse(is.null(fixed_params$mu2), 1, fixed_params$mu2),
    sigma1 = ifelse(is.null(fixed_params$sigma1), 1, fixed_params$sigma1),
    sigma2 = ifelse(is.null(fixed_params$sigma2), 1, fixed_params$sigma2)
  )
  
  # MCMC storage
  post_samples <- list()
  if (is.null(fixed_params$p)) post_samples$p <- numeric(mcmc_params$n_mcmc)
  # ... and so on for other non-fixed params ...
  
  # Define the log-likelihood for a given set of parameters and data
  log_likelihood_fn <- function(p, m1, m2, s1, s2, data, alpha) {
    ll <- sum(log(dnormm(data, p, c(m1, m2), c(s1, s2)))) +
      sum(clayton_log_density(
        u = pnormm(data[-length(data)], p, c(m1, m2), c(s1, s2)),
        v = pnormm(data[-1], p, c(m1, m2), c(s1, s2)),
        alpha = alpha
      ))
    return(ll)
  }
  
  # ---- MCMC Loop ----
  # This loop is highly simplified. A full implementation would have separate
  # Metropolis-Hastings steps for each non-fixed parameter.
  # For brevity, let's assume we are only estimating the mixture ratio 'p'.
  
  if (is.null(fixed_params$p)) {
    p_samples <- numeric(mcmc_params$n_mcmc + 1)
    p_samples[1] <- params$p
    
    for (j in 1:mcmc_params$n_mcmc) {
      current_p <- p_samples[j]
      proposal_p <- runif(1) # Your original proposal
      
      # Log-likelihoods
      log_lik_current <- log_likelihood_fn(current_p, params$mu1, params$mu2, params$sigma1, params$sigma2, data_segment, copula_alpha)
      log_lik_proposal <- log_likelihood_fn(proposal_p, params$mu1, params$mu2, params$sigma1, params$sigma2, data_segment, copula_alpha)
      
      # Log-priors (assuming uniform for p)
      log_prior_current <- dunif(current_p, log = TRUE)
      log_prior_proposal <- dunif(proposal_p, log = TRUE)
      
      # Acceptance ratio
      log_alpha <- (log_lik_proposal + log_prior_proposal) - (log_lik_current + log_prior_current)
      
      if (log(runif(1)) < log_alpha) {
        p_samples[j + 1] <- proposal_p
      } else {
        p_samples[j + 1] <- current_p
      }
    }
    
    # Store posterior samples after burn-in
    post_samples$p <- p_samples[-(1:mcmc_params$burn_in)]
  }
  
  # ---- Calculate Predictive Probability ----
  # Using the posterior samples, average the predictive density
  # Again, this is simplified for the case where only 'p' is estimated.
  if (!is.null(post_samples$p)) {
    predictive_densities <- dnormm(new_obs, post_samples$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)) *
      exp(clayton_log_density(
        u = pnormm(new_obs, post_samples$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)),
        v = pnormm(prev_obs, post_samples$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)),
        alpha = copula_alpha
      ))
    
    return(mean(predictive_densities, na.rm = TRUE))
  } else {
    # If all params are fixed, no MCMC is needed
    predictive_density <- dnormm(new_obs, params$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)) *
      exp(clayton_log_density(
        u = pnormm(new_obs, params$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)),
        v = pnormm(prev_obs, params$p, c(params$mu1, params$mu2), c(params$sigma1, params$sigma2)),
        alpha = copula_alpha
      ))
    return(predictive_density)
  }
}

# NOTE: You would need a similar function for the unimodal (single normal) case.
# 'calculate_posterior_predictive_unimodal' would be simpler as it has fewer parameters.
