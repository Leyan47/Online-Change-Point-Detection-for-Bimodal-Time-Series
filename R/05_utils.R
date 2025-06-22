# R/05_utils.R

#' Calculate the log-density of the Clayton copula
#'
#' @param u A numeric vector of probabilities (e.g., from a CDF).
#' @param v A numeric vector of probabilities (e.g., from a CDF).
#' @param alpha The Clayton copula dependency parameter.
#' @return The log-density of the copula.
clayton_log_density <- function(u, v, alpha) {
  # Constrain probabilities to avoid log(0) errors
  u <- pmin(pmax(u, 1e-4), 0.9999)
  v <- pmin(pmax(v, 1e-4), 0.9999)
  
  log_density <- log(1 + alpha) + 
    (-1 - alpha) * (log(u) + log(v)) + 
    (-1 / alpha - 2) * log(u^(-alpha) + v^(-alpha) - 1)
  
  return(log_density)
}

#' Probability density function for a two-component normal mixture
#'
#' @param x A numeric vector of quantiles.
#' @param p The mixing weight for the first component.
#' @param mu A numeric vector of the two means, c(mu1, mu2).
#' @param sigma A numeric vector of the two standard deviations, c(sigma1, sigma2).
#' @return The density of the mixture model at x.
dnormm <- function(x, p, mu, sigma) {
  return(p * dnorm(x, mu[1], sigma[1]) + (1 - p) * dnorm(x, mu[2], sigma[2]))
}

#' Cumulative distribution function for a two-component normal mixture
#'
#' @param x A numeric vector of quantiles.
#' @param p The mixing weight for the first component.
#' @param mu A numeric vector of the two means, c(mu1, mu2).
#' @param sigma A numeric vector of the two standard deviations, c(sigma1, sigma2).
#' @return The cumulative probability of the mixture model at x.
pnormm <- function(x, p, mu, sigma) {
  return(p * pnorm(x, mu[1], sigma[1]) + (1 - p) * pnorm(x, mu[2], sigma[2]))
}

#' Inverse CDF (quantile function) for a two-component normal mixture
#'
#' This function uses numerical root finding to get the quantile.
#' @param q A numeric vector of probabilities.
#' @param p The mixing weight for the first component.
#' @param mu A numeric vector of the two means, c(mu1, mu2).
#' @param sigma A numeric vector of the two standard deviations, c(sigma1, sigma2).
#' @return The quantile corresponding to probability q.
qnormm <- function(q, p, mu, sigma) {
  # Define the function whose root we want to find: pnormm(x) - q = 0
  root_fn <- function(x) {
    pnormm(x, p, mu, sigma) - q
  }
  
  # Estimate a reasonable search interval
  # A weighted average of the individual quantiles is a good starting point
  mean_est <- p * mu[1] + (1 - p) * mu[2]
  sd_est <- sqrt(p * sigma[1]^2 + (1 - p) * sigma[2]^2 + p * (1-p) * (mu[1] - mu[2])^2)
  lower_bound <- mean_est - 5 * sd_est
  upper_bound <- mean_est + 5 * sd_est
  
  # Find the root
  sapply(q, function(prob) {
    # Update the root function for the current probability
    root_fn_i <- function(x) pnormm(x, p, mu, sigma) - prob
    
    # Use uniroot to find the quantile
    uniroot(root_fn_i, interval = c(lower_bound, upper_bound))$root
  })
}
