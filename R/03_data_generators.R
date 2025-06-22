# R/03_data_generators.R
source("R/05_utils.R")

#' Generate dependent data from a unimodal distribution using Clayton Copula
#' 
#' @param n Total length of the time series.
#' @param cp The true changepoint location.
#' @param params1 A list of parameters for the first segment: list(mean, sd).
#' @param params2 A list of parameters for the second segment: list(mean, sd).
#' @param alpha The Clayton copula dependency parameter.
#' @return A numeric vector representing the generated time series.
generate_unimodal_series <- function(n, cp, params1, params2, alpha) {
  
  # Helper for generating the next value in the series
  corgen_unimodal <- function(x_prev, prev_params, next_params, alpha) {
    u <- runif(1)
    # The core copula transformation
    prob <- (1 + (u^(-alpha / (alpha + 1)) - 1) * pnorm(x_prev, prev_params$mean, prev_params$sd)^(-alpha))^(-1 / alpha)
    return(qnorm(prob, mean = next_params$mean, sd = next_params$sd))
  }
  
  series <- numeric(n)
  
  # First segment
  series[1] <- rnorm(1, params1$mean, params1$sd)
  for (i in 2:(cp - 1)) {
    series[i] <- corgen_unimodal(series[i - 1], params1, params1, alpha)
  }
  
  # Changepoint
  series[cp] <- corgen_unimodal(series[cp - 1], params1, params2, alpha)
  
  # Second segment
  for (i in (cp + 1):n) {
    series[i] <- corgen_unimodal(series[i - 1], params2, params2, alpha)
  }
  
  return(series)
}

#' Generate dependent data from a mixture normal distribution using Clayton Copula
#'
#' @param n Total length of the time series.
#' @param cp The true changepoint location.
#' @param params1 A list of parameters for the first segment: list(p, mu=c(m1,m2), sigma=c(s1,s2)).
#' @param params2 A list of parameters for the second segment: list(p, mu=c(m1,m2), sigma=c(s1,s2)).
#' @param alpha The Clayton copula dependency parameter.
#' @return A numeric vector representing the generated time series.
generate_bimodal_series <- function(n, cp, params1, params2, alpha) {
  
  # Helper for generating the next value in the series
  corgen_bimodal <- function(x_prev, prev_params, next_params, alpha) {
    u <- runif(1)
    # The core copula transformation
    prob <- (1 + (u^(-alpha / (alpha + 1)) - 1) * pnormm(x_prev, prev_params$p, prev_params$mu, prev_params$sigma)^(-alpha))^(-1 / alpha)
    return(qnormm(prob, p = next_params$p, mu = next_params$mu, sigma = next_params$sigma))
  }
  
  series <- numeric(n)
  
  # First segment
  series[1] <- qnormm(runif(1), params1$p, params1$mu, params1$sigma)
  for (i in 2:(cp - 1)) {
    series[i] <- corgen_bimodal(series[i - 1], params1, params1, alpha)
  }
  
  # Changepoint
  series[cp] <- corgen_bimodal(series[cp - 1], params1, params2, alpha)
  
  # Second segment
  for (i in (cp + 1):n) {
    series[i] <- corgen_bimodal(series[i - 1], params2, params2, alpha)
  }
  
  return(series)
}
