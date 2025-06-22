# simulations/run_simulation_study.R
library(doSNOW)
library(parallel)
# 載入所有自訂函式
source("R/01_bocd_core.R")
source("R/03_data_generators.R")
source("R/04_evaluation_metrics.R")

# ---- 1. Simulation Configuration ----
N_REPS <- 50
SIM_LENGTH <- 50
TRUE_CP <- 36
COPULA_ALPHA <- 2.0
HAZARD_RATE <- 1 / SIM_LENGTH

# MCMC and Prior settings
MCMC_PARAMS <- list(n_mcmc = 6000, burn_in = 1500, proposal_sd = list(p=0.1, mu=sqrt(10), sigma=sqrt(10)))
PRIOR_PARAMS <- list(mu = list(mean=0, sd=sqrt(100)), sigma = list(shape=1, rate=0.01))

# Define the simulation case (e.g., Case 2, Scenario I from your thesis)
# Unimodal to Bimodal transformation
PARAMS1 <- list(p = 1, mu = c(0, 0), sigma = c(1, 1)) # Effectively a single normal(0,1)
PARAMS2 <- list(p = 0.85, mu = c(5, 0), sigma = c(1, 1))

# ---- 2. Parallel Setup ----
cpu_cores <- detectCores() - 1
cl <- makeCluster(cpu_cores)
registerDoSNOW(cl)

# Progress bar setup
pb <- txtProgressBar(max = N_REPS, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# ---- 3. Run Simulation ----
start_time <- Sys.time()
simulation_results <- foreach(
  i = 1:N_REPS,
  .combine = 'cbind',
  .packages = c("truncnorm"),
  .options.snow = opts
) %dopar% {
  # Load functions inside the parallel worker
  source("R/01_bocd_core.R")
  source("R/02_posterior_predictive.R")
  source("R/03_data_generators.R")
  source("R/05_utils.R")
  
  # Generate data for this replication
  # NOTE: The generate_bimodal_series function needs to be adapted for the p=1 case,
  # or you can use generate_unimodal_series here.
  # For simplicity, let's assume generate_bimodal handles it.
  set.seed(i) # for reproducibility
  sim_data <- generate_bimodal_series(
    n = SIM_LENGTH, cp = TRUE_CP,
    params1 = PARAMS1, params2 = PARAMS2,
    alpha = COPULA_ALPHA
  )
  
  # Model parameters for BOCD
  # In this case, we estimate the mixing ratio 'p'
  model_params <- list(
    fixed_params = list(mu1=5, mu2=0, sigma1=1, sigma2=1), # Fixing everything except p
    prior_params = PRIOR_PARAMS,
    mcmc_params = MCMC_PARAMS,
    copula_alpha = COPULA_ALPHA
  )
  
  # Run BOCD
  bocd_run <- run_bocd(
    data = sim_data,
    hazard_rate = HAZARD_RATE,
    model_params = model_params,
    predictive_fn = calculate_posterior_predictive_bimodal
  )
  
  # Return the most probable run lengths
  bocd_run$max_rl
}
end_time <- Sys.time()
close(pb)
stopCluster(cl)

print(paste("Simulation took:", round(end_time - start_time, 2)))

# ---- 4. Evaluate and Save Results ----
# Summarize metrics
summary_table <- summarize_simulation_results(simulation_results, TRUE_CP)
print("Simulation Summary:")
print(summary_table)

# Create results directory if it doesn't exist
if (!dir.exists("results/simulations")) {
  dir.create("results/simulations", recursive = TRUE)
}

# Save raw results and summary
write.csv(simulation_results, "results/simulations/case2_scenario_I_raw_rl.csv", row.names = FALSE)
write.csv(summary_table, "results/simulations/case2_scenario_I_summary.csv", row.names = FALSE)

# Plot median run length
png("results/simulations/case2_scenario_I_median_rl.png", width = 800, height = 600)
median_rl <- apply(simulation_results, 1, median)
plot(1:SIM_LENGTH, median_rl, type = 'l', lwd = 2,
     main = "Median Run Length (Case 2, Scenario I)",
     xlab = "Time", ylab = "Median Run Length", ylim = c(0, SIM_LENGTH))
abline(v = TRUE_CP, col = "red", lty = 2, lwd = 2)
legend("topleft", "True Changepoint", col = "red", lty = 2, bty = 'n')
dev.off()
