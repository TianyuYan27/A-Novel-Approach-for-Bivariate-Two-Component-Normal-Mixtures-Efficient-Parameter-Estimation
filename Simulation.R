# Load necessary libraries
library(MASS)
library(ggplot2)
library(viridis)
library(mvtnorm)
library(nleqslv)
library(writexl)
library(parallel)
library(doParallel)
library(foreach)

# Define progress tracking function
progress_report <- function(message) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  cat("\n[", timestamp, "] ", message, sep = "")
  cat("\n")
}

progress_report("Starting to load function files and initialize environment")
cat("Current working directory:", getwd(), "\n")  # Display current working directory

# First load function file - will stop directly if error occurs
source("Function.R")
progress_report("Successfully loaded Function.R file")

# Setup parallel processing with 32 cores
progress_report("Setting up parallel computing environment (32 cores)")
num_cores <- 32
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Prepare cluster environment
progress_report("Preparing cluster environment")
# Load libraries and source files on each cluster node
clusterEvalQ(cl, {
  library(MASS)
  library(mvtnorm)
  library(nleqslv)
  library(writexl)
  source("Function.R")
  TRUE
})
progress_report("Cluster environment preparation completed")

# Export relevant functions to cluster
clusterExport(cl, c("SNTO_Loglik_2D", "convert_THETA_to_decomposition", "MLE_2D_Mixture_Correct"))

progress_report("Defining parameter sets")
# Define parameter sets for simulation based on the table
parameter_sets <- list(
    Mirrored = list(
    mu1 = c(2, 2), 
    Sigma1 = matrix(c(2, 0.6, 0.6, 1), 2), 
    mu2 = c(-2, -2), 
    Sigma2 = matrix(c(2, 0.6, 0.6, 1), 2), 
    pi = 0.50
  ),
    Scale = list(
    mu1 = c(2, 2), 
    Sigma1 = matrix(c(2, 0.6, 0.6, 1), 2), 
    mu2 = c(2, 2), 
    Sigma2 = matrix(c(3, 0.5, 0.5, 1), 2), 
    pi = 0.60
  ),
  Location = list(
    mu1 = c(2, 2), 
    Sigma1 = matrix(c(2, 0.6, 0.6, 1), 2), 
    mu2 = c(5, 3), 
    Sigma2 = matrix(c(2, 0.6, 0.6, 1), 2), 
    pi = 0.60
  ),
  General = list(
    mu1 = c(2, 2), 
    Sigma1 = matrix(c(2, 0.6, 0.6, 1), 2), 
    mu2 = c(5, 3), 
    Sigma2 = matrix(c(3, 0.5, 0.5, 1), 2), 
    pi = 0.60
  )
)

# Define global simulation parameters
num_repeats <- 100
sample_sizes <- c(50,100,150,200)
progress_report(paste("Setting global parameters - Number of repeats:", num_repeats, 
                      "Sample sizes:", paste(sample_sizes, collapse=", ")))

# Run simulation for all parameter sets and sample sizes
results_by_set <- list()

# Run simulation for each parameter set
total_parameter_sets <- length(names(parameter_sets))
total_sample_sizes <- length(sample_sizes)
total_iterations <- total_parameter_sets * total_sample_sizes
current_iteration <- 0

for (set_name in names(parameter_sets)) {
  progress_report(paste("Starting parameter set:", set_name, "- Progress:", 
                        current_iteration + 1, "/", total_iterations, "parameter sets"))
  cat("=============================\n")
  real_params <- parameter_sets[[set_name]]
  
  # Store results for this parameter set
  results_by_set[[set_name]] <- list()
  
  # Run simulation for all sample sizes
  for (n_idx in seq_along(sample_sizes)) {
    n <- sample_sizes[n_idx]
    current_iteration <- current_iteration + 1
    progress_report(paste("Starting sample size:", n, "simulation - Total progress:", 
                          current_iteration, "/", total_iterations, 
                          "Global completion rate:", round(current_iteration/total_iterations*100), "%"))
    
    # Use foreach for parallel execution, returning all needed data
    progress_report(paste("Starting parallel computation - Parameter set:", set_name, 
                          "Sample size:", n, "Number of repeats:", num_repeats))
    
    # Modify foreach loop to return all needed data
    results <- foreach(iter = 1:num_repeats, 
                   .combine = function(x, y) {
                     list(
                       metrics = rbind(x$metrics, y$metrics)
                     )
                   },
                   .packages = c("MASS", "mvtnorm", "nleqslv", "writexl")) %dopar% {
                         # Set random seed
                         set.seed(10 + iter)
                         
                         # Parameter extraction
                         mu1_real <- real_params$mu1
                         Sigma1_real <- real_params$Sigma1
                         mu2_real <- real_params$mu2
                         Sigma2_real <- real_params$Sigma2
                         pi1_real <- real_params$pi
                         
                         # Real parameters
                         real_params_local <- list(
                           mu1 = mu1_real,
                           Sigma1 = Sigma1_real,
                           mu2 = mu2_real,
                           Sigma2 = Sigma2_real,
                           pi = pi1_real
                         )
                         
                         # Generate mixture data
                         data <- matrix(NA, n, 2)
                         cluster_assignment <- sample(1:2, n, TRUE, c(pi1_real, 1 - pi1_real))
                         for (j in 1:n) {
                           data[j, ] <- if (cluster_assignment[j] == 1) {
                             mvrnorm(1, mu1_real, Sigma1_real)
                           } else {
                             mvrnorm(1, mu2_real, Sigma2_real)
                           }
                         }
                         
                         # ========== SNTO method (1D projection) ==========
                         projection_directions <- compute_projection_directions(data)
                         sample <- matrix(NA, n, nrow(projection_directions))
                         for (j in 1:nrow(projection_directions)) {
                           sample[, j] <- data %*% projection_directions[j, ]
                         }
                         
                         # Define cs.bound_1
                         cs.bound_1 <- c(0.1, 0.7, 0.7, 0.7, 0.7)
                         
                         # SNTO parameter estimation
                         results_list <- list()
                         for (j in 1:ncol(sample)) {
                          results_list[[j]] <- SNTO_Loglik(
                            sample[, j], 
                            plot.lb = -5, plot.ub = 10, sizeN = "large", 
                            cs.kn.pop = FALSE, cs.bound = cs.bound_1, 
                            ep = 1e-6, print.log = FALSE, qua = TRUE
                          )$est.parms
                        }
                        
                        # Process results
                        results_matrix <- do.call(cbind, results_list)
                        num_projections <- ncol(results_matrix) / 5
                        # Prepare storage space for two cases
                        hat_mu_case1 <- hat_sigma_case1 <- matrix(NA, 2, num_projections)
                        hat_mu_case2 <- hat_sigma_case2 <- matrix(NA, 2, num_projections)
                        all_cases_estimates <- list()

                        # Initialize vector to store OVL values for each projection direction
                        ovl_values <- numeric(num_projections)

                        for (j in 1:num_projections) {
                          idx <- (j - 1) * 5 + 1

                          # Extract parameters for single case (original order: alpha, mu1, var1, mu2, var2)
                          alpha <- results_matrix[1, idx]
                          mu1 <- results_matrix[1, idx + 1]
                          var1 <- results_matrix[1, idx + 2]
                          mu2 <- results_matrix[1, idx + 3]
                          var2 <- results_matrix[1, idx + 4]

                          # Convert variances to standard deviations for OVL calculation
                          sigma1 <- sqrt(var1)
                          sigma2 <- sqrt(var2)

                          # Set parameters for the case
                          hat_mu_case1[, j] <- c(mu1, mu2)      # mu1, mu2
                          hat_sigma_case1[, j] <- c(var1, var2) # var1, var2

                          # Store the case parameters
                          all_cases_estimates[[j]] <- list(
                            case1 = c(alpha = alpha, mu1 = mu1, mu2 = mu2, sigma1 = var1, sigma2 = var2)
                          )

                          # Calculate OVL for this projection direction
                          ovl_values[j] <- OVL(mu1, mu2, sigma1, sigma2)
                        }

                        # Find the projection direction with minimum OVL
                        min_ovl_index <- which.min(ovl_values)
                        min_ovl_value <- ovl_values[min_ovl_index]

                        # Get the best projection direction
                        best_projection_direction <- projection_directions[min_ovl_index, ]

                        # Get the best parameter estimates
                        best_estimates <- all_cases_estimates[[min_ovl_index]]$case1

                         projected_data <- data %*% best_projection_direction
                        # Extract optimal parameters
                        best_alpha <- best_estimates['alpha']
                        best_mu1 <- best_estimates['mu1']
                        best_mu2 <- best_estimates['mu2']
                        best_sigma1 <- sqrt(best_estimates['sigma1'])  # variance to standard deviation
                        best_sigma2 <- sqrt(best_estimates['sigma2'])  # variance to standard deviation
                        pi2 <- 1 - best_alpha

                        # Calculate Bayes classification
                        log_ratio <- log(best_alpha/pi2) - log(best_sigma1/best_sigma2)
                        a <- 1/(2*best_sigma1^2)
                        b <- 1/(2*best_sigma2^2)

                        llr <- log_ratio - a*(projected_data - best_mu1)^2 + b*(projected_data - best_mu2)^2

                        # Create indices
                        S1_indices <- which(llr > 0)
                        S2_indices <- which(llr <= 0)
                        S1 <- data[S1_indices, ]
                        S2 <- data[S2_indices, ]

                        n1 <- length(S1_indices)
                        n2 <- length(S2_indices)

                        # Calculate mean vectors for each subsample
                        mu1_vec <- colMeans(S1)  # mean vector for S1
                        mu2_vec <- colMeans(S2)  # mean vector for S2

                        # Calculate covariance matrices for each subsample
                        Sigma1_mat <- cov(S1)
                        Sigma2_mat <- cov(S2)

                        # Calculate alpha (proportion of S1)
                        alpha_estimated <- n1 / n

                        # Extract individual components for THETA construction
                        mu1_1 <- mu1_vec[1]
                        mu1_2 <- mu1_vec[2]
                        mu2_1 <- mu2_vec[1]
                        mu2_2 <- mu2_vec[2]

                        sigma_1_11 <- Sigma1_mat[1,1]
                        sigma_1_12 <- Sigma1_mat[1,2]
                        sigma_1_22 <- Sigma1_mat[2,2]

                        sigma_2_11 <- Sigma2_mat[1,1]
                        sigma_2_12 <- Sigma2_mat[1,2]
                        sigma_2_22 <- Sigma2_mat[2,2]

                        # Construct THETA_initial (13 parameters)
                        THETA_initial <- c(
                          mu1_1, mu1_2, mu2_1, mu2_2,                              # mean vectors
                          sigma_1_11, sigma_1_12, sigma_1_12, sigma_1_22,          # covariance matrix for component 1
                          sigma_2_11, sigma_2_12, sigma_2_12, sigma_2_22,          # covariance matrix for component 2
                          alpha_estimated                                           # mixing proportion
                        )

                        # Convert to 11-parameter decomposition form
                        THETA_initial_2D <- convert_THETA_to_decomposition(THETA_initial)

                        # ========== 2D SNTO method ==========
                         cs.bound_2 <- c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.1)

                        snto_2d_result <- SNTO_Loglik_2D(
                          Initial = THETA_initial_2D, sample = data,
                          plot.lb = -5, plot.ub = 10, sizeN = "normal", 
                          cs.kn.pop = FALSE, cs.bound = cs.bound_2,
                          ep = 1e-6, print.log = FALSE
                        )
                        
                        # Directly extract THETA vector
                        THETA_2D <- snto_2d_result$THETA
                        
                        # Calculate true parameter vector
                        true_params <- c(mu1_real, mu2_real, 
                                        as.vector(Sigma1_real), as.vector(Sigma2_real), pi1_real)

                        # Calculate evaluation metrics for 1D SNTO method
                        kl_value <- KL(true_params, THETA_initial, 10000)
                        l2_pdf_value <- l2distancef(true_params, THETA_initial, 10000)
                        l2_cdf_value <- l2distanceF(true_params, THETA_initial, 10000)
                        
                        # Calculate evaluation metrics for 2D SNTO method
                        kl_value_2d <- KL(true_params, THETA_2D, 10000)
                        l2_pdf_value_2d <- l2distancef(true_params, THETA_2D, 10000)
                        l2_cdf_value_2d <- l2distanceF(true_params, THETA_2D, 10000)

                        # Return results
                        result_vector <- c(
                          snto_1d_kl = kl_value, snto_1d_l2_pdf = l2_pdf_value, 
                          snto_1d_l2_cdf = l2_cdf_value,
                          snto_1d_mu1_1 = THETA_initial[1], snto_1d_mu1_2 = THETA_initial[2], 
                          snto_1d_mu2_1 = THETA_initial[3], snto_1d_mu2_2 = THETA_initial[4],
                          snto_1d_sigma_1_11 = THETA_initial[5], snto_1d_sigma_1_12 = THETA_initial[6], 
                          snto_1d_sigma_1_22 = THETA_initial[8],
                          snto_1d_sigma_2_11 = THETA_initial[9], snto_1d_sigma_2_12 = THETA_initial[10], 
                          snto_1d_sigma_2_22 = THETA_initial[12],
                          snto_1d_pi = THETA_initial[13],
                          snto_2d_kl = kl_value_2d, snto_2d_l2_pdf = l2_pdf_value_2d, 
                          snto_2d_l2_cdf = l2_cdf_value_2d,
                          snto_2d_mu1_1 = THETA_2D[1], snto_2d_mu1_2 = THETA_2D[2], 
                          snto_2d_mu2_1 = THETA_2D[3], snto_2d_mu2_2 = THETA_2D[4],
                          snto_2d_sigma_1_11 = THETA_2D[5], snto_2d_sigma_1_12 = THETA_2D[6], 
                          snto_2d_sigma_1_22 = THETA_2D[8],
                          snto_2d_sigma_2_11 = THETA_2D[9], snto_2d_sigma_2_12 = THETA_2D[10], 
                          snto_2d_sigma_2_22 = THETA_2D[12],
                          snto_2d_pi = THETA_2D[13]
                        )
                        
                        list(
                          metrics = result_vector
                        )
                       }
    
    progress_report(paste("Parallel computation completed - Parameter set:", set_name, "Sample size:", n))
    
    # Process parallel results
    metrics_results <- results$metrics
    
    # Calculate mean values
    mean_results <- colMeans(metrics_results)
    
    # Print average results for this sample size
    progress_report(paste("Sample size:", n, "average results"))
    cat("-----------------------------\n")
    cat("=== 1D SNTO Method Results ===\n")
    cat("Average KL divergence:", sprintf("%.6f", mean_results["snto_1d_kl"]), "\n")
    cat("Average L2-PDF:", sprintf("%.6f", mean_results["snto_1d_l2_pdf"]), "\n")
    cat("Average L2-CDF:", sprintf("%.6f", mean_results["snto_1d_l2_cdf"]), "\n")
    cat("=== 2D SNTO Method Results ===\n")
    cat("Average KL divergence:", sprintf("%.6f", mean_results["snto_2d_kl"]), "\n")
    cat("Average L2-PDF:", sprintf("%.6f", mean_results["snto_2d_l2_pdf"]), "\n")
    cat("Average L2-CDF:", sprintf("%.6f", mean_results["snto_2d_l2_cdf"]), "\n")
    
    # Save results for this sample size
    results_by_set[[set_name]][[n_idx]] <- list(
      snto_1d_kl = mean_results["snto_1d_kl"],
      snto_1d_l2_pdf = mean_results["snto_1d_l2_pdf"],
      snto_1d_l2_cdf = mean_results["snto_1d_l2_cdf"],
      snto_2d_kl = mean_results["snto_2d_kl"],
      snto_2d_l2_pdf = mean_results["snto_2d_l2_pdf"],
      snto_2d_l2_cdf = mean_results["snto_2d_l2_cdf"]
    )
  }
  
  # Print final results for this parameter set
  progress_report(paste("Parameter set:", set_name, "final results"))
  cat("=============================\n")
  for (i in seq_along(sample_sizes)) {
    cat("\nSample size:", sample_sizes[i], "\n")
    cat("-----------------------------\n")
    cat("=== 1D SNTO Method ===\n")
    cat("Average KL divergence:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_1d_kl), "\n")
    cat("Average L2-PDF:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_1d_l2_pdf), "\n")
    cat("Average L2-CDF:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_1d_l2_cdf), "\n")
    cat("=== 2D SNTO Method ===\n")
    cat("Average KL divergence:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_2d_kl), "\n")
    cat("Average L2-PDF:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_2d_l2_pdf), "\n")
    cat("Average L2-CDF:", sprintf("%.6f", results_by_set[[set_name]][[i]]$snto_2d_l2_cdf), "\n")
  }
}

progress_report("All simulation calculations completed, preparing to export results")

# Export results to Excel
results_for_export <- list()
for (set_name in names(results_by_set)) {
  for (i in seq_along(sample_sizes)) {
    row_name <- paste(set_name, "n=", sample_sizes[i], sep="_")
    results_for_export[[row_name]] <- results_by_set[[set_name]][[i]]
  }
}

# Create 1D SNTO results dataframe
results_1d_df <- do.call(rbind, lapply(names(results_for_export), function(name) {
  snto_1d_results <- data.frame(
    Parameter_Set = name,
    SNTO_1D_KL = results_for_export[[name]]$snto_1d_kl,
    SNTO_1D_L2_PDF = results_for_export[[name]]$snto_1d_l2_pdf,
    SNTO_1D_L2_CDF = results_for_export[[name]]$snto_1d_l2_cdf
  )
  return(snto_1d_results)
}))

# Create 2D SNTO results dataframe
results_2d_df <- do.call(rbind, lapply(names(results_for_export), function(name) {
  snto_2d_results <- data.frame(
    Parameter_Set = name,
    SNTO_2D_KL = results_for_export[[name]]$snto_2d_kl,
    SNTO_2D_L2_PDF = results_for_export[[name]]$snto_2d_l2_pdf,
    SNTO_2D_L2_CDF = results_for_export[[name]]$snto_2d_l2_cdf
  )
  return(snto_2d_results)
}))

# Combine 1D and 2D results
combined_results_df <- cbind(results_1d_df, results_2d_df[, -1])

# Save results to Excel file
write_xlsx(list(
  "Combined_Results" = combined_results_df,
  "SNTO_1D_Results" = results_1d_df,
  "SNTO_2D_Results" = results_2d_df
), "simulation_results.xlsx")

progress_report("Simulation results saved to simulation_results.xlsx")

# Close parallel cluster
progress_report("Closing parallel cluster")
stopCluster(cl)

progress_report("Simulation completed - Final result file:")
cat("simulation_results.xlsx - Contains complete SNTO simulation results for all parameter sets and sample sizes\n")
cat("  - Combined_Results: Combined 1D and 2D SNTO results table\n")
cat("  - SNTO_1D_Results: 1D SNTO method detailed results\n")
cat("  - SNTO_2D_Results: 2D SNTO method detailed results\n")
