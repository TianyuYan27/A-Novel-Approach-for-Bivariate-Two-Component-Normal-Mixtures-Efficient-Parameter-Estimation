# The Proposed Method - SNTO Parameter Estimation Algorithm
# Self-contained function that loads Function.R automatically

# Main parameter estimation function
The_proposed_method <- function(data, cs.bound_1 = c(0.1, 0.7, 0.7, 0.7, 0.7), 
                               cs.bound_2 = c(0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.1)) {
  
  # Load required function file
  source("Function.R")
  
  # ========== SNTO Method (1D Projection) ==========
  # Compute projection directions
  projection_directions <- compute_projection_directions(data)
  n <- nrow(data)
  sample <- matrix(NA, n, nrow(projection_directions))
  
  # Calculate projected data for each projection direction
  for (j in 1:nrow(projection_directions)) {
    sample[, j] <- data %*% projection_directions[j, ]
  }
  
  # Perform 1D SNTO parameter estimation for each projection direction
  results_list <- list()
  for (j in 1:ncol(sample)) {
    results_list[[j]] <- SNTO_Loglik(
      sample[, j], 
      plot.lb = -5, plot.ub = 10, sizeN = "large", 
      cs.kn.pop = FALSE, cs.bound = cs.bound_1, 
      ep = 1e-6, print.log = FALSE, qua = TRUE
    )$est.parms
  }
  
  # Process 1D estimation results
  results_matrix <- do.call(cbind, results_list)
  num_projections <- ncol(results_matrix) / 5
  all_cases_estimates <- list()
  ovl_values <- numeric(num_projections)
  
  # Calculate OVL value for each projection direction
  for (j in 1:num_projections) {
    idx <- (j - 1) * 5 + 1
    # Extract parameters (original order: alpha, mu1, var1, mu2, var2)
    alpha <- results_matrix[1, idx]
    mu1 <- results_matrix[1, idx + 1]
    var1 <- results_matrix[1, idx + 2]
    mu2 <- results_matrix[1, idx + 3]
    var2 <- results_matrix[1, idx + 4]
    
    # Convert variances to standard deviations for OVL calculation
    sigma1 <- sqrt(var1)
    sigma2 <- sqrt(var2)
    
    # Store parameters
    all_cases_estimates[[j]] <- list(
      case1 = c(alpha = alpha, mu1 = mu1, mu2 = mu2, sigma1 = var1, sigma2 = var2)
    )
    
    # Calculate OVL value
    ovl_values[j] <- OVL(mu1, mu2, sigma1, sigma2)
  }
  
  # Select projection direction with minimum OVL
  min_ovl_index <- which.min(ovl_values)
  min_ovl_value <- ovl_values[min_ovl_index]
  best_projection_direction <- projection_directions[min_ovl_index, ]
  best_estimates <- all_cases_estimates[[min_ovl_index]]$case1
  
  # ========== Bayesian Classification to Obtain Initial 2D Parameters ==========
  projected_data <- data %*% best_projection_direction
  
  # Extract optimal parameters
  best_alpha <- best_estimates['alpha']
  best_mu1 <- best_estimates['mu1']
  best_mu2 <- best_estimates['mu2']
  best_sigma1 <- sqrt(best_estimates['sigma1'])  # Convert variance to standard deviation
  best_sigma2 <- sqrt(best_estimates['sigma2'])
  pi2 <- 1 - best_alpha
  
  # Calculate Bayesian classification
  log_ratio <- log(best_alpha/pi2) - log(best_sigma1/best_sigma2)
  a <- 1/(2*best_sigma1^2)
  b <- 1/(2*best_sigma2^2)
  llr <- log_ratio - a*(projected_data - best_mu1)^2 + b*(projected_data - best_mu2)^2
  
  # Create subsample indices
  S1_indices <- which(llr > 0)
  S2_indices <- which(llr <= 0)
  S1 <- data[S1_indices, ]
  S2 <- data[S2_indices, ]
  
  n1 <- length(S1_indices)
  n2 <- length(S2_indices)
  
  # Calculate mean vectors and covariance matrices for each subsample
  mu1_vec <- colMeans(S1)
  mu2_vec <- colMeans(S2)
  Sigma1_mat <- cov(S1)
  Sigma2_mat <- cov(S2)
  alpha_estimated <- n1 / n
  
  # Construct THETA_initial (13-parameter form)
  THETA_initial <- c(
    mu1_vec[1], mu1_vec[2], mu2_vec[1], mu2_vec[2],        # Mean vectors
    Sigma1_mat[1,1], Sigma1_mat[1,2], Sigma1_mat[1,2], Sigma1_mat[2,2],  # Covariance matrix 1
    Sigma2_mat[1,1], Sigma2_mat[1,2], Sigma2_mat[1,2], Sigma2_mat[2,2],  # Covariance matrix 2
    alpha_estimated                                        # Mixing proportion
  )
  
  # ========== 2D SNTO Optimization ==========
  # Convert to 11-parameter decomposition form
  THETA_initial_2D <- convert_THETA_to_decomposition(THETA_initial)
  
  # Perform 2D SNTO optimization
  snto_2d_result <- SNTO_Loglik_2D(
    Initial = THETA_initial_2D, sample = data,
    plot.lb = -5, plot.ub = 10, sizeN = "normal", 
    cs.kn.pop = FALSE, cs.bound = cs.bound_2,
    ep = 1e-6, print.log = FALSE
  )
  
  # Extract final 2D parameter estimation
  THETA_2D <- snto_2d_result$THETA
  
  # ========== Return Results ==========
  return(list(
    # 1D SNTO results
    THETA_1D = THETA_initial,
    best_projection_direction = best_projection_direction,
    min_ovl_value = min_ovl_value,
    best_1d_estimates = best_estimates,
    
    # 2D SNTO results  
    THETA_2D = THETA_2D,
    THETA_initial_2D = THETA_initial_2D,
    
    # Classification results
    S1_indices = S1_indices,
    S2_indices = S2_indices,
    
    # Projection parameters
    projection_directions = projection_directions,
    ovl_values = ovl_values
  ))
}

# Usage example:
# Assume data is an n×2 data matrix
# result <- The_proposed_method(data)
# 
# Access results:
# result$THETA_2D      # The proposed method parameter estimation
# result$best_projection_direction  # Optimal projection direction