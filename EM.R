# Load required libraries
library(mixtools)
library(mvtnorm)
library(MASS)
library(readxl)
library(writexl)
library(parallel)
library(doParallel)
library(foreach)

# Progress reporting function
progress_report <- function(message) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  cat("\n[", timestamp, "] ", message, sep = "")
}

# Load functions from Function.R
progress_report("Loading Function.R")
source("Function.R")
progress_report("Function.R loaded successfully")

# Define parameter sets for simulation
progress_report("Defining parameter sets")
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
progress_report("Parameter sets defined")

# Setup parallel processing
progress_report("Setting up parallel processing with 32 cores")
num_cores <- 32
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Prepare cluster environment
progress_report("Preparing cluster environment")
clusterEvalQ(cl, {
  library(mixtools)
  library(mvtnorm)
  library(MASS)
  source("Function.R")
  TRUE
})

# Main simulation function with progress reporting and error handling
run_simulation_parallel <- function(param_name, params, n_values, num_runs = 100) {
  true_params_vector <- params_to_vector(params)
  all_results <- data.frame()
  
  for (n in n_values) {
    # Report progress for this sample size
    progress_report(paste0("Starting sample size n=", n, " for parameter set: ", param_name))
    
    # Create task list for parallel processing
    task_list <- list()
    
    # Only K-means initial values tasks
    for (i in 1:num_runs) {
      task_list[[i]] <- list(
        iter = i,
        n = n,
        method = "K-means",
        param_name = param_name,
        params = params,
        true_params_vector = true_params_vector
      )
    }
    
    # Report number of tasks created
    progress_report(paste0("Created ", length(task_list), " tasks for n=", n, ", starting parallel execution"))
    
    # Run parallel tasks with error handling and resampling
    results <- foreach(task = task_list, 
                      .combine = rbind,
                      .packages = c("mixtools", "mvtnorm", "MASS")) %dopar% {
      
      # Initialize counters
      resample_count <- 0
      error_count <- 0
      success <- FALSE
      result_row <- NULL
      max_retries <- 5  # Maximum resampling attempts
      
      for (attempt in 1:(max_retries + 1)) {
        tryCatch({
          # Set seed for reproducibility
          set.seed(10 + task$iter * 100 + attempt)  # Unique seed for each attempt
          
          # Generate sample
          mu1_real <- task$params$mu1
          Sigma1_real <- task$params$Sigma1
          mu2_real <- task$params$mu2
          Sigma2_real <- task$params$Sigma2
          pi1_real <- task$params$pi
          
          # Generate mixture data
          sample_data <- matrix(NA, task$n, 2)
          cluster_assignment <- sample(1:2, task$n, TRUE, c(pi1_real, 1 - pi1_real))
          for (j in 1:task$n) {
            sample_data[j, ] <- if (cluster_assignment[j] == 1) {
              mvrnorm(1, mu1_real, Sigma1_real)
            } else {
              mvrnorm(1, mu2_real, Sigma2_real)
            }
          }
          
          # Get K-means initial values
          initial_vals <- get_kmeans_initial_values(sample_data)
          
          # Fit EM model with K-means initial values
          em_result <- mvnormalmixEM(sample_data, 
                                     mu = initial_vals$mu,
                                     sigma = initial_vals$sigma,
                                     lambda = initial_vals$lambda,
                                     k = 2,
                                     maxit = 1000,
                                     epsilon = 1e-6,
                                     arbvar = TRUE, 
                                     verb = FALSE)
          
          estimated_params <- em_to_vector(em_result)
          
          # Calculate only L2PDF, L2CDF, and KL distances
          result_row <- data.frame(
            Parameter_Set = task$param_name,
            Sample_Size = task$n,
            Iteration = task$iter,
            Initial_Method = task$method,
            Resample_Count = resample_count,
            Error_Count = error_count,
            L2PDF = l2distancef(task$true_params_vector, estimated_params),
            L2CDF = l2distanceF(task$true_params_vector, estimated_params),
            KL = KL(task$true_params_vector, estimated_params),
            Status = "Success"
          )
          
          success <- TRUE
          break  # Break out of retry loop if successful
        }, error = function(e) {
          # Record error
          error_count <<- error_count + 1
          
          if (attempt <= max_retries) {
            # Resample and retry
            resample_count <<- resample_count + 1
          } else {
            # Last attempt failed, return error row
            result_row <<- data.frame(
              Parameter_Set = task$param_name,
              Sample_Size = task$n,
              Iteration = task$iter,
              Initial_Method = task$method,
              Resample_Count = resample_count,
              Error_Count = error_count,
              L2PDF = NA,
              L2CDF = NA,
              KL = NA,
              Status = "Failed"
            )
            success <<- FALSE
          }
        })
      }  # End of retry loop
      
      return(result_row)
    }
    
    # Report completion of this sample size
    progress_report(paste0("Completed sample size n=", n, " for ", param_name, 
                           " (", nrow(results), " simulations)"))
    all_results <- rbind(all_results, results)
  }
  
  return(all_results)
}

# Export functions to cluster (removed mmd_distance)
clusterExport(cl, c("run_simulation_parallel", "params_to_vector", "em_to_vector",
                    "get_kmeans_initial_values",
                    "l2distancef", "l2distanceF", "KL"))

# Set sample sizes and number of runs
n_values <- c(50,100,150,200)
num_runs <- 100

# Run simulations for all parameter sets
progress_report("Starting simulations for all parameter sets")
set.seed(123)  # For reproducibility
all_results <- data.frame()

# Calculate total tasks
total_params <- length(parameter_sets)
total_samples <- length(n_values)
total_iterations <- total_params * total_samples

# Run each parameter set with progress reporting
current_iteration <- 0
for (param_name in names(parameter_sets)) {
  progress_report(paste0("Starting simulations for parameter set: ", param_name,
                         " [", current_iteration + 1, "/", total_iterations, "]"))
  
  param_results <- run_simulation_parallel(
    param_name, parameter_sets[[param_name]], n_values, num_runs
  )
  
  # Report completion of parameter set
  progress_report(paste0("Completed all simulations for parameter set: ", param_name,
                         " (", nrow(param_results), " total simulations)"))
  all_results <- rbind(all_results, param_results)
  
  current_iteration <- current_iteration + 1
  progress_report(paste0("Overall progress: ", round(current_iteration/total_iterations*100, 1), "% complete"))
}

# Calculate summary statistics
progress_report("Calculating summary statistics")
summary_results <- data.frame()
param_sets <- unique(all_results$Parameter_Set)
methods <- unique(all_results$Initial_Method)

for (param_set in param_sets) {
  for (n in n_values) {
    for (method in methods) {
      subset_data <- all_results[all_results$Parameter_Set == param_set & 
                                   all_results$Sample_Size == n & 
                                   all_results$Initial_Method == method, ]
      
      # Calculate means for all metrics (excluding failed iterations)
      successful_runs <- subset_data[subset_data$Status == "Success", ]
      
      # Summarize retry and error statistics
      resample_stats <- data.frame(
        Mean_Resample_Count = mean(subset_data$Resample_Count, na.rm = TRUE),
        Mean_Error_Count = mean(subset_data$Error_Count, na.rm = TRUE),
        Failure_Rate = mean(subset_data$Status == "Failed")
      )
      
      # Calculate means for performance metrics (only keep L2PDF, L2CDF, KL)
      if (nrow(successful_runs) > 0) {
        perf_summary <- data.frame(
          Mean_L2PDF = mean(successful_runs$L2PDF, na.rm = TRUE),
          Mean_L2CDF = mean(successful_runs$L2CDF, na.rm = TRUE),
          Mean_KL = mean(successful_runs$KL, na.rm = TRUE)
        )
      } else {
        perf_summary <- data.frame(
          Mean_L2PDF = NA,
          Mean_L2CDF = NA,
          Mean_KL = NA
        )
      }
      
      # Combine all summary information
      summary_row <- cbind(
        data.frame(
          Parameter_Set = param_set,
          Sample_Size = n,
          Initial_Method = method
        ),
        resample_stats,
        perf_summary
      )
      
      summary_results <- rbind(summary_results, summary_row)
    }
  }
}

# Save results to Excel
progress_report("Saving results to Excel")
output_file <- "/root/autodl-tmp/simulation_results_kmeans_only_50_100_150.xlsx"
write_xlsx(list(
  "Detailed_Results" = all_results,
  "Summary_Results" = summary_results
), output_file)
progress_report(paste0("Results saved to: ", output_file))

# Close parallel cluster
progress_report("Closing parallel cluster")
stopCluster(cl)

# Final completion message
progress_report("SIMULATION COMPLETED SUCCESSFULLY")