### Functions for EM Algorithm
#####################################################################################
# Function to convert parameter list to vector format for distance functions
params_to_vector <- function(params) {
  c(params$mu1, params$mu2, 
    as.vector(t(params$Sigma1)), as.vector(t(params$Sigma2)), 
    params$pi)
}

# Function to convert EM results to vector format
em_to_vector <- function(em_result) {
  mu1 <- em_result$mu[[1]]
  mu2 <- em_result$mu[[2]]
  Sigma1 <- em_result$sigma[[1]]
  Sigma2 <- em_result$sigma[[2]]
  pi1 <- em_result$lambda[1]
  
  c(mu1, mu2, as.vector(t(Sigma1)), as.vector(t(Sigma2)), pi1)
}

# Function to get kmeans initial values
get_kmeans_initial_values <- function(sample_data) {
  kmeans_init <- kmeans(sample_data, centers = 2)
  
  mu_init <- list(kmeans_init$centers[1,], kmeans_init$centers[2,])
  sigma_init <- list(cov(sample_data), cov(sample_data))
  lambda_init <- as.numeric(table(kmeans_init$cluster) / nrow(sample_data))
  
  return(list(mu = mu_init, sigma = sigma_init, lambda = lambda_init))
}

### Functions for The Proposed Method
#####################################################################################
# MLE for 2 dimension
MLE_2D_Mixture <- function(data, theta) {
  # Parameter unpacking
  mu1 <- c(theta[1], theta[2])
  mu2 <- c(theta[3], theta[4])
  
  # Build covariance matrices (using diagonal elements and variance values)
  Sigma1 <- matrix(c(theta[5], theta[6], theta[6], theta[8]), nrow=2)
  Sigma2 <- matrix(c(theta[9], theta[10], theta[10], theta[12]), nrow=2)
  
  alpha <- theta[13]
  
  # Calculate multivariate normal densities
  dens1 <- dmvnorm(data, mean = mu1, sigma = Sigma1)
  dens2 <- dmvnorm(data, mean = mu2, sigma = Sigma2)
  
  # Calculate mixture density
  mixture_dens <- alpha * dens1 + (1 - alpha) * dens2
  
  # Calculate negative log-likelihood (avoid numerical underflow)
  nll <- -sum(log(pmax(mixture_dens, .Machine$double.eps)))
  
  return(nll)
}

# Define OVL for the case of normal distributions
OVL <- function(mu1, mu2, sigma1, sigma2) {
  # Calculate tau1 and tau2
  tau1 <- (mu1 * sigma2^2 - mu2 * sigma1^2 - sigma1 * sigma2 * sqrt((mu1 - mu2)^2 + (sigma2^2 - sigma1^2) * log(sigma2^2 / sigma1^2))) / (sigma2^2 - sigma1^2)
  
  tau2 <- (mu1 * sigma2^2 - mu2 * sigma1^2 + sigma1 * sigma2 * sqrt((mu1 - mu2)^2 + (sigma2^2 - sigma1^2) * log(sigma2^2 / sigma1^2))) / (sigma2^2 - sigma1^2)
  
  # Case 1: If sigma1 == sigma2
  if (sigma1 == sigma2) {
    return(2 * pnorm(-abs(mu1 - mu2) / (2 * sigma1)))
  }
  
  # Case 2: If sigma1 != sigma2
  return(1 + pnorm((tau1 - mu1) / sigma1) + pnorm((tau2 - mu2) / sigma2) - pnorm((tau1 - mu2) / sigma2) - pnorm((tau2 - mu1) / sigma1))
}

# Define the function to calculate 18 projection directions uniformly distributed on unit semicircle
compute_projection_directions <- function(data) {
  # Number of projection directions
  num_directions <- 18
  
  # Generate 18 angles uniformly distributed from 0 to π (semicircle)
  # This covers angles from 0° to 180°
  angles <- seq(0, pi, length.out = num_directions)
  
  # Convert angles to unit vectors on the semicircle
  # Each row represents a direction vector (cos(θ), sin(θ))
  # cos(θ) gives the x-component, sin(θ) gives the y-component
  projection_directions <- cbind(cos(angles), sin(angles))
  
  return(projection_directions)
}

### Functions for HD Estimator
# Quantile estimator function
quantile_estimator<-function(sample){
  sample<-sort(sample)
  k<-length(sample)
  p<-(2*seq(1:k)-1)/(2*k)
  i=1:k
  quantile<-c()
  for (j in 1:k){
    quantile<-c(quantile,
                sum((pbeta(i/k,p[j]*(k+1),(1-p[j])*(k+1))-
                       pbeta((i-1)/k,p[j]*(k+1),(1-p[j])*(k+1)))*sample))
  }
  return(quantile)
}

# Generate initial parameter values
Initial_thetas <- function(sample){
  xbar<-mean(sample)
  THETA0 <- c(length(sample[sample<xbar])/length(sample),
              mean(sample[sample<xbar]),
              sd(sample[sample<xbar])^2,
              mean(sample[sample>xbar]),
              sd(sample[sample>xbar])^2)
  return(THETA0)
}

# Generate GLP from Up
NT_Net_glp <- function(n,h,p){
  if (p>=n){
    print("ERROR: n should be larger than p!")
  }else{
    # Step 2: q matrix
    k <- seq(1:n)
    h <- as.matrix(h)
    k <- as.matrix(k)
    q_matrix <- k%*%t(h)%%n
    q_matrix[which(q_matrix==0)]=n
    # Step 3: glp_matrix
    glp_matrix <- (2*q_matrix-1)/(2*n)
    
    return(glp_matrix)
  }
}

# Calculate log-likelihood function
LnML <- function(sample,theta){
  w=theta[1]
  m1=theta[2]
  v1=theta[3]
  m2=theta[4]
  v2=theta[5]
  x=sample
  d=w*(1/(sqrt(2*pi*v1))) * exp((-(x-m1)^2)/(2*v1))+
    (1-w)*(1/(sqrt(2*pi*v2))) * exp((-(x-m2)^2)/(2*v2))
  return(sum(log(d)))
}

# Location-Scale Mixtures
SNTO_Loglik <- function(sample,plot.lb=-5,plot.ub=10,
                        sizeN = "large",cs.kn.pop=F,cs.bound=0.05,
                        ep=1e-10,print.log=F,qua=T){
  # Revised data
  if(qua){
    sample <- quantile_estimator(sample)
  }
  # Mixn.oral=mixn(thetas[1],thetas[2],thetas[3],thetas[4],thetas[5])
  THETA0 <- Initial_thetas(sample)
  # loglik0<-LnML(sample,THETA0)# define n1;h and n2;h for glp
  if (sizeN == "small"){
    n1<-15019
    h1<-c(1,10641,2640,6710,784)
    n2<-1543
    h2<-c(1,58,278,694,134)
  }else if (sizeN == "normal"){
    n1<-71053
    h1<-c(1,33755,65170,12740,6878)
    n2<-8191
    h2<-c(1,1386,4302,7715,3735)
  }else if (sizeN == "large"){
    n1<-374181
    h1<-c(1,343867,255381,310881,115892)
    n2<-33139
    h2<-c(1,32133,17866,21281,32247)
  }
  # set up initial half edge vector by two ways
  if(cs.kn.pop==T){
    AS0 <- thetas - cs.bound
    BS0 <- thetas + cs.bound
    AS0[1]<-ifelse(AS0[1]<0,0,AS0[1])
    AS0[3]<-ifelse(AS0[3]<0,0.01,AS0[3])
    AS0[5]<-ifelse(AS0[5]<0,0.01,AS0[5])
    BS0[1]<-ifelse(BS0[1]>1,1,BS0[1])
    CS <- (BS0-AS0)/2
  }else{
    AS0 <- THETA0 - cs.bound
    BS0 <- THETA0 + cs.bound
    AS0[1]<-ifelse(AS0[1]<0,0,AS0[1])
    AS0[3]<-ifelse(AS0[3]<=0,0.2,AS0[3])
    AS0[5]<-ifelse(AS0[5]<=0,0.2,AS0[5])
    BS0[1]<-ifelse(BS0[1]>1,1,BS0[1])
    CS <- (BS0-AS0)/2
  }
  # Create interval data frame
  interval_df <- data.frame(
    Parameter = paste0("param_", 1:length(THETA0)),
    Lower_Bound = AS0,
    Upper_Bound = BS0
  )
  # Set up initial D and P(YS0)
  #AS0 <- c(0.455,0.9,0.95,4,0.95)
  #BS0 <- c(0.545,1.2,1.05,4.8,1.05)
  #AS0 <- as.matrix(THETA0-CS)
  #BS0 <- as.matrix(THETA0+CS)
  #print(data.frame(THETA0=THETA0,CS=CS,AS0=AS0,thetas=thetas,BS0=BS0))
  # D<-data.frame(CS=CS,AS0=AS0,thetas=thetas,BS0=BS0)
  GLP0<-NT_Net_glp(n=n1,h=h1,p=5)
  YS0<-t(apply(GLP0,1,function(x){AS0+(BS0-AS0)*x}))
  # Calculate M for each P(YS0)
  M0<-apply(YS0,1,function(x){LnML(sample,x)})
  # Select the maximum M and THETA 
  maxM <- max(M0)
  THETA<-YS0[which(M0==(max(M0))),]
  # Preparation for iteration
  iter=0
  GLP<-NT_Net_glp(n=n2,h=h2,p=5)
  while(max(CS)>ep){
    iter=iter+1
    # Reduce the edge vector by half and begin the iteration
    CS <- CS/2
    # Update A(lower bound) and B(upper bound) of D
    AS <- apply(cbind(THETA-CS,AS0),1,max)
    BS <- apply(cbind(THETA+CS,BS0),1,min)
    AS <- AS+(THETA-AS)/2
    BS <- BS-(BS-THETA)/2
    YS<-t(apply(GLP,1,function(x){AS+(BS-AS)*x}))
    M.temp<-apply(YS,1,function(x){LnML(sample,x)})
    # If the larger M found, update the THETA
    if (maxM<max(M.temp)){
      M <- M.temp
      maxM <- max(M)
      THETA <- YS[which(M==(max(M)))[1],]
    }
    if(print.log){
      print(data.frame(iteration=iter,thetas=thetas,THETA=THETA,Loglik=maxM))
    }
  }
  # Generate two output tables
  # Mixn.est.snto=mixn(THETA[1],THETA[2],THETA[3],THETA[4],THETA[5])
  #Mixn.est.ini=mixn(THETA0[1],THETA0[2],THETA0[3],THETA0[4],THETA0[5])
  tb_parm <-as.data.frame(t(THETA))
  names(tb_parm)=c("alpha","mu1","var1","mu2","var2")
  print(paste("number of iteration is", iter))
  return(list(est.parms=tb_parm))
}

# Convert 13-parameter THETA to 11-parameter decomposition form (add eps to avoid singular values)
convert_THETA_to_decomposition <- function(THETA_13, eps = 1e-8) {
  # Extract parameters
  mu1_1 <- THETA_13[1]
  mu1_2 <- THETA_13[2] 
  mu2_1 <- THETA_13[3]
  mu2_2 <- THETA_13[4]
  
  # First covariance matrix
  sigma_1_11 <- THETA_13[5]
  sigma_1_12 <- THETA_13[6]
  sigma_1_22 <- THETA_13[8]
  
  # Second covariance matrix
  sigma_2_11 <- THETA_13[9]
  sigma_2_12 <- THETA_13[10]
  sigma_2_22 <- THETA_13[12]
  
  alpha <- THETA_13[13]
  
  # Decompose first covariance matrix (add eps to avoid singular values)
  sigma_1_11 <- pmax(sigma_1_11, eps)
  d_11_1 <- log(sigma_1_11)
  u_21_1 <- sigma_1_12 / sigma_1_11
  d_22_1 <- log(pmax(sigma_1_22 - u_21_1^2 * sigma_1_11, eps))
  
  # Decompose second covariance matrix (add eps to avoid singular values)
  sigma_2_11 <- pmax(sigma_2_11, eps)
  d_11_2 <- log(sigma_2_11)
  u_21_2 <- sigma_2_12 / sigma_2_11
  d_22_2 <- log(pmax(sigma_2_22 - u_21_2^2 * sigma_2_11, eps))
  
  # Return 11-parameter vector
  return(c(mu1_1, mu1_2, mu2_1, mu2_2, 
           d_11_1, u_21_1, d_22_1, d_11_2, u_21_2, d_22_2, alpha))
}

# Simplified likelihood function based on covariance decomposition
MLE_2D_Mixture_Correct <- function(data,theta) {
  # Parameter description:
  # theta: parameter vector [mu1_x, mu1_y, mu2_x, mu2_y, 
  #          d11_1, u_1, d22_1,   # First distribution covariance parameters
  #          d11_2, u_2, d22_2,   # Second distribution covariance parameters
  #          alpha]                # Mixture weight
  
  # Extract means
  mu1 <- theta[1:2]
  mu2 <- theta[3:4]
  
  # Function to build covariance matrix
  build_Sigma <- function(d11, u, d22) {
    # Parameterize covariance matrix
    # Σ = [ exp(d11)               , exp(d11)*u  ]
    #     [ exp(d11)*u    , u^2*exp(d11) + exp(d22) ]
    matrix(c(exp(d11), 
             u * exp(d11),
             u * exp(d11),
             u^2 * exp(d11) + exp(d22)), 
           nrow = 2)
  }
  
  # Build covariance matrices for both distributions
  Sigma1 <- build_Sigma(theta[5], theta[6], theta[7])
  Sigma2 <- build_Sigma(theta[8], theta[9], theta[10])
  
  # Mixture weight 
  alpha <- theta[11]
  
  # Calculate densities
  dens1 <- dmvnorm(data, mu1, Sigma1)
  dens2 <- dmvnorm(data, mu2, Sigma2)
  
  # Calculate mixture density
  mixture_dens <- alpha * dens1 + (1 - alpha) * dens2
  
  # Return log likelihood
  sum(log(mixture_dens)) 
}

# Location-Scale Mixtures (2D Version)
SNTO_Loglik_2D <- function(Initial= THETA_initial,sample, plot.lb = -5, plot.ub = 10,
                           sizeN = "normal", cs.kn.pop = FALSE, cs.bound = 0.05,
                           ep = 1e-10, print.log = FALSE) {
  
  # Inverse transformation function: convert from 11 parameters to 13 parameters
  convert_decomposition_to_THETA <- function(THETA_11) {
    mu1_1 <- THETA_11[1]
    mu1_2 <- THETA_11[2]
    mu2_1 <- THETA_11[3]
    mu2_2 <- THETA_11[4]
    d_11_1 <- THETA_11[5]
    u_21_1 <- THETA_11[6]
    d_22_1 <- THETA_11[7]
    d_11_2 <- THETA_11[8]
    u_21_2 <- THETA_11[9]
    d_22_2 <- THETA_11[10]
    alpha <- THETA_11[11]
    
    # Reconstruct covariance matrices
    sigma_1_11 <- exp(d_11_1)
    sigma_1_12 <- exp(d_11_1) * u_21_1
    sigma_1_22 <- u_21_1^2 * exp(d_11_1) + exp(d_22_1)
    
    sigma_2_11 <- exp(d_11_2)
    sigma_2_12 <- exp(d_11_2) * u_21_2
    sigma_2_22 <- u_21_2^2 * exp(d_11_2) + exp(d_22_2)
    
    
    # Return 13 parameters
    return(c(mu1_1, mu1_2, mu2_1, mu2_2,
             sigma_1_11, sigma_1_12, sigma_1_12, sigma_1_22,
             sigma_2_11, sigma_2_12, sigma_2_12, sigma_2_22,
             alpha))
  }
  # Get initial parameters for 2D mixture (11 parameters)
  THETA0 <- Initial
  if (sizeN == "normal") {
   n1 <- 698047
   h1 <- c(1, 685041, 646274, 582461, 494796, 384914, 254860, 107051, 642292, 467527, 284044)
   n2 <- 58358
   h2 <- c(1, 57271, 54030, 48695, 41366, 32180, 21307, 8950, 53697, 39086, 23747)
  } 
    else if (sizeN == "middle") {
   n1 <- 1243423
   h1 <- c(1,1228845,1185282,1113244,1013577,887449,736338,562016,366527,152163,1164860)
   n2 <- 698047
   h2 <- c(1, 685041, 646274, 582461, 494796, 384914, 254860, 107051, 642292, 467527, 284044)
  } 
    else if (sizeN == "large") {
   n1 <- 7494007
   h1 <- c(1, 7354408,6838211,6253169,5312043,4132365,2736109,1149286,6895461,5019180,3049402)
   n2 <- 297974
   h2 <- c(1, 294481, 284041,266778,242894,212668,176456,134682,87835,36464,279147)
  } 
  
  # Set up initial bounds
  AS0 <- THETA0 - cs.bound
  BS0 <- THETA0 + cs.bound
  
  # Apply constraints to parameters
  # 1. Mixing weight (alpha) in [0,1]
  AS0[11] <- ifelse(AS0[11] < 0, 0, AS0[11])  # alpha_logit constraint
  BS0[11] <- ifelse(BS0[11] > 1, 1, BS0[11])
  # Compute initial parameter ranges
  CS <- (BS0 - AS0) / 2
  
  # Generate initial grid points
  GLP0 <- NT_Net_glp(n = n1, h = h1, p = 11)
  YS0 <- t(apply(GLP0, 1, function(x) {
    AS0 + (BS0 - AS0) * x
  }))
  
  # Evaluate at grid points
  M0 <- apply(YS0, 1, function(x) {
    MLE_2D_Mixture_Correct(sample, x)  # Use new 2D likelihood function
  })
  
  # Find best parameter set
  maxM <- max(M0)
  THETA <- YS0[which.max(M0), ]
  
  # Iterative refinement
  iter <- 0
  GLP <- NT_Net_glp(n = n2, h = h2, p = 11)
  
  # Refine until convergence
  while(max(CS) > ep) {
    iter <- iter + 1
    CS <- CS / 2  # Halve search radius
    
    # Update parameter bounds
    AS <- pmax(THETA - CS, AS0)
    BS <- pmin(THETA + CS, BS0)
    AS <- AS + (THETA - AS) / 2
    BS <- BS - (BS - THETA) / 2
    
    # Generate new grid points
    YS <- t(apply(GLP, 1, function(x) {
      AS + (BS - AS) * x
    }))
    
    # Evaluate at new grid points
    M_temp <- apply(YS, 1, function(x) {
      MLE_2D_Mixture_Correct(sample, x)  # Use new 2D likelihood function
    })
    
    # Update if better solution found
    if (max(M_temp) > maxM) {
      maxM <- max(M_temp)
      THETA <- YS[which.max(M_temp), ]
    }
    
    # Print progress if requested
    if(print.log) {
      cat(sprintf("Iteration %d, Log-likelihood: %.4f\n", iter, maxM))
    }
  }
  
  # Convert 11-parameter THETA to 13 parameters
  params <- convert_decomposition_to_THETA(THETA)
  
  # Extract components and build matrices
  mu1 <- c(params[1], params[2])        # mu1_x, mu1_y
  mu2 <- c(params[3], params[4])        # mu2_x, mu2_y
  Sigma1 <- matrix(c(params[5], params[6], params[7], params[8]), nrow = 2)  # sigma_1_11, sigma_1_12, sigma_1_21, sigma_1_22
  Sigma2 <- matrix(c(params[9], params[10], params[11], params[12]), nrow = 2) # sigma_2_11, sigma_2_12, sigma_2_21, sigma_2_22
  alpha <- params[13]                   # mixture weight
  
  # Format results in a named data frame
  result_df <- data.frame(
    Parameter = c("mu1_1", "mu1_2", "mu2_1", "mu2_2",
                  "sigma_1_11", "sigma_1_12", "sigma_1_21", "sigma_1_22",
                  "sigma_2_11", "sigma_2_12", "sigma_2_21", "sigma_2_22",
                  "alpha"),
    Value = params,
    stringsAsFactors = FALSE
  )
  
# In the SNTO_2D function's end, modify the return value section
return(list(
  # Keep other original return values unchanged
  parameters = result_df,
  mu1 = mu1,
  mu2 = mu2,
  Sigma1 = Sigma1,
  Sigma2 = Sigma2,
  alpha = alpha,
  log_likelihood = maxM,
  
  # New addition: THETA vector consistent with true_params format
  THETA = c(
    mu1[1], mu1[2],           # Two components of mu1
    mu2[1], mu2[2],           # Two components of mu2
    as.vector(Sigma1),        # 4 elements of Sigma1 (column-wise expansion)
    as.vector(Sigma2),        # 4 elements of Sigma2 (column-wise expansion)
    alpha                     # Mixture weight
  )
))
}

### Functions for Accuracy Measures
#####################################################################################
# L2 Distance for the density comparison (Mixture of Two Normals) using Monte Carlo integration
l2distancef<- function(thetas, THETA, n_samples = 10000) {
  # Unpack the true parameters from the true vector `thetas`
  mu1_true <- thetas[1:2]   # mu1_real(1,1), mu1_real(1,2)
  mu2_true <- thetas[3:4]   # mu2_real(1,1), mu2_real(1,2)
  Sigma1_true <- matrix(thetas[5:8], nrow = 2, byrow = TRUE)   # Sigma1_real(1,1), Sigma1_real(1,2), Sigma1_real(2,1), Sigma1_real(2,2)
  Sigma2_true <- matrix(thetas[9:12], nrow = 2, byrow = TRUE)  # Sigma2_real(1,1), Sigma2_real(1,2), Sigma2_real(2,1), Sigma2_real(2,2)
  pi1_true <- thetas[13]    # pi1_real
  
  # Unpack the estimated parameters from the estimated vector `THETA`
  mu1_est <- THETA[1:2]   # Predicted mu1
  mu2_est <- THETA[3:4]   # Predicted mu2
  Sigma1_est <- matrix(THETA[5:8], nrow = 2, byrow = TRUE)   # Predicted Sigma1
  Sigma2_est <- matrix(THETA[9:12], nrow = 2, byrow = TRUE)  # Predicted Sigma2
  pi1_est <- THETA[13]    # Predicted pi1
  
  # Generate random samples for the 2D space (x1, x2) for Monte Carlo integration
  x_samples <- matrix(runif(n_samples * 2, -20, 20), ncol = 2)  # Generate 2D samples between -10 and 10
  
  # Calculate the mixture density for the true model at each sample point
  density_true_samples <- pi1_true * apply(x_samples, 1, function(x) dmvnorm(x, mean = mu1_true, sigma = Sigma1_true)) +
    (1 - pi1_true) * apply(x_samples, 1, function(x) dmvnorm(x, mean = mu2_true, sigma = Sigma2_true))
  
  # Calculate the mixture density for the estimated model at each sample point
  density_est_samples <- pi1_est * apply(x_samples, 1, function(x) dmvnorm(x, mean = mu1_est, sigma = Sigma1_est)) +
    (1 - pi1_est) * apply(x_samples, 1, function(x) dmvnorm(x, mean = mu2_est, sigma = Sigma2_est))
  
  # Compute the squared difference between the true and estimated densities for each sample
  squared_diff <- (density_true_samples - density_est_samples)^2
  
  # Compute the average squared difference (integral approximation)
  avg_squared_diff <- mean(squared_diff)
  
  # Return the square root of the average squared difference (L2 distance)
  return(sqrt(avg_squared_diff))
}

# L2 Distance for the cumulative distribution comparison (Mixture of Two Normals) using Monte Carlo integration
l2distanceF <- function(thetas, THETA, n_samples = 10000) {
  # Unpack the true parameters from the true vector `thetas`
  mu1_true <- thetas[1:2]   # mu1_real(1,1), mu1_real(1,2)
  mu2_true <- thetas[3:4]   # mu2_real(1,1), mu2_real(1,2)
  Sigma1_true <- matrix(thetas[5:8], nrow = 2, byrow = TRUE)   # Sigma1_real(1,1), Sigma1_real(1,2), Sigma1_real(2,1), Sigma1_real(2,2)
  Sigma2_true <- matrix(thetas[9:12], nrow = 2, byrow = TRUE)  # Sigma2_real(1,1), Sigma2_real(1,2), Sigma2_real(2,1), Sigma2_real(2,2)
  pi1_true <- thetas[13]    # pi1_real
  
  # Unpack the estimated parameters from the estimated vector `THETA`
  mu1_est <- THETA[1:2]   # Predicted mu1
  mu2_est <- THETA[3:4]   # Predicted mu2
  Sigma1_est <- matrix(THETA[5:8], nrow = 2, byrow = TRUE)   # Predicted Sigma1
  Sigma2_est <- matrix(THETA[9:12], nrow = 2, byrow = TRUE)  # Predicted Sigma2
  pi1_est <- THETA[13]    # Predicted pi1
  
  # Generate random samples for the 2D space (x1, x2) for Monte Carlo integration
  x_samples <- matrix(runif(n_samples * 2, -20, 20), ncol = 2)  # Generate 2D samples between -50 and 50
  
  # Function to compute the CDF of a 2D normal mixture
  cdf_function <- function(x, mu1, Sigma1, mu2, Sigma2, pi1) {
    # Compute the CDF for the mixture of two normals at point x
    cdf1 <- pmvnorm(lower = rep(-Inf, 2), upper = x, mean = mu1, sigma = Sigma1)[1]
    cdf2 <- pmvnorm(lower = rep(-Inf, 2), upper = x, mean = mu2, sigma = Sigma2)[1]
    
    # Mixture CDF: weighted sum of the two CDFs
    return(pi1 * cdf1 + (1 - pi1) * cdf2)
  }
  
  # Compute the mixture CDF for the true model at each sample point
  cdf_true_samples <- apply(x_samples, 1, function(x) cdf_function(x, mu1_true, Sigma1_true, mu2_true, Sigma2_true, pi1_true))
  
  # Compute the mixture CDF for the estimated model at each sample point
  cdf_est_samples <- apply(x_samples, 1, function(x) cdf_function(x, mu1_est, Sigma1_est, mu2_est, Sigma2_est, pi1_est))
  
  # Compute the squared difference between the CDFs for each sample
  squared_diff <- (cdf_true_samples - cdf_est_samples)^2
  
  # Compute the average squared difference (integral approximation)
  avg_squared_diff <- mean(squared_diff)
  
  # Return the square root of the average squared difference (L2 distance)
  return(sqrt(avg_squared_diff))
}

# KL Divergence using sampling from the true distribution
KL <- function(thetas, THETA, n_samples = 10000) {
  # Unpack the true parameters
  mu1_true <- thetas[1:2]
  mu2_true <- thetas[3:4]
  Sigma1_true <- matrix(thetas[5:8], nrow = 2, byrow = TRUE)
  Sigma2_true <- matrix(thetas[9:12], nrow = 2, byrow = TRUE)
  pi1_true <- thetas[13]
  
  # Unpack the estimated parameters (handling the duplicated off-diagonal elements)
  mu1_est <- THETA[1:2]
  mu2_est <- THETA[3:4]
  Sigma1_est <- matrix(c(THETA[5], THETA[6], THETA[6], THETA[8]), nrow = 2)
  Sigma2_est <- matrix(c(THETA[9], THETA[10], THETA[10], THETA[12]), nrow = 2)
  pi1_est <- THETA[13]
  
  # Sample from the true distribution
  n1_samples <- round(n_samples * pi1_true)
  n2_samples <- n_samples - n1_samples
  
  # Generate samples from each component of the true distribution
  samples1 <- rmvnorm(n1_samples, mean = mu1_true, sigma = Sigma1_true)
  samples2 <- rmvnorm(n2_samples, mean = mu2_true, sigma = Sigma2_true)
  all_samples <- rbind(samples1, samples2)
  
  # Compute densities under true distribution
  p_true <- pi1_true * apply(all_samples, 1, function(x) dmvnorm(x, mean = mu1_true, sigma = Sigma1_true)) +
           (1 - pi1_true) * apply(all_samples, 1, function(x) dmvnorm(x, mean = mu2_true, sigma = Sigma2_true))
  
  # Compute densities under estimated distribution
  p_est <- pi1_est * apply(all_samples, 1, function(x) dmvnorm(x, mean = mu1_est, sigma = Sigma1_est)) +
          (1 - pi1_est) * apply(all_samples, 1, function(x) dmvnorm(x, mean = mu2_est, sigma = Sigma2_est))

  # Compute KL divergence
  kl_divergence <- mean(log(p_true / p_est))
  
  return(kl_divergence)
}