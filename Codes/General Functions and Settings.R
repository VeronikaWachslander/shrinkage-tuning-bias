################################################################################
################################################################################
###                     General Functions and Settings                       ###
################################################################################
################################################################################

### Definition of Required Functions 

################################################################################

### Function to Load Required Packages

fun_LoadLibraries <- function(){
  library(glmnetUtils)
  library(OpenMx)
  library(MASS)
  library(ggplot2)
  library(tidyverse)
  library(cvTools)
  require(caret)
}

set.seed(1234)

### Function to Generate Scenario Parameters (Initial Settings; Later Modified)

fun_GenerateScenario <- function(){
  scenario <- list();
  scenario$corr <- 0.6 # Correlation of FC
  scenario$variance <- 2.0
  scenario$populationsize <- 20000
  scenario$train <- n <- 20
  scenario$test <- scenario$populationsize - scenario$train
  scenario$repetitions = 20
  scenario$k <- 4 # Number of FC
  scenario$tau_grid <- expand.grid(replicate(scenario$k, c(0,1), simplify=FALSE))
  scenario$tau_grid <- cbind(scenario$tau_grid, scenario$tau_grid)
  scenario$I <- 1 # Number of Runs
  scenario$no_steps <- 101 # Define Shrinkage Steps
  scenario$seed <- runif(scenario$repetitions, 1, 1000000) # Reproducibility
  scenario$source_vector <- rep(1, scenario$k)
  scenario$target_vector <- scenario$source_vector; 
  return (scenario)
}

### Functions to Generate Covariance Matrix Sigma 

# Sigma with Predefined Variance Vector

fun_MakeSigma <- function(this.k, this.corr, this.varvec){
  variance_vector <- rep(NA, this.k)
  for(v in 1:this.k){
    variance_vector[v] <- this.varvec[v]
  }
  mu=rep(0,this.k); # Means of Individual Distributions
  TempMatrix <- matrix(rep(this.corr,this.k^2), ncol = this.k, byrow = TRUE) 
  TempMatrix <- TempMatrix - diag(rep(this.corr, this.k)) + diag(variance_vector);
  Sigma <- TempMatrix + diag(rep(1, this.k)) - diag( as.vector(t(diag2vec(TempMatrix)) ))
  for (v in 1:this.k){
    Sigma[v,] <- (Sigma[v,]*sqrt(variance_vector[v])); Sigma[,v] <- (Sigma[,v]*sqrt(variance_vector[v]) )
    Sigma[v,v] <-variance_vector[v]
  }
  return (Sigma)
} 

# Sigma with Linearly Increasing Variances

fun_MakeSigmaLin  <- function(this.k, this.corr, this.maxvar){
  variance_vector <- rep(NA, this.k)
  variance_vector[1] <- 1.00 # 1-st FC
  variance_vector[this.k] <- this.maxvar # k-th FC
  #  Variances grow linearly from FC 1 to FC k 
  for (v in 1:(this.k-2)) { 
    variance_vector[v+1] <- ( ((this.k-v-1)*variance_vector[1]) + (v*variance_vector[this.k]) )/ (this.k-1)
  }
  mu=rep(0,this.k); # Means of Individual Distributions
  TempMatrix <- matrix(rep(this.corr,this.k^2), ncol = this.k, byrow = TRUE) 
  TempMatrix <- TempMatrix - diag(rep(this.corr, this.k)) + diag(variance_vector);
  Sigma <- TempMatrix + diag(rep(1, this.k)) - diag( as.vector(t(diag2vec(TempMatrix)) ))
  for (v in 1:this.k){
    Sigma[v,] <- (Sigma[v,]*sqrt(variance_vector[v])); Sigma[,v] <- (Sigma[,v]*sqrt(variance_vector[v]) )
    Sigma[v,v] <- variance_vector[v]
  }
  return (Sigma)
} 


# Sigma with Quadratically Increasing Variances

fun_MakeSigmaExp  <- function(this.k, this.corr, this.maxvar){
  variance_vector <- rep(NA, this.k)
  variance_vector[1] <- 1.00 # 1-st FC
  variance_vector[this.k] <- this.maxvar # k-th FC
  #  Variances grow quadratically from FC 1 to FC k
  for (v in 1:(this.k-2)) { 
    variance_vector[v+1] <- (( ((this.k-v-1)*variance_vector[1]) + (v*sqrt(variance_vector[this.k])) )/(this.k-1))^2
  }
  mu=rep(0,this.k); # Means of Individual Distributions
  TempMatrix <- matrix(rep(this.corr,this.k^2), ncol = this.k, byrow = TRUE) 
  TempMatrix <- TempMatrix - diag(rep(this.corr, this.k)) + diag(variance_vector);
  Sigma <- TempMatrix + diag(rep(1, this.k)) - diag( as.vector(t(diag2vec(TempMatrix)) ))
  for (v in 1:this.k){
    Sigma[v,] <- (Sigma[v,]*sqrt(variance_vector[v])); Sigma[,v] <- (Sigma[,v]*sqrt(variance_vector[v]) )
    Sigma[v,v] <- variance_vector[v]
  }
  return (Sigma)
} 


# Sigma with Linearly Increasing Variances and Deviating Correlations

fun_MakeSigmaLinDev <- function(this.k, this.corr, this.maxvar){
  variance_vector <- rep(NA, this.k)
  variance_vector[1] <- 1.00 # 1-st FC
  variance_vector[this.k] <- this.maxvar # k-th FC
  #  Variances grow linearly from FC 1 to FC k 
  for (v in 1:(this.k-2)) { 
    variance_vector[v+1] <- ( ((this.k-v-1)*variance_vector[1]) + (v*variance_vector[this.k]) )/ (this.k-1)
  }
  
  mu=rep(0,this.k); # Means of Individual Distributions
  
  TempMatrix <- matrix(NA, ncol = this.k, nrow = this.k)
  
  for(i in 1:this.k){
    for(j in 1:this.k){
      if(i <= (this.k+1)/2 & j <= (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr-0.1)
      }else if(i > (this.k+1)/2 & j > (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr+0.1)
      }else{
        TempMatrix[i,j] <- this.corr
      }
    }
  }
  
  for(i in 1:this.k){
    for(j in 1:this.k){
      if(i == j){
        TempMatrix[i,j] <- 1
      }
    }
  } 
  Sigma <- TempMatrix
  
  for (v in 1:this.k){
    Sigma[v,] <- (Sigma[v,]*sqrt(variance_vector[v])); Sigma[,v] <- (Sigma[,v]*sqrt(variance_vector[v]) )
    Sigma[v,v] <-variance_vector[v]
  }
  return(Sigma)
}


# Sigma with Quadratically Increasing Variances and Deviating Correlations

fun_MakeSigmaExpDev <- function(this.k, this.corr, this.maxvar){
  variance_vector <- rep(NA, this.k)
  variance_vector[1] <- 1.00 # 1-st FC
  variance_vector[this.k] <- this.maxvar # k-th FC
  #  Variances grow quadratically from FC 1 to FC k
  for (v in 1:(this.k-2)) { 
    variance_vector[v+1] <- (( ((this.k-v-1)*variance_vector[1]) + (v*sqrt(variance_vector[this.k])) )/(this.k-1))^2
  }
  
  mu=rep(0,this.k); # Means of Individual Distributions
  
  TempMatrix <- matrix(NA, ncol = this.k, nrow = this.k)
  
  for(i in 1:this.k){
    for(j in 1:this.k){
      if(i <= (this.k+1)/2 & j <= (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr-0.1)
      }else if(i > (this.k+1)/2 & j > (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr+0.1)
      }else{
        TempMatrix[i,j] <- this.corr
      }
    }
  }
  
  for(i in 1:this.k){
    for(j in 1:this.k){
      if(i == j){
        TempMatrix[i,j] <- 1
      }
    }
  } 
  Sigma <- TempMatrix
  
  for (v in 1:this.k){
    Sigma[v,] <- (Sigma[v,]*sqrt(variance_vector[v])); Sigma[,v] <- (Sigma[,v]*sqrt(variance_vector[v]) )
    Sigma[v,v] <-variance_vector[v]
  }
  return(Sigma)
}


### Function to Generate Error Data

fun_MakeErrors <- function(this.Sigma, this.populationsize, this.n_train, this.n_test, this.seed){
  this.k <- ncol(this.Sigma)
  set.seed(this.seed)
  mu = rep(0, this.k)
  population_errors <- mvrnorm(this.populationsize, mu = mu, Sigma = this.Sigma, empirical = TRUE);# generate errors
  Sample_Indices <- sample(seq(1:this.populationsize), this.n_test + this.n_train, replace = FALSE);
  C_sample <- population_errors[Sample_Indices,]
  C_train <- C_sample[1:this.n_train,];
  C_test <- C_sample[(this.n_train + 1):(this.n_train + this.n_test),];
  returnlist <- list();
  returnlist[[1]] <- C_train;
  returnlist[[2]] <- C_test;
  returnlist[[3]] <- C_sample;
  names(returnlist) <- c("Training", "Test", "Population")
  return <- returnlist;
}


### Function to Compute Optimal Weights (OW)

fun_ComputeOW <- function(this.data, this.source_vector){
  this.source_vector_colToUse <- (this.source_vector*seq(1:length(this.source_vector)))
  this.source_vector_colToUse <- this.source_vector_colToUse[this.source_vector_colToUse>0]
  
  this.data_source <- this.data[, this.source_vector_colToUse]
  
  # Regress ek ~ di to Find Optimal Weights
  this.ks <- ncol(this.data_source)
  ErrorDiv <- array(dim = c(nrow(this.data_source),(this.ks-1))); 
  for (fc in 1:(this.ks-1)){
    ErrorDiv[,fc] <- this.data_source[,this.ks] - this.data_source[,fc];  
  }
  head(this.data_source)
  head(ErrorDiv)
  OW.model <- lm( this.data_source[,this.ks] ~ ErrorDiv[,]-1 )
  OW.weights <- c(OW.model$coefficients[1:(this.ks-1)], (1-sum(OW.model$coefficients[1:(this.ks-1)]))) 
  
  where_to_set_OW <- this.source_vector > 0 
  OW.position = 1;      
  OW.computed_nr = 1;      
  OW_weights_in_long_vector <- rep(0, length(this.source_vector))
  for (i in 1:length(this.source_vector)) {
    if (where_to_set_OW[i]){
      OW_weights_in_long_vector[i] <- OW.weights[OW.computed_nr];
      OW.computed_nr <- OW.computed_nr + 1;
    } 
  }   
  OW_weights_in_long_vector
  return(OW_weights_in_long_vector)
}


### Function to Calculate Mean Squared Error (MSE)

fun_DetermineMSE <- function(this.weights, this.data){
  this.weights <- as.matrix(this.weights) 
  if (ncol(this.weights) == 1){
    this.weights <- t(this.weights)
  }
  
  howMANYWeightVectors <- nrow(this.weights) 
  
  MSE <- array(NA, dim = howMANYWeightVectors)  
  data.weighted <- this.data
  
  if(is.vector(data.weighted) == TRUE){
    for (i in 1: howMANYWeightVectors ){        
      data.weighted <- this.weights[i,] * this.data
      MSE[i] <- (sum(data.weighted))^2
    }
  }else{
    for (i in 1: howMANYWeightVectors ){        
      for (j in 1:nrow(data.weighted)){ 
        data.weighted[j,] <- this.weights[i,] * this.data[j,]
      }
      MSE[i] <- mean((rowSums(data.weighted))^2)
    }
  }
  return(MSE)
}

### Function to Calculate Shrinkage Path (here from OW to EW for All Forecasters)

fun_ComputeTTSWeights <- function(this.OW, this.target_vector, this.steps) {
  this.kt <- length(this.target_vector);
  weight_matrix <- matrix(data = rep(0, (this.steps*this.kt )), ncol = this.kt)
  subAweights <- this.target_vector/sum(this.target_vector);
  
  for (step in 0: (nrow(weight_matrix) - 1)){
    weight_matrix[(step + 1),] <- 
      (
        (  (1-(step/(this.steps-1))) * this.OW)    
        +
          ( (step/(this.steps-1))* subAweights)  
      ) 
  }
  return(weight_matrix);
}


################################################################################

### Script               

################################################################################

### Load Libraries

fun_LoadLibraries() 

### Generate Scenario; return: scenario (list)

scenario <-fun_GenerateScenario() 
