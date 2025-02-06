################################################################################
################################################################################
###                  Generation of Analytical Datasets                       ###
################################################################################
################################################################################

### Generation of Synthetic Datasets for Different Parameter Value 

### Constellations (Scenarios) to Analyze Shrinkage Tuning Bias 

################################################################################


source("General Functions and Settings.R")


# Four Separate Procedures (Pairwise Constant vs. Differing Error Correlations
# and General Scenarios vs. Special Case n = 20, J = {12,15})
# (Special Case as for n = 20, J = {12,15} a 2-Fold CV is not possible)

# Do Each of the Procedures for All Combinations of Parameter Values 


### 1. Procedure (Pairwise Constant Correlations; General Scenarios)

# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (5,8,10,12,15) and Choosing Respective Number
# of Training Observations (Ensuring Sufficient Amount of
# Training Observations to Calculate OW) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.2,
# 1.5, 2, 4, 9)
# -> Resulting in 50 Separate Runs


scenario$k <- 5 # 5, 8, 10, 12, 15

scenario$train <- c(10,20,30,40,50,60,70,80,90,100,125,150,175,200) # for 5 FC
# scenario$train <- c(20,30,40,50,60,70,80,90,100,125,150,175,200) # for 8,10 FC
# scenario$train <- c(30,40,50,60,70,80,90,100,125,150,175,200) # for 12,15 FC

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # Pairwise Correlation


set.seed(11111111)

# 250 Repetitions of Each Parameter Constellation (Scenario)

for(running in 1:250){
  
  # Data Containers
  
  no_folds <- rep(NA, 5)
  
  test_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.2, 1.5, 2, 4, 9)
    
    Sigma <- fun_MakeSigmaLin(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    #Sigma <- fun_MakeSigmaExp(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character)
      
      no_folds <- c(2, 5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      no_folds_chr <- c(2, 5, 10, "LTO", "LOO")
      
      
      # Data Containers
      
      weights_on_test <- list()
      MSE_on_test <- list()
      
      
      # Generate New Error Data in Each Run
      
      scenario$test <- scenario$populationsize - scenario$train[ntrain]
      
      ErrorData  <- fun_MakeErrors(Sigma, scenario$populationsize, scenario$train[ntrain], scenario$test, runif(1, 1, 100000000))
      
      
      # Estimation of Correlation-Related Variables on Test Data
      
      true_cormat <- cor(ErrorData$Test)
      
      diag(true_cormat) <- NA 
      
      true_mean_corr <- mean(true_cormat, na.rm = TRUE)
      
      true_median_corr <- median(true_cormat, na.rm = TRUE)
      
      true_min_corr <- min(true_cormat, na.rm = TRUE)
      
      true_max_corr <- max(true_cormat, na.rm = TRUE)
      
      true_low_quantile_corr <- quantile(true_cormat, probs = 0.2, na.rm = TRUE)
      
      true_high_quantile_corr <- quantile(true_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- quantile(train_cormat, probs = 0.2, na.rm = TRUE)
      
      high_quantile_corr_est <- quantile(train_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Variance-Related Variables on Test Data
      
      true_Sigma <- 1/scenario$test * t(ErrorData$Test) %*% ErrorData$Test
      
      true_av_var <- (max(diag(true_Sigma))-min(diag(true_Sigma)))/(scenario$k-1)
      
      true_mean_var <- mean(diag(true_Sigma))
      
      true_min_var <- min(diag(true_Sigma))
      
      true_max_var <- max(diag(true_Sigma))
      
      
      # Estimation of Variance-Related Variables on Training Data
      
      Sigma_est <- 1/scenario$train[ntrain] * t(ErrorData$Training) %*% ErrorData$Training
      
      av_var_est <- (max(diag(Sigma_est))-min(diag(Sigma_est)))/(scenario$k-1)
      
      mean_var_est <- mean(diag(Sigma_est))
      
      min_var_est <- min(diag(Sigma_est))
      
      max_var_est <- max(diag(Sigma_est))
      
      
      
      # Calculate OW on Full Training Set and Determine Shrinkage Path to EW 
      
      OW_emp <- fun_ComputeOW(ErrorData$Training, scenario$source_vector) 
      
      weights_on_test <- fun_ComputeTTSWeights(OW_emp, target_vector, granularity)  
      
      # Calculate MSE Values for Shrinkage Path on Test Set
      
      MSE_on_test <- fun_DetermineMSE(this.weights = weights_on_test, this.data = ErrorData$Test)  
      
      
      # Data Containers
      
      cv_val_result_per_fold_agg <- list()
      
      
      for(folds in 1:length(no_folds)){
        
        this.no_folds <- no_folds[folds]
        
        
        # Split Training Data in Folds for CV
        
        flds <- cvFolds(nrow(ErrorData$Training), K = this.no_folds)
        
        # Data Containers
        
        cv_weight.Candidates <- list()
        cv_val_result_per_fold <- list()
        
        
        for (this.fold in 1: this.no_folds){
          
          # Split Folds in Validation and Train (Calibration) Set for CV
          
          cv_this.val   <- ErrorData$Training[flds$subsets[flds$which == this.fold], ]   
          cv_this.train <- ErrorData$Training[flds$subsets[flds$which != this.fold], ]
          
          
          # Calculate OW on CV-Training Set and Determine Shrinkage Path to EW
          
          cv_this.OW <- fun_ComputeOW(this.data = cv_this.train, this.source_vector = source_vector) 
          
          cv_weight.Candidates[[this.fold]] <- fun_ComputeTTSWeights(cv_this.OW , target_vector, granularity)
          
          # Calculate MSE Values for Shrinkage Path on Respective Validation Set
          
          cv_val_result_per_fold[[this.fold]] <- fun_DetermineMSE(this.weights = cv_weight.Candidates[[this.fold]], this.data = cv_this.val) 
          
        }
        
        # Average MSE Values over Folds
        
        cv_val_result_per_fold_agg[[folds]] <- Reduce("+", cv_val_result_per_fold) / length(cv_val_result_per_fold) 
        
        # Identify CV-optimal Shrinkage Level 
        
        cv_lambda_agg[ntrain,folds,corr] <- which.min(cv_val_result_per_fold_agg[[folds]])-1
        
        # Identify Resulting MSE Value on Test Data with CV-optimal Shrinkage Level
        
        MSE_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[cv_lambda_agg[ntrain,folds,corr]+1]
        
      }
      
      
      test_av_var[ntrain,corr]<- true_av_var
      
      train_av_var[ntrain,corr] <- av_var_est
      
      test_mean_var[ntrain,corr] <-  true_mean_var
      test_min_var[ntrain,corr] <- true_min_var
      test_max_var[ntrain,corr] <-  true_max_var
      
      train_mean_var[ntrain,corr] <- mean_var_est
      train_min_var[ntrain,corr] <- min_var_est
      train_max_var[ntrain,corr] <- max_var_est
      
      test_mean_corr[ntrain,corr] <- true_mean_corr
      test_median_corr[ntrain,corr] <- true_median_corr
      test_min_corr[ntrain,corr] <- true_min_corr
      test_max_corr[ntrain,corr] <- true_max_corr
      test_low_quantile_corr[ntrain,corr] <- true_low_quantile_corr
      test_high_quantile_corr[ntrain,corr] <- true_high_quantile_corr
      
      train_mean_corr[ntrain,corr] <- mean_corr_est
      train_median_corr[ntrain,corr] <- median_corr_est
      train_min_corr[ntrain,corr] <- min_corr_est
      train_max_corr[ntrain,corr] <- max_corr_est
      train_low_quantile_corr[ntrain,corr] <- low_quantile_corr_est
      train_high_quantile_corr[ntrain,corr] <- high_quantile_corr_est
      
      
      # Identify Truly Optimal Shrinkage Level and Resulting MSE on Test Data
      
      which_min_MSE_on_test[ntrain,corr] <- which.min(MSE_on_test)-1
      
      min_MSE_on_test[ntrain,corr] <- min(MSE_on_test)
      
      
      # MSE for OW on Test Data
      
      OW_MSE_on_test[ntrain,corr] <- MSE_on_test[1]
      
      # MSE for EW on Test Data
      
      EW_MSE_on_test[ntrain,corr] <- MSE_on_test[101]
      
      
      print(c("finished ntrain =", ntrain))
    }
    print(c("finished corr =", corr))
  }
  
  
  # Organize and Store Results 
  
  analysis_data_corr0.1 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                               test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                               train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                               test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                               test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                               train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                               train_low_quantile_corr[,1], train_high_quantile_corr[,1], cv_lambda_agg[,,1], MSE_cv_on_test_agg[,,1],
                                               which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                               OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(analysis_data_corr0.1) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.2 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                               test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                               train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                               test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                               test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                               train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                               train_low_quantile_corr[,2], train_high_quantile_corr[,2], cv_lambda_agg[,,2], MSE_cv_on_test_agg[,,2],
                                               which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                               OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(analysis_data_corr0.2) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.3 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                               test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                               train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                               test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                               test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                               train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                               train_low_quantile_corr[,3], train_high_quantile_corr[,3], cv_lambda_agg[,,3], MSE_cv_on_test_agg[,,3],
                                               which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                               OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(analysis_data_corr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.4 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                               test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                               train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                               test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                               test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                               train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                               train_low_quantile_corr[,4], train_high_quantile_corr[,4], cv_lambda_agg[,,4], MSE_cv_on_test_agg[,,4],
                                               which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                               OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(analysis_data_corr0.4) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.5 <- as.data.frame(cbind(rep(scenario$corr[5],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,5], train_av_var[,5], 
                                               test_mean_var[,5], test_min_var[,5], test_max_var[,5],
                                               train_mean_var[,5], train_min_var[,5], train_max_var[,5],
                                               test_mean_corr[,5], test_median_corr[,5], test_min_corr[,5], test_max_corr[,5],
                                               test_low_quantile_corr[,5], test_high_quantile_corr[,5], 
                                               train_mean_corr[,5], train_median_corr[,5], train_min_corr[,5], train_max_corr[,5],
                                               train_low_quantile_corr[,5], train_high_quantile_corr[,5], cv_lambda_agg[,,5], MSE_cv_on_test_agg[,,5],
                                               which_min_MSE_on_test[,5], min_MSE_on_test[,5],
                                               OW_MSE_on_test[,5], EW_MSE_on_test[,5]))
  
  colnames(analysis_data_corr0.5) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.6 <- as.data.frame(cbind(rep(scenario$corr[6],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,6], train_av_var[,6],
                                               test_mean_var[,6], test_min_var[,6], test_max_var[,6],
                                               train_mean_var[,6], train_min_var[,6], train_max_var[,6],
                                               test_mean_corr[,6], test_median_corr[,6], test_min_corr[,6], test_max_corr[,6],
                                               test_low_quantile_corr[,6], test_high_quantile_corr[,6], 
                                               train_mean_corr[,6], train_median_corr[,6], train_min_corr[,6], train_max_corr[,6],
                                               train_low_quantile_corr[,6], train_high_quantile_corr[,6], cv_lambda_agg[,,6], MSE_cv_on_test_agg[,,6],
                                               which_min_MSE_on_test[,6], min_MSE_on_test[,6],
                                               OW_MSE_on_test[,6], EW_MSE_on_test[,6]))
  
  colnames(analysis_data_corr0.6) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                        "test_mean_var", "test_min_var", "test_max_var",
                                        "train_mean_var", "train_min_var", "train_max_var",
                                        "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                        "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                        "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.7 <- as.data.frame(cbind(rep(scenario$corr[7],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,7], train_av_var[,7], 
                                               test_mean_var[,7], test_min_var[,7], test_max_var[,7],
                                               train_mean_var[,7], train_min_var[,7], train_max_var[,7],
                                               test_mean_corr[,7], test_median_corr[,7], test_min_corr[,7], test_max_corr[,7],
                                               test_low_quantile_corr[,7], test_high_quantile_corr[,7], 
                                               train_mean_corr[,7], train_median_corr[,7], train_min_corr[,7], train_max_corr[,7],
                                               train_low_quantile_corr[,7], train_high_quantile_corr[,7], cv_lambda_agg[,,7], MSE_cv_on_test_agg[,,7],
                                               which_min_MSE_on_test[,7], min_MSE_on_test[,7],
                                               OW_MSE_on_test[,7], EW_MSE_on_test[,7]))
  
  colnames(analysis_data_corr0.7) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.8 <- as.data.frame(cbind(rep(scenario$corr[8],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,8], train_av_var[,8], 
                                               test_mean_var[,8], test_min_var[,8], test_max_var[,8],
                                               train_mean_var[,8], train_min_var[,8], train_max_var[,8],
                                               test_mean_corr[,8], test_median_corr[,8], test_min_corr[,8], test_max_corr[,8],
                                               test_low_quantile_corr[,8], test_high_quantile_corr[,8], 
                                               train_mean_corr[,8], train_median_corr[,8], train_min_corr[,8], train_max_corr[,8],
                                               train_low_quantile_corr[,8], train_high_quantile_corr[,8], cv_lambda_agg[,,8], MSE_cv_on_test_agg[,,8],
                                               which_min_MSE_on_test[,8], min_MSE_on_test[,8],
                                               OW_MSE_on_test[,8], EW_MSE_on_test[,8]))
  
  colnames(analysis_data_corr0.8) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.9 <- as.data.frame(cbind(rep(scenario$corr[9],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,9], train_av_var[,9],
                                               test_mean_var[,9], test_min_var[,9], test_max_var[,9],
                                               train_mean_var[,9], train_min_var[,9], train_max_var[,9],
                                               test_mean_corr[,9], test_median_corr[,9], test_min_corr[,9], test_max_corr[,9],
                                               test_low_quantile_corr[,9], test_high_quantile_corr[,9], 
                                               train_mean_corr[,9], train_median_corr[,9], train_min_corr[,9], train_max_corr[,9],
                                               train_low_quantile_corr[,9], train_high_quantile_corr[,9], cv_lambda_agg[,,9], MSE_cv_on_test_agg[,,9],
                                               which_min_MSE_on_test[,9], min_MSE_on_test[,9],
                                               OW_MSE_on_test[,9], EW_MSE_on_test[,9]))
  
  colnames(analysis_data_corr0.9) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                        "test_mean_var", "test_min_var", "test_max_var",
                                        "train_mean_var", "train_min_var", "train_max_var",
                                        "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                        "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                        "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data <- rbind(analysis_data_corr0.1,analysis_data_corr0.2,
                         analysis_data_corr0.3,analysis_data_corr0.4,analysis_data_corr0.5,
                         analysis_data_corr0.6,analysis_data_corr0.7,analysis_data_corr0.8,
                         analysis_data_corr0.9)
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  analysis_data <- data.frame(rep("linVar1to1.2", nrow(analysis_data)), rep(scenario$k, nrow(analysis_data)), analysis_data)
  
  colnames(analysis_data)[1] <- c("VarType")
  colnames(analysis_data)[2] <- c("NumFC")
  
  
  
  if(running == 1){
    analysis_datapart <- analysis_data}else{
      analysis_datapart <- rbind(analysis_datapart, analysis_data)
    }
  
  
  print(c("finished running =", running))
}


# table(duplicated(analysis_datapart))

# save(analysis_datapart, file ="linearVariance1-2_5FC.RData")


################################################################################

### 2. Procedure (Pairwise Constant Correlations; Special Case)

# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (12,15) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.2,
# 1.5, 2, 4, 9)
# -> Resulting in 20 Separate Runs


scenario$k <- 15  # 12, 15

scenario$train <- c(20,20) # for 12, 15 FC

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # Pairwise Correlation


set.seed(11111111)

# 250 Repetitions of Each Parameter Constellation (Scenario) in Total 
# (Due to scenario$train <- c(20,20))

for(running in 1:125){
  
  # Data Containers
  
  no_folds <- rep(NA, 4)
  
  test_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.2, 1.5, 2, 4, 9)
    
    Sigma <- fun_MakeSigmaLin(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    #Sigma <- fun_MakeSigmaExp(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character)
      
      no_folds <- c(5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      no_folds_chr <- c(5, 10, "LTO", "LOO")
      
      
      # Data Containers
      
      weights_on_test <- list()
      MSE_on_test <- list()
      
      
      # Generate New Error Data in Each Run
      
      scenario$test <- scenario$populationsize - scenario$train[ntrain]
      
      ErrorData  <- fun_MakeErrors(Sigma, scenario$populationsize, scenario$train[ntrain], scenario$test, runif(1, 1, 100000000))
      
      
      # Estimation of Correlation-Related Variables on Test Data
      
      true_cormat <- cor(ErrorData$Test)
      
      diag(true_cormat) <- NA 
      
      true_mean_corr <- mean(true_cormat, na.rm = TRUE)
      
      true_median_corr <- median(true_cormat, na.rm = TRUE)
      
      true_min_corr <- min(true_cormat, na.rm = TRUE)
      
      true_max_corr <- max(true_cormat, na.rm = TRUE)
      
      true_low_quantile_corr <- quantile(true_cormat, probs = 0.2, na.rm = TRUE)
      
      true_high_quantile_corr <- quantile(true_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- quantile(train_cormat, probs = 0.2, na.rm = TRUE)
      
      high_quantile_corr_est <- quantile(train_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Variance-Related Variables on Test Data
      
      true_Sigma <- 1/scenario$test * t(ErrorData$Test) %*% ErrorData$Test
      
      true_av_var <- (max(diag(true_Sigma))-min(diag(true_Sigma)))/(scenario$k-1)
      
      true_mean_var <- mean(diag(true_Sigma))
      
      true_min_var <- min(diag(true_Sigma))
      
      true_max_var <- max(diag(true_Sigma))
      
      
      # Estimation of Variance-Related Variables on Training Data
      
      Sigma_est <- 1/scenario$train[ntrain] * t(ErrorData$Training) %*% ErrorData$Training
      
      av_var_est <- (max(diag(Sigma_est))-min(diag(Sigma_est)))/(scenario$k-1)
      
      mean_var_est <- mean(diag(Sigma_est))
      
      min_var_est <- min(diag(Sigma_est))
      
      max_var_est <- max(diag(Sigma_est))
      
      
      # Calculate OW on Full Training Set and Determine Shrinkage Path to EW 
      
      OW_emp <- fun_ComputeOW(ErrorData$Training, scenario$source_vector) 
      
      weights_on_test <- fun_ComputeTTSWeights(OW_emp, target_vector, granularity)  
      
      # Calculate MSE Values for Shrinkage Path on Test Set
      
      MSE_on_test <- fun_DetermineMSE(this.weights = weights_on_test, this.data = ErrorData$Test)  
      
      
      # Data Containers
      
      cv_val_result_per_fold_agg <- list()
      
      
      for(folds in 1:length(no_folds)){
        
        this.no_folds <- no_folds[folds]
        
        
        # Split Training Data in Folds for CV
        
        flds <- cvFolds(nrow(ErrorData$Training), K = this.no_folds)
        
        # Data Containers
        
        cv_weight.Candidates <- list()
        cv_val_result_per_fold <- list()
        
        
        for (this.fold in 1: this.no_folds){
          
          # Split Folds in Validation and Train (Calibration) Set for CV
          
          cv_this.val   <- ErrorData$Training[flds$subsets[flds$which == this.fold], ]   
          cv_this.train <- ErrorData$Training[flds$subsets[flds$which != this.fold], ]
          
          
          # Calculate OW on CV-Training Set and Determine Shrinkage Path to EW
          
          cv_this.OW <- fun_ComputeOW(this.data = cv_this.train, this.source_vector = source_vector) 
          
          cv_weight.Candidates[[this.fold]] <- fun_ComputeTTSWeights(cv_this.OW , target_vector, granularity)
          
          # Calculate MSE Values for Shrinkage Path on Respective Validation Set
          
          cv_val_result_per_fold[[this.fold]] <- fun_DetermineMSE(this.weights = cv_weight.Candidates[[this.fold]], this.data = cv_this.val) 
          
        }
        
        # Average MSE Values over Folds
        
        cv_val_result_per_fold_agg[[folds]] <- Reduce("+", cv_val_result_per_fold) / length(cv_val_result_per_fold) 
        
        # Identify CV-optimal Shrinkage Level 
        
        cv_lambda_agg[ntrain,folds,corr] <- which.min(cv_val_result_per_fold_agg[[folds]])-1
        
        # Identify Resulting MSE Value on Test Data with CV-optimal Shrinkage Level
        
        MSE_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[cv_lambda_agg[ntrain,folds,corr]+1]
        
      }
      
      
      test_av_var[ntrain,corr]<- true_av_var
      
      train_av_var[ntrain,corr] <- av_var_est
      
      
      test_mean_var[ntrain,corr] <-  true_mean_var
      test_min_var[ntrain,corr] <- true_min_var
      test_max_var[ntrain,corr] <-  true_max_var
      
      train_mean_var[ntrain,corr] <- mean_var_est
      train_min_var[ntrain,corr] <- min_var_est
      train_max_var[ntrain,corr] <- max_var_est
      
      test_mean_corr[ntrain,corr] <- true_mean_corr
      test_median_corr[ntrain,corr] <- true_median_corr
      test_min_corr[ntrain,corr] <- true_min_corr
      test_max_corr[ntrain,corr] <- true_max_corr
      test_low_quantile_corr[ntrain,corr] <- true_low_quantile_corr
      test_high_quantile_corr[ntrain,corr] <- true_high_quantile_corr
      
      train_mean_corr[ntrain,corr] <- mean_corr_est
      train_median_corr[ntrain,corr] <- median_corr_est
      train_min_corr[ntrain,corr] <- min_corr_est
      train_max_corr[ntrain,corr] <- max_corr_est
      train_low_quantile_corr[ntrain,corr] <- low_quantile_corr_est
      train_high_quantile_corr[ntrain,corr] <- high_quantile_corr_est
      
      
      # Identify Truly Optimal Shrinkage Level and Resulting MSE on Test Data
      
      which_min_MSE_on_test[ntrain,corr] <- which.min(MSE_on_test)-1
      
      min_MSE_on_test[ntrain,corr] <- min(MSE_on_test)
      
      
      # MSE for OW on Test Data
      
      OW_MSE_on_test[ntrain,corr] <- MSE_on_test[1]
      
      # MSE for EW on Test Data
      
      EW_MSE_on_test[ntrain,corr] <- MSE_on_test[101]
      
      
      print(c("finished ntrain =", ntrain))
    }
    print(c("finished corr =", corr))
  }
  
  
  # Organize and Store Results 
  
  analysis_data_corr0.1 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                               test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                               train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                               test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                               test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                               train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                               train_low_quantile_corr[,1], train_high_quantile_corr[,1], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,1],  rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,1],
                                               which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                               OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(analysis_data_corr0.1) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.2 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                               test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                               train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                               test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                               test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                               train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                               train_low_quantile_corr[,2], train_high_quantile_corr[,2],  rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,2], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,2],
                                               which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                               OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(analysis_data_corr0.2) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.3 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                               test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                               train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                               test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                               test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                               train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                               train_low_quantile_corr[,3], train_high_quantile_corr[,3], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,3], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,3],
                                               which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                               OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(analysis_data_corr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.4 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                               test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                               train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                               test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                               test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                               train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                               train_low_quantile_corr[,4], train_high_quantile_corr[,4], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,4], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,4],
                                               which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                               OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(analysis_data_corr0.4) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                         "test_mean_var", "test_min_var", "test_max_var",
                                         "train_mean_var", "train_min_var", "train_max_var",
                                         "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                         "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                         "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                         "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.5 <- as.data.frame(cbind(rep(scenario$corr[5],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,5], train_av_var[,5], 
                                               test_mean_var[,5], test_min_var[,5], test_max_var[,5],
                                               train_mean_var[,5], train_min_var[,5], train_max_var[,5],
                                               test_mean_corr[,5], test_median_corr[,5], test_min_corr[,5], test_max_corr[,5],
                                               test_low_quantile_corr[,5], test_high_quantile_corr[,5], 
                                               train_mean_corr[,5], train_median_corr[,5], train_min_corr[,5], train_max_corr[,5],
                                               train_low_quantile_corr[,5], train_high_quantile_corr[,5], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,5], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,5],
                                               which_min_MSE_on_test[,5], min_MSE_on_test[,5],
                                               OW_MSE_on_test[,5], EW_MSE_on_test[,5]))
  
  colnames(analysis_data_corr0.5) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.6 <- as.data.frame(cbind(rep(scenario$corr[6],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,6], train_av_var[,6],
                                               test_mean_var[,6], test_min_var[,6], test_max_var[,6],
                                               train_mean_var[,6], train_min_var[,6], train_max_var[,6],
                                               test_mean_corr[,6], test_median_corr[,6], test_min_corr[,6], test_max_corr[,6],
                                               test_low_quantile_corr[,6], test_high_quantile_corr[,6], 
                                               train_mean_corr[,6], train_median_corr[,6], train_min_corr[,6], train_max_corr[,6],
                                               train_low_quantile_corr[,6], train_high_quantile_corr[,6], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,6], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,6],
                                               which_min_MSE_on_test[,6], min_MSE_on_test[,6],
                                               OW_MSE_on_test[,6], EW_MSE_on_test[,6]))
  
  colnames(analysis_data_corr0.6) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                        "test_mean_var", "test_min_var", "test_max_var",
                                        "train_mean_var", "train_min_var", "train_max_var",
                                        "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                        "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                        "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.7 <- as.data.frame(cbind(rep(scenario$corr[7],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,7], train_av_var[,7], 
                                               test_mean_var[,7], test_min_var[,7], test_max_var[,7],
                                               train_mean_var[,7], train_min_var[,7], train_max_var[,7],
                                               test_mean_corr[,7], test_median_corr[,7], test_min_corr[,7], test_max_corr[,7],
                                               test_low_quantile_corr[,7], test_high_quantile_corr[,7], 
                                               train_mean_corr[,7], train_median_corr[,7], train_min_corr[,7], train_max_corr[,7],
                                               train_low_quantile_corr[,7], train_high_quantile_corr[,7], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,7], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,7],
                                               which_min_MSE_on_test[,7], min_MSE_on_test[,7],
                                               OW_MSE_on_test[,7], EW_MSE_on_test[,7]))
  
  colnames(analysis_data_corr0.7) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.8 <- as.data.frame(cbind(rep(scenario$corr[8],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,8], train_av_var[,8], 
                                               test_mean_var[,8], test_min_var[,8], test_max_var[,8],
                                               train_mean_var[,8], train_min_var[,8], train_max_var[,8],
                                               test_mean_corr[,8], test_median_corr[,8], test_min_corr[,8], test_max_corr[,8],
                                               test_low_quantile_corr[,8], test_high_quantile_corr[,8], 
                                               train_mean_corr[,8], train_median_corr[,8], train_min_corr[,8], train_max_corr[,8],
                                               train_low_quantile_corr[,8], train_high_quantile_corr[,8], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,8], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,8],
                                               which_min_MSE_on_test[,8], min_MSE_on_test[,8],
                                               OW_MSE_on_test[,8], EW_MSE_on_test[,8]))
  
  colnames(analysis_data_corr0.8) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_corr0.9 <- as.data.frame(cbind(rep(scenario$corr[9],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,9], train_av_var[,9],
                                               test_mean_var[,9], test_min_var[,9], test_max_var[,9],
                                               train_mean_var[,9], train_min_var[,9], train_max_var[,9],
                                               test_mean_corr[,9], test_median_corr[,9], test_min_corr[,9], test_max_corr[,9],
                                               test_low_quantile_corr[,9], test_high_quantile_corr[,9], 
                                               train_mean_corr[,9], train_median_corr[,9], train_min_corr[,9], train_max_corr[,9],
                                               train_low_quantile_corr[,9], train_high_quantile_corr[,9], rep(NA,length(scenario$train)), 
                                               cv_lambda_agg[,,9], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,9],
                                               which_min_MSE_on_test[,9], min_MSE_on_test[,9],
                                               OW_MSE_on_test[,9], EW_MSE_on_test[,9]))
  
  colnames(analysis_data_corr0.9) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                        "test_mean_var", "test_min_var", "test_max_var",
                                        "train_mean_var", "train_min_var", "train_max_var",
                                        "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                        "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                        "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                        "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data <- rbind(analysis_data_corr0.1,analysis_data_corr0.2,
                         analysis_data_corr0.3,analysis_data_corr0.4,analysis_data_corr0.5,
                         analysis_data_corr0.6,analysis_data_corr0.7,analysis_data_corr0.8,
                         analysis_data_corr0.9)
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  analysis_data <- data.frame(rep("linVar1to1.2", nrow(analysis_data)), rep(scenario$k, nrow(analysis_data)), analysis_data)
  
  colnames(analysis_data)[1] <- c("VarType")
  colnames(analysis_data)[2] <- c("NumFC")
  
  
  if(running == 1){
    analysis_datapart <- analysis_data}else{
      analysis_datapart <- rbind(analysis_datapart, analysis_data)
    }
  
  print(c("finished running =", running))
}


# table(duplicated(analysis_datapart))

# save(analysis_datapart, file ="linearVariance1-2_15FCn20.RData")


################################################################################

### 3. Procedure (Differing Correlations; General Scenarios)


# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (5,8,10,12,15) and Choosing Respective Number
# of Training Observations (Ensuring Sufficient Amount of
# Training Observations to Calculate OW) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.2,
# 1.5, 2, 4, 9)
# -> Resulting in 50 Separate Runs


scenario$k <- 5 # 5, 8, 10, 12, 15

scenario$train <- c(10,20,30,40,50,60,70,80,90,100,125,150,175,200) # for 5 FC
# scenario$train <- c(20,30,40,50,60,70,80,90,100,125,150,175,200) # for 8,10 FC
# scenario$train <- c(30,40,50,60,70,80,90,100,125,150,175,200) # for 12,15 FC

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8) # Pairwise Correlation


set.seed(11111111)

# 250 Repetitions of Each Parameter Constellation (Scenario)

for(running in 1:250){
  
  # Data Containers
  
  no_folds <- rep(NA, 5)
  
  test_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.2, 1.5, 2, 4, 9)
    
    Sigma <- fun_MakeSigmaLinDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    #Sigma <- fun_MakeSigmaExpDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character)
      
      no_folds <- c(2, 5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      no_folds_chr <- c(2, 5, 10, "LTO", "LOO")
      
      
      # Data Containers
      
      weights_on_test <- list()
      MSE_on_test <- list()
      
      
      # Generate New Error Data in Each Run
      
      scenario$test <- scenario$populationsize - scenario$train[ntrain]
      
      ErrorData  <- fun_MakeErrors(Sigma, scenario$populationsize, scenario$train[ntrain], scenario$test, runif(1, 1, 100000000))
      
      
      # Estimation of Correlation-Related Variables on Test Data
      
      true_cormat <- cor(ErrorData$Test)
      
      diag(true_cormat) <- NA 
      
      true_mean_corr <- mean(true_cormat, na.rm = TRUE)
      
      true_median_corr <- median(true_cormat, na.rm = TRUE)
      
      true_min_corr <- min(true_cormat, na.rm = TRUE)
      
      true_max_corr <- max(true_cormat, na.rm = TRUE)
      
      true_low_quantile_corr <- quantile(true_cormat, probs = 0.2, na.rm = TRUE)
      
      true_high_quantile_corr <- quantile(true_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- quantile(train_cormat, probs = 0.2, na.rm = TRUE)
      
      high_quantile_corr_est <- quantile(train_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Variance-Related Variables on Test Data
      
      true_Sigma <- 1/scenario$test * t(ErrorData$Test) %*% ErrorData$Test
      
      true_av_var <- (max(diag(true_Sigma))-min(diag(true_Sigma)))/(scenario$k-1)
      
      true_mean_var <- mean(diag(true_Sigma))
      
      true_min_var <- min(diag(true_Sigma))
      
      true_max_var <- max(diag(true_Sigma))
      
      
      # Estimation of Variance-Related Variables on Training Data
      
      Sigma_est <- 1/scenario$train[ntrain] * t(ErrorData$Training) %*% ErrorData$Training
      
      av_var_est <- (max(diag(Sigma_est))-min(diag(Sigma_est)))/(scenario$k-1)
      
      mean_var_est <- mean(diag(Sigma_est))
      
      min_var_est <- min(diag(Sigma_est))
      
      max_var_est <- max(diag(Sigma_est))
      
      
      # Calculate OW on Full Training Set and Determine Shrinkage Path to EW 
      
      OW_emp <- fun_ComputeOW(ErrorData$Training, scenario$source_vector) 
      
      weights_on_test <- fun_ComputeTTSWeights(OW_emp, target_vector, granularity)  
      
      # Calculate MSE Values for Shrinkage Path on Test Set
      
      MSE_on_test <- fun_DetermineMSE(this.weights = weights_on_test, this.data = ErrorData$Test)  
      
      
      # Data Containers
      
      cv_val_result_per_fold_agg <- list()
      
      
      for(folds in 1:length(no_folds)){
        
        this.no_folds <- no_folds[folds]
        
        
        # Split Training Data in Folds for CV
        
        flds <- cvFolds(nrow(ErrorData$Training), K = this.no_folds)
        
        # Data Containers
        
        cv_weight.Candidates <- list()
        cv_val_result_per_fold <- list()
        
        
        for (this.fold in 1: this.no_folds){
          
          # Split Folds in Validation and Train (Calibration) Set for CV
          
          cv_this.val   <- ErrorData$Training[flds$subsets[flds$which == this.fold], ]   
          cv_this.train <- ErrorData$Training[flds$subsets[flds$which != this.fold], ]
          
          
          # Calculate OW on CV-Training Set and Determine Shrinkage Path to EW
          
          cv_this.OW <- fun_ComputeOW(this.data = cv_this.train, this.source_vector = source_vector) 
          
          cv_weight.Candidates[[this.fold]] <- fun_ComputeTTSWeights(cv_this.OW , target_vector, granularity)
          
          # Calculate MSE Values for Shrinkage Path on Respective Validation Set
          
          cv_val_result_per_fold[[this.fold]] <- fun_DetermineMSE(this.weights = cv_weight.Candidates[[this.fold]], this.data = cv_this.val) 
          
        }
        
        # Average MSE Values over Folds
        
        cv_val_result_per_fold_agg[[folds]] <- Reduce("+", cv_val_result_per_fold) / length(cv_val_result_per_fold) 
        
        # Identify CV-optimal Shrinkage Level 
        
        cv_lambda_agg[ntrain,folds,corr] <- which.min(cv_val_result_per_fold_agg[[folds]])-1
        
        # Identify Resulting MSE Value on Test Data with CV-optimal Shrinkage Level
        
        MSE_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[cv_lambda_agg[ntrain,folds,corr]+1]
        
      }
      
      
      test_av_var[ntrain,corr]<- true_av_var
      
      train_av_var[ntrain,corr] <- av_var_est
      
      test_mean_var[ntrain,corr] <-  true_mean_var
      test_min_var[ntrain,corr] <- true_min_var
      test_max_var[ntrain,corr] <-  true_max_var
      
      train_mean_var[ntrain,corr] <- mean_var_est
      train_min_var[ntrain,corr] <- min_var_est
      train_max_var[ntrain,corr] <- max_var_est
      
      test_mean_corr[ntrain,corr] <- true_mean_corr
      test_median_corr[ntrain,corr] <- true_median_corr
      test_min_corr[ntrain,corr] <- true_min_corr
      test_max_corr[ntrain,corr] <- true_max_corr
      test_low_quantile_corr[ntrain,corr] <- true_low_quantile_corr
      test_high_quantile_corr[ntrain,corr] <- true_high_quantile_corr
      
      train_mean_corr[ntrain,corr] <- mean_corr_est
      train_median_corr[ntrain,corr] <- median_corr_est
      train_min_corr[ntrain,corr] <- min_corr_est
      train_max_corr[ntrain,corr] <- max_corr_est
      train_low_quantile_corr[ntrain,corr] <- low_quantile_corr_est
      train_high_quantile_corr[ntrain,corr] <- high_quantile_corr_est
      
      
      # Identify Truly Optimal Shrinkage Level and Resulting MSE on Test Data
      
      which_min_MSE_on_test[ntrain,corr] <- which.min(MSE_on_test)-1
      
      min_MSE_on_test[ntrain,corr] <- min(MSE_on_test)
      
      
      # MSE for OW on Test Data
      
      OW_MSE_on_test[ntrain,corr] <- MSE_on_test[1]
      
      # MSE for EW on Test Data
      
      EW_MSE_on_test[ntrain,corr] <- MSE_on_test[101]
      
      
      print(c("finished ntrain =", ntrain))
    }
    print(c("finished corr =", corr))
  }
  
  
  # Organize and Store Results 
  
  analysis_data_devcorr0.2 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                                  test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                                  train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                                  test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                                  test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                                  train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                                  train_low_quantile_corr[,1], train_high_quantile_corr[,1], cv_lambda_agg[,,1], MSE_cv_on_test_agg[,,1],
                                                  which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                                  OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(analysis_data_devcorr0.2) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.3 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                                  test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                                  train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                                  test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                                  test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                                  train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                                  train_low_quantile_corr[,2], train_high_quantile_corr[,2], cv_lambda_agg[,,2], MSE_cv_on_test_agg[,,2],
                                                  which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                                  OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(analysis_data_devcorr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.4 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                                  test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                                  train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                                  test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                                  test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                                  train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                                  train_low_quantile_corr[,3], train_high_quantile_corr[,3], cv_lambda_agg[,,3], MSE_cv_on_test_agg[,,3],
                                                  which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                                  OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(analysis_data_devcorr0.4) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.5 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                                  test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                                  train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                                  test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                                  test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                                  train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                                  train_low_quantile_corr[,4], train_high_quantile_corr[,4], cv_lambda_agg[,,4], MSE_cv_on_test_agg[,,4],
                                                  which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                                  OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(analysis_data_devcorr0.5) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.6 <- as.data.frame(cbind(rep(scenario$corr[5],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,5], train_av_var[,5], 
                                                  test_mean_var[,5], test_min_var[,5], test_max_var[,5],
                                                  train_mean_var[,5], train_min_var[,5], train_max_var[,5],
                                                  test_mean_corr[,5], test_median_corr[,5], test_min_corr[,5], test_max_corr[,5],
                                                  test_low_quantile_corr[,5], test_high_quantile_corr[,5], 
                                                  train_mean_corr[,5], train_median_corr[,5], train_min_corr[,5], train_max_corr[,5],
                                                  train_low_quantile_corr[,5], train_high_quantile_corr[,5], cv_lambda_agg[,,5], MSE_cv_on_test_agg[,,5],
                                                  which_min_MSE_on_test[,5], min_MSE_on_test[,5],
                                                  OW_MSE_on_test[,5], EW_MSE_on_test[,5]))
  
  colnames(analysis_data_devcorr0.6) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                             "test_mean_var", "test_min_var", "test_max_var",
                                             "train_mean_var", "train_min_var", "train_max_var",
                                             "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                             "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                             "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.7 <- as.data.frame(cbind(rep(scenario$corr[6],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,6], train_av_var[,6],
                                                  test_mean_var[,6], test_min_var[,6], test_max_var[,6],
                                                  train_mean_var[,6], train_min_var[,6], train_max_var[,6],
                                                  test_mean_corr[,6], test_median_corr[,6], test_min_corr[,6], test_max_corr[,6],
                                                  test_low_quantile_corr[,6], test_high_quantile_corr[,6], 
                                                  train_mean_corr[,6], train_median_corr[,6], train_min_corr[,6], train_max_corr[,6],
                                                  train_low_quantile_corr[,6], train_high_quantile_corr[,6], cv_lambda_agg[,,6], MSE_cv_on_test_agg[,,6],
                                                  which_min_MSE_on_test[,6], min_MSE_on_test[,6],
                                                  OW_MSE_on_test[,6], EW_MSE_on_test[,6]))
  
  colnames(analysis_data_devcorr0.7) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                           "test_mean_var", "test_min_var", "test_max_var",
                                           "train_mean_var", "train_min_var", "train_max_var",
                                           "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                           "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                           "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.8 <- as.data.frame(cbind(rep(scenario$corr[7],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,7], train_av_var[,7], 
                                                  test_mean_var[,7], test_min_var[,7], test_max_var[,7],
                                                  train_mean_var[,7], train_min_var[,7], train_max_var[,7],
                                                  test_mean_corr[,7], test_median_corr[,7], test_min_corr[,7], test_max_corr[,7],
                                                  test_low_quantile_corr[,7], test_high_quantile_corr[,7], 
                                                  train_mean_corr[,7], train_median_corr[,7], train_min_corr[,7], train_max_corr[,7],
                                                  train_low_quantile_corr[,7], train_high_quantile_corr[,7], cv_lambda_agg[,,7], MSE_cv_on_test_agg[,,7],
                                                  which_min_MSE_on_test[,7], min_MSE_on_test[,7],
                                                  OW_MSE_on_test[,7], EW_MSE_on_test[,7]))
  
  colnames(analysis_data_devcorr0.8) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                             "test_mean_var", "test_min_var", "test_max_var",
                                             "train_mean_var", "train_min_var", "train_max_var",
                                             "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                             "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                             "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data <- rbind(analysis_data_devcorr0.2,analysis_data_devcorr0.3,analysis_data_devcorr0.4,analysis_data_devcorr0.5,
                         analysis_data_devcorr0.6,analysis_data_devcorr0.7,analysis_data_devcorr0.8)
  
  analysis_data$gen_corr <- paste("Dev", analysis_data$gen_corr, sep = "")
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  analysis_data <- data.frame(rep("linVar1to1.2", nrow(analysis_data)), rep(scenario$k, nrow(analysis_data)), analysis_data)
  
  colnames(analysis_data)[1] <- c("VarType")
  colnames(analysis_data)[2] <- c("NumFC")
  
  
  if(running == 1){
    analysis_datapart <- analysis_data}else{
      analysis_datapart <- rbind(analysis_datapart, analysis_data)
    }
  
  
  print(c("finished running =", running))
}

# table(duplicated(analysis_datapart))

# save(analysis_datapart, file ="linearVariance1-2_5FCDevCorr.RData")

################################################################################

### 4. Procedure (Differing Correlations; Special Case)

# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (12,15) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.2,
# 1.5, 2, 4, 9)
# -> Resulting in 20 Separate Runs

scenario$k <- 15  # 12, 15

scenario$train <- c(20,20) # for 12, 15 FC

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8) # Pairwise Correlation


set.seed(11111111)

# 250 Repetitions of Each Parameter Constellation (Scenario) in Total 
# (Due to scenario$train <- c(20,20))

for(running in 1:125){
  
  # Data Containers
  
  no_folds <- rep(NA, 4)
  
  test_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_av_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_var <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  test_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  test_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  train_mean_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_median_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_min_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_max_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_low_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  train_high_quantile_corr <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.2, 1.5, 2, 4, 9)
    
    Sigma <- fun_MakeSigmaLinDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    #Sigma <- fun_MakeSigmaExpDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.2)
    
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character)
      
      no_folds <- c(5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      no_folds_chr <- c(5, 10, "LTO", "LOO")
      
      
      # Data Containers
      
      weights_on_test <- list()
      MSE_on_test <- list()
      
      
      # Generate New Error Data in Each Run
      
      scenario$test <- scenario$populationsize - scenario$train[ntrain]
      
      ErrorData  <- fun_MakeErrors(Sigma, scenario$populationsize, scenario$train[ntrain], scenario$test, runif(1, 1, 100000000))
      
      
      # Estimation of Correlation-Related Variables on Test Data
      
      true_cormat <- cor(ErrorData$Test)
      
      diag(true_cormat) <- NA 
      
      true_mean_corr <- mean(true_cormat, na.rm = TRUE)
      
      true_median_corr <- median(true_cormat, na.rm = TRUE)
      
      true_min_corr <- min(true_cormat, na.rm = TRUE)
      
      true_max_corr <- max(true_cormat, na.rm = TRUE)
      
      true_low_quantile_corr <- quantile(true_cormat, probs = 0.2, na.rm = TRUE)
      
      true_high_quantile_corr <- quantile(true_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- quantile(train_cormat, probs = 0.2, na.rm = TRUE)
      
      high_quantile_corr_est <- quantile(train_cormat, probs = 0.8, na.rm = TRUE)
      
      
      # Estimation of Variance-Related Variables on Test Data
      
      true_Sigma <- 1/scenario$test * t(ErrorData$Test) %*% ErrorData$Test
      
      true_av_var <- (max(diag(true_Sigma))-min(diag(true_Sigma)))/(scenario$k-1)
      
      true_mean_var <- mean(diag(true_Sigma))
      
      true_min_var <- min(diag(true_Sigma))
      
      true_max_var <- max(diag(true_Sigma))
      
      
      # Estimation of Variance-Related Variables on Training Data
      
      Sigma_est <- 1/scenario$train[ntrain] * t(ErrorData$Training) %*% ErrorData$Training
      
      av_var_est <- (max(diag(Sigma_est))-min(diag(Sigma_est)))/(scenario$k-1)
      
      mean_var_est <- mean(diag(Sigma_est))
      
      min_var_est <- min(diag(Sigma_est))
      
      max_var_est <- max(diag(Sigma_est))
      
      
      # Calculate OW on Full Training Set and Determine Shrinkage Path to EW 
      
      OW_emp <- fun_ComputeOW(ErrorData$Training, scenario$source_vector)
      
      weights_on_test <- fun_ComputeTTSWeights(OW_emp, target_vector, granularity)  
      
      # Calculate MSE Values for Shrinkage Path on Test Set
      
      MSE_on_test <- fun_DetermineMSE(this.weights = weights_on_test, this.data = ErrorData$Test)  
      
      
      # Data Containers
      
      cv_val_result_per_fold_agg <- list()
      
      
      for(folds in 1:length(no_folds)){
        
        this.no_folds <- no_folds[folds]
        
        
        # Split Training Data in Folds for CV
        
        flds <- cvFolds(nrow(ErrorData$Training), K = this.no_folds)
        
        # Data Containers
        
        cv_weight.Candidates <- list()
        cv_val_result_per_fold <- list()
        
        
        for (this.fold in 1: this.no_folds){
          
          # Split Folds in Validation and Train (Calibration) Set for CV
          
          cv_this.val   <- ErrorData$Training[flds$subsets[flds$which == this.fold], ]   
          cv_this.train <- ErrorData$Training[flds$subsets[flds$which != this.fold], ]
          
          
          # Calculate OW on CV-Training Set and Determine Shrinkage Path to EW
          
          cv_this.OW <- fun_ComputeOW(this.data = cv_this.train, this.source_vector = source_vector) 
          
          cv_weight.Candidates[[this.fold]] <- fun_ComputeTTSWeights(cv_this.OW , target_vector, granularity)
          
          # Calculate MSE Values for Shrinkage Path on Respective Validation Set
          
          cv_val_result_per_fold[[this.fold]] <- fun_DetermineMSE(this.weights = cv_weight.Candidates[[this.fold]], this.data = cv_this.val) 
          
        }
        
        # Average MSE Values over Folds
        
        cv_val_result_per_fold_agg[[folds]] <- Reduce("+", cv_val_result_per_fold) / length(cv_val_result_per_fold) 
        
        # Identify CV-optimal Shrinkage Level 
        
        cv_lambda_agg[ntrain,folds,corr] <- which.min(cv_val_result_per_fold_agg[[folds]])-1
        
        # Identify Resulting MSE Value on Test Data with CV-optimal Shrinkage Level
        
        MSE_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[cv_lambda_agg[ntrain,folds,corr]+1]
        
      }
      
      
      test_av_var[ntrain,corr]<- true_av_var
      
      train_av_var[ntrain,corr] <- av_var_est
      
      test_mean_var[ntrain,corr] <-  true_mean_var
      test_min_var[ntrain,corr] <- true_min_var
      test_max_var[ntrain,corr] <-  true_max_var
      
      train_mean_var[ntrain,corr] <- mean_var_est
      train_min_var[ntrain,corr] <- min_var_est
      train_max_var[ntrain,corr] <- max_var_est
      
      test_mean_corr[ntrain,corr] <- true_mean_corr
      test_median_corr[ntrain,corr] <- true_median_corr
      test_min_corr[ntrain,corr] <- true_min_corr
      test_max_corr[ntrain,corr] <- true_max_corr
      test_low_quantile_corr[ntrain,corr] <- true_low_quantile_corr
      test_high_quantile_corr[ntrain,corr] <- true_high_quantile_corr
      
      train_mean_corr[ntrain,corr] <- mean_corr_est
      train_median_corr[ntrain,corr] <- median_corr_est
      train_min_corr[ntrain,corr] <- min_corr_est
      train_max_corr[ntrain,corr] <- max_corr_est
      train_low_quantile_corr[ntrain,corr] <- low_quantile_corr_est
      train_high_quantile_corr[ntrain,corr] <- high_quantile_corr_est
      
      
      # Identify Truly Optimal Shrinkage Level and Resulting MSE on Test Data
      
      which_min_MSE_on_test[ntrain,corr] <- which.min(MSE_on_test)-1
      
      min_MSE_on_test[ntrain,corr] <- min(MSE_on_test)
      
      
      # MSE for OW on Test Data
      
      OW_MSE_on_test[ntrain,corr] <- MSE_on_test[1]
      
      # MSE for EW on Test Data
      
      EW_MSE_on_test[ntrain,corr] <- MSE_on_test[101]
      
      
      print(c("finished ntrain =", ntrain))
    }
    print(c("finished corr =", corr))
  }
  
  
  # Organize and Store Results 
  
  analysis_data_devcorr0.2 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                                  test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                                  train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                                  test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                                  test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                                  train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                                  train_low_quantile_corr[,1], train_high_quantile_corr[,1], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,1], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,1],
                                                  which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                                  OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(analysis_data_devcorr0.2) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.3 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                                  test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                                  train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                                  test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                                  test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                                  train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                                  train_low_quantile_corr[,2], train_high_quantile_corr[,2], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,2], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,2],
                                                  which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                                  OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(analysis_data_devcorr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.4 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                                  test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                                  train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                                  test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                                  test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                                  train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                                  train_low_quantile_corr[,3], train_high_quantile_corr[,3], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,3], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,3],
                                                  which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                                  OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(analysis_data_devcorr0.4) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.5 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                                  test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                                  train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                                  test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                                  test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                                  train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                                  train_low_quantile_corr[,4], train_high_quantile_corr[,4], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,4], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,4],
                                                  which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                                  OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(analysis_data_devcorr0.5) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.6 <- as.data.frame(cbind(rep(scenario$corr[5],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,5], train_av_var[,5], 
                                                  test_mean_var[,5], test_min_var[,5], test_max_var[,5],
                                                  train_mean_var[,5], train_min_var[,5], train_max_var[,5],
                                                  test_mean_corr[,5], test_median_corr[,5], test_min_corr[,5], test_max_corr[,5],
                                                  test_low_quantile_corr[,5], test_high_quantile_corr[,5], 
                                                  train_mean_corr[,5], train_median_corr[,5], train_min_corr[,5], train_max_corr[,5],
                                                  train_low_quantile_corr[,5], train_high_quantile_corr[,5], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,5], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,5],
                                                  which_min_MSE_on_test[,5], min_MSE_on_test[,5],
                                                  OW_MSE_on_test[,5], EW_MSE_on_test[,5]))
  
  colnames(analysis_data_devcorr0.6) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                             "test_mean_var", "test_min_var", "test_max_var",
                                             "train_mean_var", "train_min_var", "train_max_var",
                                             "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                             "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                             "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.7 <- as.data.frame(cbind(rep(scenario$corr[6],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,6], train_av_var[,6],
                                                  test_mean_var[,6], test_min_var[,6], test_max_var[,6],
                                                  train_mean_var[,6], train_min_var[,6], train_max_var[,6],
                                                  test_mean_corr[,6], test_median_corr[,6], test_min_corr[,6], test_max_corr[,6],
                                                  test_low_quantile_corr[,6], test_high_quantile_corr[,6], 
                                                  train_mean_corr[,6], train_median_corr[,6], train_min_corr[,6], train_max_corr[,6],
                                                  train_low_quantile_corr[,6], train_high_quantile_corr[,6], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,6], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,6],
                                                  which_min_MSE_on_test[,6], min_MSE_on_test[,6],
                                                  OW_MSE_on_test[,6], EW_MSE_on_test[,6]))
  
  colnames(analysis_data_devcorr0.7) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                           "test_mean_var", "test_min_var", "test_max_var",
                                           "train_mean_var", "train_min_var", "train_max_var",
                                           "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                           "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                           "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data_devcorr0.8 <- as.data.frame(cbind(rep(scenario$corr[7],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,7], train_av_var[,7], 
                                                  test_mean_var[,7], test_min_var[,7], test_max_var[,7],
                                                  train_mean_var[,7], train_min_var[,7], train_max_var[,7],
                                                  test_mean_corr[,7], test_median_corr[,7], test_min_corr[,7], test_max_corr[,7],
                                                  test_low_quantile_corr[,7], test_high_quantile_corr[,7], 
                                                  train_mean_corr[,7], train_median_corr[,7], train_min_corr[,7], train_max_corr[,7],
                                                  train_low_quantile_corr[,7], train_high_quantile_corr[,7], rep(NA,length(scenario$train)),
                                                  cv_lambda_agg[,,7], rep(NA,length(scenario$train)), MSE_cv_on_test_agg[,,7],
                                                  which_min_MSE_on_test[,7], min_MSE_on_test[,7],
                                                  OW_MSE_on_test[,7], EW_MSE_on_test[,7]))
  
  colnames(analysis_data_devcorr0.8) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                             "test_mean_var", "test_min_var", "test_max_var",
                                             "train_mean_var", "train_min_var", "train_max_var",
                                             "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                             "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                             "train_low_quantile_corr", "train_high_quantile_corr", paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  analysis_data <- rbind(analysis_data_devcorr0.2,analysis_data_devcorr0.3,analysis_data_devcorr0.4,analysis_data_devcorr0.5,
                         analysis_data_devcorr0.6,analysis_data_devcorr0.7,analysis_data_devcorr0.8)
  
  analysis_data$gen_corr <- paste("Dev", analysis_data$gen_corr, sep = "")
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  analysis_data <- data.frame(rep("linVar1to1.2", nrow(analysis_data)), rep(scenario$k, nrow(analysis_data)), analysis_data)
  
  colnames(analysis_data)[1] <- c("VarType")
  colnames(analysis_data)[2] <- c("NumFC")
  
  
  if(running == 1){
    analysis_datapart <- analysis_data}else{
      analysis_datapart <- rbind(analysis_datapart, analysis_data)
    }
  
  
  print(c("finished running =", running))
}


# table(duplicated(analysis_datapart))

# save(analysis_datapart, file ="linearVariance1-1.2_15FCn20DevCorr.RData")
