################################################################################
################################################################################
###                  Evaluation of Shrinkage Correction                      ###
################################################################################
################################################################################

### Generation of Synthetic Datasets for Different Parameter Value 

### Constellations (Scenarios) to Evaluate Shrinkage Bias Correction

################################################################################

source("General Functions and Settings.R")

# Reconstructing the Regression Tree Model (Resulting from "Regression Tree Model.R")
# Needed for the Evaluation

library(rpart)
library(rpart.plot)
library(mlr)

load("Datasets/Regression Dataset/final_regression_complete.RData")

overshrinkTib <- mutate_if(regression_complete, is.character, as.factor)

overshrinkTask <- makeRegrTask(data = overshrinkTib, target = "overshrinkage")

tree <- makeLearner("regr.rpart")

tunedTree <- setHyperPars(tree, par.vals = list(cp = 0.00000028))

tunedTreeModel <- train(tunedTree, overshrinkTask)

treeModelData <- getLearnerModel(tunedTreeModel)

################################################################################

# Two Separate Procedures (Pairwise Constant vs. Differing Error Correlations)

# Do Each of The Procedures for All Combinations of Parameter Values 

### 1. Procedure (Pairwise Constant Correlations)


# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (4,9) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.1,
# 1.25, 1.7, 3, 6, 8)
# -> Resulting in 24 Separate Runs

# Determine the Number of Forecasters (4,9)

scenario$k <- 4 # 4,9

scenario$train <- c(25,50,75,100,150,200) 

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.15,0.3,0.45,0.6,0.75,0.9) # Pairwise Correlation


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
  
  mean_corr_train_group <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  var_range_train_group <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  corr_differ_train <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  corr_factor <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  corr_cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_corr_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.1, 1.25, 1.7, 3, 6, 8)
    
    Sigma <- fun_MakeSigmaLin(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.1)
    #Sigma <- fun_MakeSigmaExp(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.1)
    
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character) and 
      # Number (and Percentage) of Training Observations in CV-Calibration Set
      
      no_folds <- c(2, 5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      n_calib <- c(floor(scenario$train[ntrain]/2), round(0.8*scenario$train[ntrain],0),
                   round(0.9*scenario$train[ntrain],0), (scenario$train[ntrain]-2), (scenario$train[ntrain]-1))
      
      perc_calib <- c(0.5,0.8,0.9,round((scenario$train[ntrain]-2)/scenario$train[ntrain],3), 
                      round((scenario$train[ntrain]-1)/scenario$train[ntrain],3))
      
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
      
      true_low_quantile_corr <- unname(quantile(true_cormat, probs = 0.2, na.rm = TRUE))
      
      true_high_quantile_corr <- unname(quantile(true_cormat, probs = 0.8, na.rm = TRUE))
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- unname(quantile(train_cormat, probs = 0.25, na.rm = TRUE))
      
      high_quantile_corr_est <- unname(quantile(train_cormat, probs = 0.75, na.rm = TRUE))
      
      
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
      
      range_var_est <- (max_var_est - min_var_est)
      
      
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
        
        
        # Rename and Generate Variables for Application of Regression Tree
        
        J <- scenario$k
        n <- scenario$train[ntrain]
        K <- no_folds_chr[folds]
        n_cal <- n_calib[folds]
        perc_cal <- perc_calib[folds]
        
        
        if(mean_corr_est < 0.25){
          mean_corr <- "weak"} else if(0.25 <= mean_corr_est & mean_corr_est < 0.55){
            mean_corr <- "moderate"} else if(0.55 <= mean_corr_est & mean_corr_est < 0.75){
              mean_corr <- "strong"} else{
                mean_corr <- "extreme"  
              } 
        
        
        corr_diff <- (high_quantile_corr_est - low_quantile_corr_est)
        
        
        if(range_var_est < 0.45){
          var_range <- "tiny"} else if(0.45 <= range_var_est & range_var_est < 0.95){
            var_range <- "low"} else if(0.95 <= range_var_est & range_var_est < 2.5){
              var_range <- "medium"} else if(2.5 <= range_var_est & range_var_est < 6){
                var_range <- "high"} else{
                  var_range <- "extreme"  
                } 
        
        
        params <- data.frame(J, n, K, n_cal, perc_cal, mean_corr,
                             corr_diff, var_range)
        
        
        # Application of Regression Tree to Determine Correction Factors
        
        corr_factor[ntrain,folds,corr] <- round(predict(treeModelData, params))
        
        # Apply Shrinkage Correction and Determine Respective MSE Value
        
        corr_cv_lambda_agg[ntrain,folds,corr] <- cv_lambda_agg[ntrain,folds,corr] - corr_factor[ntrain,folds,corr]
        if(corr_cv_lambda_agg[ntrain,folds,corr]<0)(corr_cv_lambda_agg[ntrain,folds,corr] <- 0)
        if(corr_cv_lambda_agg[ntrain,folds,corr]>100)(corr_cv_lambda_agg[ntrain,folds,corr] <- 100)
        MSE_corr_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[corr_cv_lambda_agg[ntrain,folds,corr]+1]
        
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
      
      mean_corr_train_group[ntrain,corr] <- mean_corr
      var_range_train_group[ntrain,corr] <- var_range
      corr_differ_train[ntrain,corr] <- corr_diff
      
      
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
  
  evaluation_data_corr0.15 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                                  test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                                  train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                                  test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                                  test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                                  train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                                  train_low_quantile_corr[,1], train_high_quantile_corr[,1], 
                                                  mean_corr_train_group[,1], var_range_train_group[,1], corr_differ_train[,1],
                                                  cv_lambda_agg[,,1], MSE_cv_on_test_agg[,,1],
                                                  corr_factor[,,1], corr_cv_lambda_agg[,,1], MSE_corr_cv_on_test_agg[,,1],
                                                  which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                                  OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(evaluation_data_corr0.15) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr",
                                            "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                            paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_corr0.3 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                                 test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                                 train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                                 test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                                 test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                                 train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                                 train_low_quantile_corr[,2], train_high_quantile_corr[,2], 
                                                 mean_corr_train_group[,2], var_range_train_group[,2], corr_differ_train[,2],
                                                 cv_lambda_agg[,,2], MSE_cv_on_test_agg[,,2],
                                                 corr_factor[,,2], corr_cv_lambda_agg[,,2], MSE_corr_cv_on_test_agg[,,2],
                                                 which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                                 OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(evaluation_data_corr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                           "test_mean_var", "test_min_var", "test_max_var",
                                           "train_mean_var", "train_min_var", "train_max_var",
                                           "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                           "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                           "train_low_quantile_corr", "train_high_quantile_corr",
                                           "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                           paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_corr0.45 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                                  test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                                  train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                                  test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                                  test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                                  train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                                  train_low_quantile_corr[,3], train_high_quantile_corr[,3], 
                                                  mean_corr_train_group[,3], var_range_train_group[,3], corr_differ_train[,3],
                                                  cv_lambda_agg[,,3], MSE_cv_on_test_agg[,,3],
                                                  corr_factor[,,3], corr_cv_lambda_agg[,,3], MSE_corr_cv_on_test_agg[,,3],
                                                  which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                                  OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(evaluation_data_corr0.45) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                            "test_mean_var", "test_min_var", "test_max_var",
                                            "train_mean_var", "train_min_var", "train_max_var",
                                            "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                            "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                            "train_low_quantile_corr", "train_high_quantile_corr",
                                            "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                            paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                            "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_corr0.6 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                                 test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                                 train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                                 test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                                 test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                                 train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                                 train_low_quantile_corr[,4], train_high_quantile_corr[,4], 
                                                 mean_corr_train_group[,4], var_range_train_group[,4], corr_differ_train[,4],
                                                 cv_lambda_agg[,,4], MSE_cv_on_test_agg[,,4],
                                                 corr_factor[,,4], corr_cv_lambda_agg[,,4], MSE_corr_cv_on_test_agg[,,4],
                                                 which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                                 OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(evaluation_data_corr0.6) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                           "test_mean_var", "test_min_var", "test_max_var",
                                           "train_mean_var", "train_min_var", "train_max_var",
                                           "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                           "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                           "train_low_quantile_corr", "train_high_quantile_corr",
                                           "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                           paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                           "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_corr0.75 <- as.data.frame(cbind(rep(scenario$corr[5],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,5], train_av_var[,5], 
                                                  test_mean_var[,5], test_min_var[,5], test_max_var[,5],
                                                  train_mean_var[,5], train_min_var[,5], train_max_var[,5],
                                                  test_mean_corr[,5], test_median_corr[,5], test_min_corr[,5], test_max_corr[,5],
                                                  test_low_quantile_corr[,5], test_high_quantile_corr[,5], 
                                                  train_mean_corr[,5], train_median_corr[,5], train_min_corr[,5], train_max_corr[,5],
                                                  train_low_quantile_corr[,5], train_high_quantile_corr[,5],
                                                  mean_corr_train_group[,5], var_range_train_group[,5], corr_differ_train[,5],
                                                  cv_lambda_agg[,,5], MSE_cv_on_test_agg[,,5],
                                                  corr_factor[,,5], corr_cv_lambda_agg[,,5], MSE_corr_cv_on_test_agg[,,5],
                                                  which_min_MSE_on_test[,5], min_MSE_on_test[,5],
                                                  OW_MSE_on_test[,5], EW_MSE_on_test[,5]))
  
  colnames(evaluation_data_corr0.75) <-    c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                             "test_mean_var", "test_min_var", "test_max_var",
                                             "train_mean_var", "train_min_var", "train_max_var",
                                             "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                             "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                             "train_low_quantile_corr", "train_high_quantile_corr",
                                             "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                             paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                             "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_corr0.9 <- as.data.frame(cbind(rep(scenario$corr[6],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,6], train_av_var[,6],
                                                 test_mean_var[,6], test_min_var[,6], test_max_var[,6],
                                                 train_mean_var[,6], train_min_var[,6], train_max_var[,6],
                                                 test_mean_corr[,6], test_median_corr[,6], test_min_corr[,6], test_max_corr[,6],
                                                 test_low_quantile_corr[,6], test_high_quantile_corr[,6], 
                                                 train_mean_corr[,6], train_median_corr[,6], train_min_corr[,6], train_max_corr[,6],
                                                 train_low_quantile_corr[,6], train_high_quantile_corr[,6], 
                                                 mean_corr_train_group[,6], var_range_train_group[,6], corr_differ_train[,6],
                                                 cv_lambda_agg[,,6], MSE_cv_on_test_agg[,,6],
                                                 corr_factor[,,6], corr_cv_lambda_agg[,,6], MSE_corr_cv_on_test_agg[,,6],
                                                 which_min_MSE_on_test[,6], min_MSE_on_test[,6],
                                                 OW_MSE_on_test[,6], EW_MSE_on_test[,6]))
  
  colnames(evaluation_data_corr0.9) <-  c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                          "test_mean_var", "test_min_var", "test_max_var",
                                          "train_mean_var", "train_min_var", "train_max_var",
                                          "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                          "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                          "train_low_quantile_corr", "train_high_quantile_corr",
                                          "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                          paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                          "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data <- rbind(evaluation_data_corr0.15, evaluation_data_corr0.3,
                           evaluation_data_corr0.45, evaluation_data_corr0.6,
                           evaluation_data_corr0.75, evaluation_data_corr0.9)
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  evaluation_data <- data.frame(rep("linVar1to1.1", nrow(evaluation_data)), rep(scenario$k, nrow(evaluation_data)), evaluation_data)
  
  colnames(evaluation_data)[1] <- c("VarType")
  colnames(evaluation_data)[2] <- c("NumFC")
  
  
  if(running == 1){
    evaluation_datapart <- evaluation_data}else{
      evaluation_datapart <- rbind(evaluation_datapart, evaluation_data)
    }
  
  # For k = 4,9 and Constant Correlations, the Condition Avoids Duplicates
  
  if(running == 200){
    set.seed(23474545)
  }
  
  print(c("finished running =", running))
}

# table(duplicated(evaluation_datapart))

# save(evaluation_datapart, file ="linearVariance1-1.1_4FCeval.RData")

################################################################################

### 2. Procedure (Differing Correlations)

# Do the Procedure for All Combinations of Parameter Values by 
# Determining Number of Forecasters (4,9) as well as Choosing Function for Sigma 
# (fun_MakeSigmaLin or fun_MakeSigmaExp) and Determining Variance of k-th FC (1.1,
# 1.25, 1.7, 3, 6, 8)
# -> Resulting in 24 Separate Runs

# Modify Sigma Calculation 

fun_MakeSigmaLinDev <- function(this.k, this.corr, this.maxvar){
  
  variance_vector <- rep(NA, this.k)
  variance_vector[1] <- 1.00 # 1-st FC
  variance_vector[this.k] <- this.maxvar # k-th FC
  #  Variances grow linearly from FC 1 to FC k 
  for (v in 1:(this.k-2)) { 
    variance_vector[v+1] <- ( ((this.k-v-1)*variance_vector[1]) + (v*variance_vector[this.k]) )/ (this.k-1)
  }
  
  mu=rep(0,this.k);# Means of Individual Distributions
  
  TempMatrix <- matrix(NA, ncol = this.k, nrow = this.k)
  
  for(i in 1:this.k){
    for(j in 1:this.k){
      if(i <= (this.k+1)/2 & j <= (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr-0.15)
      }else if(i > (this.k+1)/2 & j > (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr+0.15)
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
        TempMatrix[i,j] <- (this.corr-0.15)
      }else if(i > (this.k+1)/2 & j > (this.k+1)/2){
        TempMatrix[i,j] <- (this.corr+0.15)
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



# Determine the Number of Forecasters (4,9)

scenario$k <- 4 # 4,9

scenario$train <- c(25,50,75,100,150,200) 

scenario$source_vector <- rep(1, scenario$k)
source_vector <- scenario$source_vector

scenario$target_vector <- rep(1, scenario$k)
target_vector <- scenario$target_vector

scenario$populationsize <- 20000

granularity <- scenario$no_steps

scenario$corr <- c(0.3,0.45,0.6,0.75) # Pairwise Correlation


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
  
  mean_corr_train_group <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  var_range_train_group <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  corr_differ_train <-  matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  corr_factor <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  corr_cv_lambda_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  MSE_corr_cv_on_test_agg <- array(NA, dim = c(length(scenario$train), length(no_folds), length(scenario$corr)))
  
  which_min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  min_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  EW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  OW_MSE_on_test <- matrix(NA, nrow = length(scenario$train), ncol = length(scenario$corr))
  
  
  # Calculation
  
  
  for(corr in 1:length(scenario$corr)){
    
    # Choose Type of Sigma and Set Value for Highest Variance (1.1, 1.25, 1.7, 3, 6, 8)
    
    Sigma <- fun_MakeSigmaLinDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.1)
    #Sigma <- fun_MakeSigmaExpDev(this.k = scenario$k, this.corr = scenario$corr[corr], this.maxvar = 1.1)
    
    
    # Average Variance Increase
    
    gen_av_var <- (max(diag(Sigma))-min(diag(Sigma)))/(scenario$k-1)
    
    
    for(ntrain in 1:length(scenario$train)){
      
      # Number of Folds (numerical and as character) and 
      # Number (and Percentage) of Training Observations in CV-Calibration Set
      
      no_folds <- c(2, 5, 10, floor(scenario$train[ntrain]/2), scenario$train[ntrain])
      
      n_calib <- c(floor(scenario$train[ntrain]/2), round(0.8*scenario$train[ntrain],0),
                   round(0.9*scenario$train[ntrain],0), (scenario$train[ntrain]-2), (scenario$train[ntrain]-1))
      
      perc_calib <- c(0.5,0.8,0.9,round((scenario$train[ntrain]-2)/scenario$train[ntrain],3), 
                      round((scenario$train[ntrain]-1)/scenario$train[ntrain],3))
      
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
      
      true_low_quantile_corr <- unname(quantile(true_cormat, probs = 0.2, na.rm = TRUE))
      
      true_high_quantile_corr <- unname(quantile(true_cormat, probs = 0.8, na.rm = TRUE))
      
      
      # Estimation of Correlation-Related Variables on Training Data
      
      train_cormat <- cor(ErrorData$Training)
      
      diag(train_cormat) <- NA 
      
      mean_corr_est <- mean(train_cormat, na.rm = TRUE)
      
      median_corr_est <- median(train_cormat, na.rm = TRUE)
      
      min_corr_est <- min(train_cormat, na.rm = TRUE)
      
      max_corr_est <- max(train_cormat, na.rm = TRUE)
      
      low_quantile_corr_est <- unname(quantile(train_cormat, probs = 0.25, na.rm = TRUE))
      
      high_quantile_corr_est <- unname(quantile(train_cormat, probs = 0.75, na.rm = TRUE))
      
      
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
      
      range_var_est <- (max_var_est - min_var_est)
      
      
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
        
        
        # Rename and Generate Variables for Application of Regression Tree
        
        J <- scenario$k
        n <- scenario$train[ntrain]
        K <- no_folds_chr[folds]
        n_cal <- n_calib[folds]
        perc_cal <- perc_calib[folds]
        
        
        if(mean_corr_est < 0.25){
          mean_corr <- "weak"} else if(0.25 <= mean_corr_est & mean_corr_est < 0.55){
            mean_corr <- "moderate"} else if(0.55 <= mean_corr_est & mean_corr_est < 0.75){
              mean_corr <- "strong"} else{
                mean_corr <- "extreme"  
              } 
        
        
        corr_diff <- (high_quantile_corr_est - low_quantile_corr_est)
        
        
        if(range_var_est < 0.45){
          var_range <- "tiny"} else if(0.45 <= range_var_est & range_var_est < 0.95){
            var_range <- "low"} else if(0.95 <= range_var_est & range_var_est < 2.5){
              var_range <- "medium"} else if(2.5 <= range_var_est & range_var_est < 6){
                var_range <- "high"} else{
                  var_range <- "extreme"  
                } 
        
        
        params <- data.frame(J, n, K, n_cal, perc_cal, mean_corr,
                             corr_diff, var_range)
        
        
        # Application of Regression Tree to Determine Correction Factors
        
        corr_factor[ntrain,folds,corr] <- round(predict(treeModelData, params))
        
        # Apply Shrinkage Correction and Determine Respective MSE Value
        
        corr_cv_lambda_agg[ntrain,folds,corr] <- cv_lambda_agg[ntrain,folds,corr] - corr_factor[ntrain,folds,corr]
        if(corr_cv_lambda_agg[ntrain,folds,corr]<0)(corr_cv_lambda_agg[ntrain,folds,corr] <- 0)
        if(corr_cv_lambda_agg[ntrain,folds,corr]>100)(corr_cv_lambda_agg[ntrain,folds,corr] <- 100)
        MSE_corr_cv_on_test_agg[ntrain,folds,corr]  <- MSE_on_test[corr_cv_lambda_agg[ntrain,folds,corr]+1]
        
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
      
      mean_corr_train_group[ntrain,corr] <- mean_corr
      var_range_train_group[ntrain,corr] <- var_range
      corr_differ_train[ntrain,corr] <- corr_diff
      
      
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
  
  evaluation_data_devcorr0.3 <- as.data.frame(cbind(rep(scenario$corr[1],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,1], train_av_var[,1], 
                                                    test_mean_var[,1], test_min_var[,1], test_max_var[,1],
                                                    train_mean_var[,1], train_min_var[,1], train_max_var[,1],
                                                    test_mean_corr[,1], test_median_corr[,1], test_min_corr[,1], test_max_corr[,1],
                                                    test_low_quantile_corr[,1], test_high_quantile_corr[,1], 
                                                    train_mean_corr[,1], train_median_corr[,1], train_min_corr[,1], train_max_corr[,1],
                                                    train_low_quantile_corr[,1], train_high_quantile_corr[,1], 
                                                    mean_corr_train_group[,1], var_range_train_group[,1], corr_differ_train[,1],
                                                    cv_lambda_agg[,,1], MSE_cv_on_test_agg[,,1],
                                                    corr_factor[,,1], corr_cv_lambda_agg[,,1], MSE_corr_cv_on_test_agg[,,1],
                                                    which_min_MSE_on_test[,1], min_MSE_on_test[,1],
                                                    OW_MSE_on_test[,1], EW_MSE_on_test[,1]))
  
  colnames(evaluation_data_devcorr0.3) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                              "test_mean_var", "test_min_var", "test_max_var",
                                              "train_mean_var", "train_min_var", "train_max_var",
                                              "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                              "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                              "train_low_quantile_corr", "train_high_quantile_corr",
                                              "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                              paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_devcorr0.45 <- as.data.frame(cbind(rep(scenario$corr[2],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,2], train_av_var[,2],
                                                     test_mean_var[,2], test_min_var[,2], test_max_var[,2],
                                                     train_mean_var[,2], train_min_var[,2], train_max_var[,2],
                                                     test_mean_corr[,2], test_median_corr[,2], test_min_corr[,2], test_max_corr[,2],
                                                     test_low_quantile_corr[,2], test_high_quantile_corr[,2], 
                                                     train_mean_corr[,2], train_median_corr[,2], train_min_corr[,2], train_max_corr[,2],
                                                     train_low_quantile_corr[,2], train_high_quantile_corr[,2], 
                                                     mean_corr_train_group[,2], var_range_train_group[,2], corr_differ_train[,2],
                                                     cv_lambda_agg[,,2], MSE_cv_on_test_agg[,,2],
                                                     corr_factor[,,2], corr_cv_lambda_agg[,,2], MSE_corr_cv_on_test_agg[,,2],
                                                     which_min_MSE_on_test[,2], min_MSE_on_test[,2],
                                                     OW_MSE_on_test[,2], EW_MSE_on_test[,2]))
  
  colnames(evaluation_data_devcorr0.45) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                               "test_mean_var", "test_min_var", "test_max_var",
                                               "train_mean_var", "train_min_var", "train_max_var",
                                               "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                               "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                               "train_low_quantile_corr", "train_high_quantile_corr",
                                               "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                               paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_devcorr0.6 <- as.data.frame(cbind(rep(scenario$corr[3],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,3], train_av_var[,3], 
                                                    test_mean_var[,3], test_min_var[,3], test_max_var[,3],
                                                    train_mean_var[,3], train_min_var[,3], train_max_var[,3],
                                                    test_mean_corr[,3], test_median_corr[,3], test_min_corr[,3], test_max_corr[,3],
                                                    test_low_quantile_corr[,3], test_high_quantile_corr[,3], 
                                                    train_mean_corr[,3], train_median_corr[,3], train_min_corr[,3], train_max_corr[,3],
                                                    train_low_quantile_corr[,3], train_high_quantile_corr[,3], 
                                                    mean_corr_train_group[,3], var_range_train_group[,3], corr_differ_train[,3],
                                                    cv_lambda_agg[,,3], MSE_cv_on_test_agg[,,3],
                                                    corr_factor[,,3], corr_cv_lambda_agg[,,3], MSE_corr_cv_on_test_agg[,,3],
                                                    which_min_MSE_on_test[,3], min_MSE_on_test[,3],
                                                    OW_MSE_on_test[,3], EW_MSE_on_test[,3]))
  
  colnames(evaluation_data_devcorr0.6) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                              "test_mean_var", "test_min_var", "test_max_var",
                                              "train_mean_var", "train_min_var", "train_max_var",
                                              "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                              "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                              "train_low_quantile_corr", "train_high_quantile_corr",
                                              "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                              paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                              "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  evaluation_data_devcorr0.75 <- as.data.frame(cbind(rep(scenario$corr[4],length(scenario$train)), scenario$train, gen_av_var, test_av_var[,4], train_av_var[,4], 
                                                     test_mean_var[,4], test_min_var[,4], test_max_var[,4],
                                                     train_mean_var[,4], train_min_var[,4], train_max_var[,4],
                                                     test_mean_corr[,4], test_median_corr[,4], test_min_corr[,4], test_max_corr[,4],
                                                     test_low_quantile_corr[,4], test_high_quantile_corr[,4], 
                                                     train_mean_corr[,4], train_median_corr[,4], train_min_corr[,4], train_max_corr[,4],
                                                     train_low_quantile_corr[,4], train_high_quantile_corr[,4], 
                                                     mean_corr_train_group[,4], var_range_train_group[,4], corr_differ_train[,4],
                                                     cv_lambda_agg[,,4], MSE_cv_on_test_agg[,,4],
                                                     corr_factor[,,4], corr_cv_lambda_agg[,,4], MSE_corr_cv_on_test_agg[,,4],
                                                     which_min_MSE_on_test[,4], min_MSE_on_test[,4],
                                                     OW_MSE_on_test[,4], EW_MSE_on_test[,4]))
  
  colnames(evaluation_data_devcorr0.75) <-   c("gen_corr", "ntrain", "gen_av_var", "test_av_var", "train_av_var", 
                                               "test_mean_var", "test_min_var", "test_max_var",
                                               "train_mean_var", "train_min_var", "train_max_var",
                                               "test_mean_corr", "test_median_corr", "test_min_corr", "test_max_corr", "test_low_quantile_corr",
                                               "test_high_quantile_corr", "train_mean_corr", "train_median_corr", "train_min_corr", "train_max_corr",
                                               "train_low_quantile_corr", "train_high_quantile_corr",
                                               "mean_corr_train_group", "var_range_train_group", "corr_differ_train",
                                               paste("ShrinkCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("MSEWithCV_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("Corr_Factor_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("corr_ShrinkCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               paste("corr_MSEWithCVAgg_folds", c("2","5","10","LTO","LOO"), sep = ""),
                                               "TrulyOptShrinkOnTest","MinMSEOnTest","OWMSEOnTest", "EWMSEOnTest")
  
  
  
  evaluation_data <- rbind(evaluation_data_devcorr0.3,evaluation_data_devcorr0.45,
                           evaluation_data_devcorr0.6,evaluation_data_devcorr0.75)
  
  evaluation_data$gen_corr <- paste("Dev",  evaluation_data$gen_corr, sep = "")
  
  
  # Modify Type of Variance Increase and Highest Variance Value
  
  evaluation_data <- data.frame(rep("linVar1to1.1", nrow(evaluation_data)), rep(scenario$k, nrow(evaluation_data)), evaluation_data)
  
  colnames(evaluation_data)[1] <- c("VarType")
  colnames(evaluation_data)[2] <- c("NumFC")
  
  
  if(running == 1){
    evaluation_datapart <- evaluation_data}else{
      evaluation_datapart <- rbind(evaluation_datapart, evaluation_data)
    }
  
  
  
  print(c("finished running =", running))
}


# table(duplicated(evaluation_datapart))

# save(evaluation_datapart, file ="linearVariance1-1.1_4FCDeveval.RData")
