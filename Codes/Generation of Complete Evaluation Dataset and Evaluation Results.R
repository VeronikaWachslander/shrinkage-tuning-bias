################################################################################
################################################################################
###     Generation of Complete Evaluation Dataset and Evaluation Results     ###
################################################################################
################################################################################

################################################################################

### Combine Single Datasets to Complete One

################################################################################


source("General Functions and Settings.R")


# Note: The Single Generated Datasets are Contained in Folder "Datasets/Evaluation Datasets"

# Set Respective Directory

path <- "~/"

name <- list.files(path = path, pattern = "*.RData", full.names = TRUE)


for(i in 1:length(name)){
  
  evaluation_datapart <- list()
  
  load(paste(name[i], sep = ""))
  
  dataset <- evaluation_datapart
  
  
  if(i == 1){
    evaluation_parts <- dataset
  }else{
    evaluation_parts <- rbind(evaluation_parts, dataset)
  }
}


colnames(evaluation_parts)

# Generate and Mutate Variables

evaluation_complete <- pivot_longer(evaluation_parts, cols = -c(VarType, NumFC, gen_corr, ntrain, gen_av_var,test_av_var,
                                                                train_av_var, test_mean_var, test_min_var, test_max_var,
                                                                train_mean_var, train_min_var, train_max_var,
                                                                test_mean_corr, test_median_corr, test_min_corr,
                                                                test_max_corr,test_low_quantile_corr,test_high_quantile_corr,
                                                                train_mean_corr,train_median_corr,train_min_corr,
                                                                train_max_corr,train_low_quantile_corr,train_high_quantile_corr,
                                                                mean_corr_train_group, var_range_train_group, corr_differ_train,
                                                                TrulyOptShrinkOnTest,MinMSEOnTest,OWMSEOnTest,EWMSEOnTest),
                                    names_to = c(".value", "NumberFolds"), 
                                    names_sep = "_folds")


evaluation_complete$overshrinkage <- (as.numeric(evaluation_complete$ShrinkCV) - as.numeric(evaluation_complete$TrulyOptShrinkOnTest))

evaluation_complete$corrovershrinkage <- (as.numeric(evaluation_complete$corr_ShrinkCVAgg) - as.numeric(evaluation_complete$TrulyOptShrinkOnTest))



evaluation_complete <- evaluation_complete %>%
  mutate(test_var_range = as.numeric(test_av_var) * (as.numeric(NumFC)-1))


evaluation_complete <- evaluation_complete %>% 
  mutate(var_range = case_when( 
    test_var_range < 0.45  ~ "tiny",
    0.45 <= test_var_range  & test_var_range < 0.95 ~ "low",
    0.95 <= test_var_range  & test_var_range < 2.5 ~ "medium",
    2.5 <= test_var_range  & test_var_range < 6 ~ "high",
    6 <= test_var_range ~ "extreme"))


# Correlation Groups and Differences

evaluation_complete <- evaluation_complete %>% 
  mutate(mean_corr = case_when(          
    test_mean_corr < 0.25  ~ "weak",
    0.25 <= test_mean_corr  & test_mean_corr < 0.55 ~ "moderate",
    0.55 <= test_mean_corr  & test_mean_corr < 0.75 ~ "strong",
    0.75 <= test_mean_corr ~ "extreme"))

# save(evaluation_complete, file ="final_evaluation_complete.RData")


evaluation_parts <- evaluation_parts %>%
  mutate(test_var_range = as.numeric(test_av_var) * (as.numeric(NumFC)-1))


evaluation_parts <- evaluation_parts %>% 
  mutate(var_range = case_when( 
    test_var_range < 0.45  ~ "tiny",
    0.45 <= test_var_range  & test_var_range < 0.95 ~ "low",
    0.95 <= test_var_range  & test_var_range < 2.5 ~ "medium",
    2.5 <= test_var_range  & test_var_range < 6 ~ "high",
    6 <= test_var_range ~ "extreme"))



# Correlation Groups and Differences

evaluation_parts <- evaluation_parts %>% 
  mutate(mean_corr = case_when(          
    test_mean_corr < 0.25  ~ "weak",
    0.25 <= test_mean_corr  & test_mean_corr < 0.55 ~ "moderate",
    0.55 <= test_mean_corr  & test_mean_corr < 0.75 ~ "strong",
    0.75 <= test_mean_corr ~ "extreme"))


# save(evaluation_parts, file ="final_evaluation_parts.RData")

################################################################################

### Evaluation Results of Paper

################################################################################

# Load the Generated Complete Dataset 

path <- "~/"

load("~/final_evaluation_complete.RData")

load("~/final_evaluation_parts.RData")


# Select Parameter Values Presented in Paper

evaluation_complete <- evaluation_complete %>%
  filter(NumberFolds %in% c(2,5,"LOO")) %>%
  filter(ntrain %in% c(25,50,100,200))

dim(evaluation_complete) # 720,000 Cases



evaluation_parts <- evaluation_parts %>%
  filter(ntrain %in% c(25,50,100,200))




# Evaluation Results of Table 3 in Paper (CVMSEPerc)

View(evaluation_complete %>%
       group_by(ntrain, NumberFolds, mean_corr, var_range) %>%
       summarise(CVMSEPerc = round(100*(mean(as.numeric(corr_MSEWithCVAgg)) - mean(as.numeric(MSEWithCV)))/mean(as.numeric(MSEWithCV)),2),
                 corr_MSECV = mean(as.numeric(corr_MSEWithCVAgg)), MSECV = mean(as.numeric(MSEWithCV))))


# MSE Comparison Differentiated by Number of Folds

MSE_Comparison <- evaluation_parts %>%
  filter(ntrain %in% c(25,50,100,200)) %>%
  group_by(ntrain, mean_corr, var_range) %>%
  summarise(CVMSEPerc2 = round(100*(mean(as.numeric(corr_MSEWithCVAgg_folds2)) - mean(as.numeric(MSEWithCV_folds2)))/mean(as.numeric(MSEWithCV_folds2)),4),
            CVMSEPerc5 = round(100*(mean(as.numeric(corr_MSEWithCVAgg_folds5)) - mean(as.numeric(MSEWithCV_folds5)))/mean(as.numeric(MSEWithCV_folds5)),4),
            CVMSEPercLOO = round(100*(mean(as.numeric(corr_MSEWithCVAgg_foldsLOO)) - mean(as.numeric(MSEWithCV_foldsLOO)))/mean(as.numeric(MSEWithCV_foldsLOO)),4),
            corr_MSECV2 = mean(as.numeric(corr_MSEWithCVAgg_folds2)), MSECV2 = mean(as.numeric(MSEWithCV_folds2)),
            corr_MSECV5 = mean(as.numeric(corr_MSEWithCVAgg_folds5)), MSECV5 = mean(as.numeric(MSEWithCV_folds5)),
            corr_MSECVLOO = mean(as.numeric(corr_MSEWithCVAgg_foldsLOO)), MSECVLOO = mean(as.numeric(MSEWithCV_foldsLOO)),
            MSECheck = (min(corr_MSECVLOO,MSECVLOO) < min(corr_MSECV2,MSECV2,corr_MSECV5,MSECV5)))


table(MSE_Comparison$MSECheck) # In 69 of 80 Evaluated Scenarios, LOO (Corrected or Uncorrected) Results in Lowest Mean MSE.


# Characteristics of Scenarios For Which LOO Does Not Result in Lowest Mean MSE:

MSE_Comparison[MSE_Comparison$MSECheck == FALSE, c("ntrain","mean_corr","var_range")]


# Comparison of Uncorrected 2-Fold CV (as this Method is Mostly the Best One in These Scenarios)
# and Uncorrected LOO for These Scenarios (in %):

CV_MSE_2vsLOO <- evaluation_complete %>%
  filter(NumberFolds %in% c(2, "LOO")) %>%
  filter(var_range == "tiny" & ntrain == "25"|
           var_range == "tiny" & mean_corr == "weak"|
           var_range == "low" & ntrain == "25" & mean_corr %in% c("weak", "moderate")|
           var_range == "low" & ntrain == "50" & mean_corr == "weak"|
           var_range == "medium" & ntrain == "200" & mean_corr == "extreme") %>%
  group_by(NumberFolds) %>%
  summarise(CVMSE = mean(as.numeric(MSEWithCV)))


round(100 * (CV_MSE_2vsLOO[CV_MSE_2vsLOO$NumberFolds == "2", "CVMSE"]/CV_MSE_2vsLOO[CV_MSE_2vsLOO$NumberFolds == "LOO", "CVMSE"] - 1),2)



# Comparison of Mean MSE of LOO with 2- and 5-Fold CV for Larger Amounts of Training Observations (in %)

as_tibble(evaluation_parts) %>%
  filter(ntrain %in% c("100", "200")) %>%
  filter(var_range %in% c("medium", "high", "extreme")) %>%
  group_by(ntrain, var_range, mean_corr) %>%
  summarise(corr_MSECV2 = mean(as.numeric(corr_MSEWithCVAgg_folds2)), MSECV2 = mean(as.numeric(MSEWithCV_folds2)),
            corr_MSECV5 = mean(as.numeric(corr_MSEWithCVAgg_folds5)), MSECV5 = mean(as.numeric(MSEWithCV_folds5)),
            corr_MSECVLOO = mean(as.numeric(corr_MSEWithCVAgg_foldsLOO)), MSECVLOO = mean(as.numeric(MSEWithCV_foldsLOO)),
            minCV2 = min(corr_MSECV2,MSECV2),minCV5 = min(corr_MSECV5,MSECV5),minCVLOO = min(corr_MSECVLOO,MSECVLOO)) %>%
  ungroup() %>%
  summarise(comp2LOO = 100*(mean(as.numeric(minCV2))/mean(as.numeric(minCVLOO))-1),
            comp5LOO = 100*(mean(as.numeric(minCV5))/mean(as.numeric(minCVLOO))-1))
