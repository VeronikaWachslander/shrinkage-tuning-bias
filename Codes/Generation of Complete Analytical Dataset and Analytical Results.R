################################################################################
################################################################################
###     Generation of Complete Analytical Dataset and Analytical Results     ###
################################################################################
################################################################################

################################################################################

### Combine Single Datasets to Complete One

################################################################################


source("General Functions and Settings.R")


# Note: The Single Generated Datasets are Contained in Folder "Datasets/Analytical Datasets"

# Set Respective Directory

path <- "Datasets/Analytical Datasets/"

name <- list.files(path = path, pattern = "*.RData", full.names = TRUE)


for(i in 1:length(name)){
  
  analysis_datapart <- list()
  
  load(paste(name[i], sep = ""))
  
  dataset <- analysis_datapart
  
  
  if(i == 1){
    analysis_parts <- dataset
  }else{
    analysis_parts <- rbind(analysis_parts, dataset)
  }
}


colnames(analysis_parts)


# Generate and Mutate Variables

analysis_parts$gen_var_range <- round(analysis_parts$gen_av_var * (analysis_parts$NumFC-1),2)


analysis_parts <- analysis_parts %>% 
  mutate(gen_var_range_group = case_when(
    gen_var_range < 0.45  ~ "tiny",
    0.45 <= gen_var_range  & gen_var_range < 0.95 ~ "low",
    0.95 <= gen_var_range  & gen_var_range < 2.5 ~ "medium",
    2.5 <= gen_var_range  & gen_var_range < 6 ~ "high",
    6 <= gen_var_range ~ "extreme"))



analysis_complete <- pivot_longer(analysis_parts, cols = -c(VarType, NumFC, gen_corr, ntrain, gen_av_var,test_av_var,
                                                            train_av_var, gen_var_range, test_mean_var, test_min_var, test_max_var,
                                                            train_mean_var, train_min_var, train_max_var,
                                                            test_mean_corr, test_median_corr, test_min_corr,
                                                            test_max_corr,test_low_quantile_corr,test_high_quantile_corr,
                                                            train_mean_corr,train_median_corr,train_min_corr,
                                                            train_max_corr,train_low_quantile_corr,train_high_quantile_corr,
                                                            TrulyOptShrinkOnTest,MinMSEOnTest,OWMSEOnTest,EWMSEOnTest,
                                                            gen_var_range_group),
                                  names_to = c(".value", "NumberFolds"), 
                                  names_sep = "_folds")



analysis_complete$overshrinkage <- (analysis_complete$ShrinkCV - analysis_complete$TrulyOptShrinkOnTest)



analysis_complete <- analysis_complete %>%
  mutate(test_var_range = test_av_var * (NumFC-1),
         train_var_range = train_av_var * (NumFC-1))


analysis_complete <- analysis_complete %>% 
  mutate(var_range = case_when( 
    test_var_range < 0.45  ~ "tiny",
    0.45 <= test_var_range  & test_var_range < 0.95 ~ "low",
    0.95 <= test_var_range  & test_var_range < 2.5 ~ "medium",
    2.5 <= test_var_range  & test_var_range < 6 ~ "high",
    6 <= test_var_range ~ "extreme"))


analysis_complete <- analysis_complete %>% 
  mutate(var_range_train = case_when(
    train_var_range < 0.45  ~ "tiny",
    0.45 <= train_var_range  & train_var_range < 0.95 ~ "low",
    0.95 <= train_var_range  & train_var_range < 2.5 ~ "medium",
    2.5 <= train_var_range  & train_var_range < 6 ~ "high",
    6 <= train_var_range ~ "extreme"))


analysis_complete <- analysis_complete %>% 
  mutate(n_cal = case_when(         # Number of Observations in Calibration Sets
    NumberFolds == 2  ~ round(ntrain/2,0),
    NumberFolds == 5 ~ round(0.8*ntrain,0),
    NumberFolds == 10 ~ round(0.9*ntrain,0),
    NumberFolds == "LTO" ~ (ntrain-2),
    NumberFolds == "LOO" ~ (ntrain-1)))


analysis_complete <- analysis_complete %>% 
  mutate(perc_cal = case_when(   # Proportion of Training Set in Calibration Sets
    NumberFolds == 2  ~ 0.5,
    NumberFolds == 5 ~ 0.8,
    NumberFolds == 10 ~ 0.9,
    NumberFolds == "LTO" ~ round((ntrain-2)/ntrain,3),
    NumberFolds == "LOO" ~ round((ntrain-1)/ntrain,3)))


# Correlation Groups and Differences

analysis_complete <- analysis_complete %>% 
  mutate(mean_corr = case_when(           
    test_mean_corr < 0.25  ~ "weak",
    0.25 <= test_mean_corr  & test_mean_corr < 0.55 ~ "moderate",
    0.55 <= test_mean_corr  & test_mean_corr < 0.75 ~ "strong",
    0.75 <= test_mean_corr ~ "extreme"))

analysis_complete <- analysis_complete %>% 
  mutate(mean_corr_train = case_when(
    train_mean_corr < 0.25  ~ "weak",
    0.25 <= train_mean_corr  & train_mean_corr < 0.55 ~ "moderate",
    0.55 <= train_mean_corr  & train_mean_corr < 0.75 ~ "strong",
    0.75 <= train_mean_corr ~ "extreme"))


analysis_complete <- analysis_complete %>%  
  mutate(corr_diff = round(test_max_corr - test_min_corr,1))


analysis_complete <- analysis_complete %>% 
  mutate(corr_diff_train = round(train_max_corr - train_min_corr,1))


# save(analysis_parts, file ="final_analysis_parts.RData")

# save(analysis_complete, file ="final_analysis_complete.RData")


################################################################################

### Analytical Results of Paper

################################################################################

# Load the Generated Complete Analytical Dataset (if not done yet)

# path <- "~/"

# load("~/final_analysis_complete.RData")

# load("~/final_analysis_parts.RData")  # Needed for Generation of Figure 1

dim(analysis_complete) # Over 13 Million Cases

mean(analysis_complete$overshrinkage, na.rm = TRUE)


# Analytical Results of Table 1 in Paper (meanOverShrink)

# Delta Rho < 0.1

View(analysis_complete %>%
       filter(NumFC %in% c(5,15)) %>%
       filter(gen_corr %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) %>%
       filter(NumberFolds %in% c(2,5,"LOO")) %>%
       group_by(var_range, NumFC, mean_corr, NumberFolds) %>%
       summarise(meanOverShrink = round(mean(overshrinkage, na.rm = TRUE),2)))


# Delta Rho >= 0.1

View(analysis_complete %>%
       filter(NumFC %in% c(5,15)) %>%
       filter(gen_corr %in% c("Dev0.2","Dev0.3","Dev0.4","Dev0.5","Dev0.6","Dev0.7","Dev0.8")) %>%
       filter(NumberFolds %in% c(2,5,"LOO")) %>%
       group_by(var_range, NumFC, mean_corr, NumberFolds) %>%
       summarise(meanOverShrink = round(mean(overshrinkage, na.rm = TRUE),2)))


################################################################################

# Generation of Figure 1 in Paper (Labeling Slightly Modified)

dataset_graphic <- analysis_parts %>%
  filter(gen_corr %in% c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) %>%
  filter(ntrain %in% c(20, 30, 40, 50, 60, 70, 80, 100,125,150,175,200)) %>%
  select(NumFC, gen_var_range_group, ntrain, gen_corr, ShrinkCV_folds5, TrulyOptShrinkOnTest)


corr_weak <- dataset_graphic %>%
  filter(gen_corr %in% c(0.1, 0.2)) %>%
  group_by(gen_var_range_group, ntrain) %>%
  summarise(meanCVShrink = mean(ShrinkCV_folds5), meanOptShrink = mean(TrulyOptShrinkOnTest))

corr_moderate <- dataset_graphic %>%
  filter(gen_corr %in% c(0.3, 0.4, 0.5)) %>%
  group_by(gen_var_range_group, ntrain) %>%
  summarise(meanCVShrink = mean(ShrinkCV_folds5), meanOptShrink = mean(TrulyOptShrinkOnTest))

corr_strong <- dataset_graphic %>%
  filter(gen_corr %in% c(0.6, 0.7)) %>%
  group_by(gen_var_range_group, ntrain) %>%
  summarise(meanCVShrink = mean(ShrinkCV_folds5), meanOptShrink = mean(TrulyOptShrinkOnTest))

corr_extreme <- dataset_graphic %>%
  filter(gen_corr %in% c(0.8, 0.9)) %>%
  group_by(gen_var_range_group, ntrain) %>%
  summarise(meanCVShrink = mean(ShrinkCV_folds5), meanOptShrink = mean(TrulyOptShrinkOnTest))


gg_design <- list(geom_line(aes(color="Var1.2", linetype = "CV"), linewidth = 1),
                  geom_line(aes(ntrain, meanCVShrink_low, color="Var1.5", linetype = "CV"), linewidth = 1),
                  geom_line(aes(ntrain, meanCVShrink_medium, color="Var2", linetype = "CV"), linewidth = 1),
                  geom_line(aes(ntrain, meanCVShrink_high, color="Var4", linetype = "CV"), linewidth = 1),
                  geom_line(aes(ntrain, meanCVShrink_extreme, color="Var9", linetype = "CV"), linewidth = 1),
                  geom_line(aes(ntrain, meanOptShrink_tiny, color="Var1.2", linetype = "TrulyOpt"), linewidth = 1),
                  geom_line(aes(ntrain, meanOptShrink_low, color="Var1.5", linetype = "TrulyOpt"), linewidth = 1),
                  geom_line(aes(ntrain, meanOptShrink_medium, color="Var2", linetype = "TrulyOpt"), linewidth = 1),
                  geom_line(aes(ntrain, meanOptShrink_high, color="Var4", linetype = "TrulyOpt"), linewidth = 1),
                  geom_line(aes(ntrain, meanOptShrink_extreme, color="Var9", linetype = "TrulyOpt"), linewidth = 1),
                  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100)),
                  scale_x_continuous(breaks = c(20,40,60,80,100,120,140,160,180,200), limits = c(20,200)),
                  scale_linetype_manual(values = c(CV = "solid", TrulyOpt = "dashed"), name = "Shrinkage",
                                        labels = c(TrulyOpt = bquote(italic("\u03BB"[true]^"*")), CV = bquote(italic("\u03BB"[CV]^"*"))),
                                        limits = c("CV", "TrulyOpt")),
                  scale_color_grey(start = 0.85, end = 0.1, labels = c(Var1.2 = "tiny", Var1.5 = "low", Var2 = "medium", Var4 = "high", Var9 = "extreme")),
                  theme(legend.position = "none"),
                  theme_bw(),
                  theme(text = element_text(size = 12), #change font size of all text
                        axis.text = element_text(size = 12), #change font size of axis text
                        axis.title = element_text(size = 12), #change font size of axis titles
                        axis.title.x.bottom = element_blank(), #change font size of axis titles
                        axis.title.y.left = element_blank(), #change font size of axis titles
                        plot.title = element_text(size = 12, hjust = 0.5), #change font size of plot title
                        plot.subtitle = element_text(size = 12, hjust = 0.5), #change font size of plot title
                        plot.margin = margin(0.1,0.3,0.6,0.3, 'cm'), #c(top, right, bottom, left)
                        legend.text = element_text(size = 11), #change font size of legend text
                        legend.title = element_text(size = 11))) #change font size of legend title 


corr_weak <- corr_weak %>% 
  pivot_wider(names_from = gen_var_range_group, values_from = c(meanCVShrink,meanOptShrink))

corr_moderate <- corr_moderate %>% 
  pivot_wider(names_from = gen_var_range_group, values_from = c(meanCVShrink,meanOptShrink))

corr_strong <- corr_strong %>% 
  pivot_wider(names_from = gen_var_range_group, values_from = c(meanCVShrink,meanOptShrink))

corr_extreme <- corr_extreme %>% 
  pivot_wider(names_from = gen_var_range_group, values_from = c(meanCVShrink,meanOptShrink))


corr_weak  <- ggplot(corr_weak, aes(x = ntrain, y = meanCVShrink_tiny, showlabels = TRUE )) +
  ggtitle(expression(italic("\u03C1 = weak"))) +
  gg_design 

corr_moderate  <- ggplot(corr_moderate, aes(x = ntrain, y = meanCVShrink_tiny, showlabels = TRUE )) +
  ggtitle(expression(italic("\u03C1 = moderate"))) +
  gg_design


corr_strong  <- ggplot(corr_strong, aes(x = ntrain, y = meanCVShrink_tiny, showlabels = TRUE )) +
  ggtitle(expression(italic("\u03C1 = strong"))) +
  gg_design

corr_extreme  <- ggplot(corr_extreme, aes(x = ntrain, y = meanCVShrink_tiny, showlabels = TRUE )) +
  ggtitle(expression(italic("\u03C1 = extreme"))) +
  gg_design


library(ggpubr)

ggarrange(corr_weak, corr_moderate, corr_strong, corr_extreme, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(left = text_grob("Shrinkage (%)", rot = 90, vjust = 1, hjust = 0.2),
                  bottom = text_grob("n", hjust = -1, vjust = -4.5)) 

# 10 x 7

