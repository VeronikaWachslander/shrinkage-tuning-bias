################################################################################
################################################################################
###        Regression Tree for Predicting the Shrinkage Tuning Bias          ###
################################################################################
################################################################################

### Generating and Tuning the Regression Tree 

################################################################################

source("General Functions and Settings.R")

# Load the Generated Complete Analytical Dataset (Resulting from "Generation of Complete 
# Analytical Dataset and Analytical Results.R")

path <- "~/"

load("~/final_analysis_complete.RData")

# Rename Variables

analysis_complete <- analysis_complete %>%
  rename(n = ntrain,
         K = NumberFolds,
         J = NumFC)


# Additional Packages Required for Regression Tree

library(rpart)
library(rpart.plot)

#install.packages("mlr")
library(mlr)


# Tuning the Tree

regression_complete <- subset(analysis_complete, select = c(overshrinkage,J, n, K, n_cal, perc_cal,
                                                            var_range, mean_corr, corr_diff))

regression_complete <- regression_complete %>%
  drop_na(overshrinkage)


# save(regression_complete, file ="final_regression_complete.RData")


overshrinkTib <- mutate_if(regression_complete, is.character, as.factor)

overshrinkTask <- makeRegrTask(data = overshrinkTib, target = "overshrinkage")

tree <- makeLearner("regr.rpart")

getParamSet(tree)

treeParamSpace <- makeParamSet(
  makeDiscreteParam("cp", values = seq(0.00000001,0.0000005,0.00000001)))

control.grid <- makeTuneControlGrid() 

cvForTuning <- makeResampleDesc("CV", iters = 5)

tunedTreePars <- tuneParams(tree, task = overshrinkTask,
                            resampling = cvForTuning,
                            par.set = treeParamSpace,
                            control = control.grid)

# Best Value for cp: 2.8 x 10^-7

tunedTree <- setHyperPars(tree, par.vals = tunedTreePars$x)

################################################################################

### Learning the Regression Tree 

################################################################################

# Load the Generated Regression Dataset and Train the Regression Tree (if not done yet)

load("Datasets/Regression Dataset/final_regression_complete.RData")

dim(regression_complete)

overshrinkTib <- mutate_if(regression_complete, is.character, as.factor)

overshrinkTask <- makeRegrTask(data = overshrinkTib, target = "overshrinkage")

tree <- makeLearner("regr.rpart")

tunedTree <- setHyperPars(tree, par.vals = list(cp = 0.00000028))

tunedTreeModel <- train(tunedTree, overshrinkTask)

treeModelData <- getLearnerModel(tunedTreeModel)

rpart.plot(treeModelData, roundint = FALSE)


# Tree Depth

trunc(log(max(as.integer(row.names(treeModelData$frame))), base = 2))


# Generation of Figure 2 in Paper (Slightly Modified in Paper due to Conciseness)

tunedTreePaper <- setHyperPars(tree, par.vals = list(maxdepth = 5, cp = 0.00000028))

tunedTreeModelPaper <- train(tunedTreePaper, overshrinkTask)

treeModelDataPaper <- getLearnerModel(tunedTreeModelPaper)

rpart.plot(treeModelDataPaper, roundint = FALSE, digits = 4)

