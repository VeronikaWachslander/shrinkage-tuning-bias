# Weight Shrinkage Tuning Bias

This repository contains the source code for the paper *Analyzing and Correcting Biased Machine Learning-Based Tuning of Weight Shrinkage in Forecast Combination*. 

The structure of the repository is as follows: The folder **Codes** contains all source codes for the data generation and the reproduction of all results shown in the paper. The respective files are:

- General Functions and Settings
-> Required Functions and Settings. 
  
- Generation of Analytical Datasets
-> Code for Generation of Synthetic Datasets Used to Analyze Shrinkage Tuning Bias (generated Datasets contained in Datasets/Analytical Datasets).
  
- Generation of Complete Analytical Dataset and Analytical Results
-> Code to Combine Single Analytical Datasets and Reproduce the Analytical Results.
  
- Regression Tree Model
-> Code for Generation and Tuning of Regression Tree to Predict Shrinkage Bias.
  
- Generation of Evaluation Datasets
-> Code for Generation of Synthetic Datasets Used to Evaluate Shrinkage Correction (generated Datasets contained in Datasets/Evaluation Datasets).

  
- Generation of Complete Evaluation Dataset and Evaluation Results
-> Code to Combine Single Evaluation Datasets and Reproduce the Evaluation Results.


The folder **Datasets** contains the generated synthetic datasets, with two subfolders:

- Analytical Datasets
-> All 140 Datasets generated for Analysis of Shrinkage Tuning Bias.

- Evaluation Datasets
-> All 48 Datasets generated for Evaluation of Shrinkage Correction. 
