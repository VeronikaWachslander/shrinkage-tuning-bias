# Weight Shrinkage Tuning Bias

This repository contains the source code for the paper *Analyzing and Correcting Biased Machine Learning-Based Tuning of Weight Shrinkage in Forecast Combination*. 

The repository contains two folders **Codes** and **Datasets**, which are structured as follows.

<ul>
<li> <b> Codes </b> <i> (containing all source codes for data generation & reproduction of results shown in the paper) </i> <br> 
    <ul>
    <li> <b> General Functions and Settings </b> 
    <i> (required functions and settings) </i>  </li> 
    <li> <b> Generation of Analytical Datasets </b> 
    <i> (generating synthetic datasets to analyze shrinkage tuning bias) </i> </li>
    <li> <b> Generation of Complete Analytical Dataset and Analytical Results </b> 
    <i> (combining analytical datasets & reproducing analytical results) </i> </li>
    <li> <b> Regression Tree Model </b> 
    <i> (generating and tuning regression tree to predict shrinkage bias) </i> </li>
    <li> <b> Generation of Evaluation Datasets </b> 
    <i> (generating synthetic datasets to evaluate shrinkage correction) </i> </li>
    <li> <b> Generation of Complete Evaluation Dataset and Evaluation Results </b> 
    <i> (combining evaluation datasets & reproducing evaluation results) </i> </li>
    </ul>
<li> <b> Datasets </b> <i> (containing generated synthetic datasets) </i>
   <ul>
   <li> <b> Analytical Datasets </b> </li>
   <li> <b> Evaluation Datasets </b> </li>
</li>
</ul>


The folder **Datasets** contains the generated synthetic datasets, with two subfolders:

- **Analytical Datasets**
-> All 140 Datasets generated for Analysis of Shrinkage Tuning Bias.

- **Evaluation Datasets**
-> All 48 Datasets generated for Evaluation of Shrinkage Correction. 
