# Weight Shrinkage Tuning Bias

This repository contains the source code for the paper *Analyzing and Correcting Biased Machine Learning-Based Tuning of Weight Shrinkage in Forecast Combination*.
% which is accepted for presentation at the *European Conference on Machine Learning and Principles and Practice of Knowledge Discovery in Databases (ECML PKDD) 2025* in Porto, Portugal. 

The repository contains two folders **Codes** and **Datasets**, which are structured as follows: <br>

<ul>
<li> <b> Codes </b> <i> (containing all source codes for data generation & reproduction of results shown in the paper) </i> <br> <br>
    <ul>
    <li> <b> General Functions and Settings.R </b> 
    <i> (required functions & settings) </i>  </li> <br>
    <li> <b> Generation of Analytical Datasets.R </b> 
    <i> (generating synthetic datasets to analyze shrinkage tuning bias) </i> </li> <br>
    <li> <b> Generation of Complete Analytical Dataset and Analytical Results.R </b> 
    <i> (combining analytical datasets & reproducing analytical results) </i> </li> <br>
    <li> <b> Regression Tree Model.R </b> 
    <i> (generating & tuning regression tree to predict shrinkage bias) </i> </li> <br>
    <li> <b> Generation of Evaluation Datasets.R </b> 
    <i> (generating synthetic datasets to evaluate shrinkage correction) </i> </li> <br>
    <li> <b> Generation of Complete Evaluation Dataset and Evaluation Results.R </b> 
    <i> (combining evaluation datasets & reproducing evaluation results) </i> </li>
    </ul> <br>
<li> <b> Datasets </b> <i> (containing generated synthetic datasets) </i> <br> <br>
   <ul>
   <li> <b> Analytical Datasets </b> 
   <i> (140 datasets generated with <b> Generation of Analytical Datasets.R </b>) </i> </li> <br>
   <li> <b> Evaluation Datasets </b>
   <i> (48 datasets generated with <b> Generation of Evaluation Datasets.R </b>) </i> </li> <br>
        <li> <b> Regression Dataset </b> 
   <i> (1 dataset generated with <b> Regression Tree Model.R </b>) </i> </li>
</li>
</ul>
