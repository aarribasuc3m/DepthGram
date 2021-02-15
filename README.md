# DepthGram
THESE ARE THE FILES ACCOMPANYING THE PAPER: "Depthgram: Visualizing Outliers in High
Dimensional Functional Data with application to Task fMRI data exploration"

List of files:

1. depthGram.R: R function to compute the depthgrams configurations (not to plot them)

2. depthGramPlot.R: R function to plot depthgrams configurations

3. Examples.R: Examples of the use of 1 and 2 in real data sets

4. Figures_Code.R: R code to reproduce all the figures in the paper

5. other_functions.R: auxiliary functions to plot FOM and ms-plot matching aesthetics with DepthGram

6. OutDetectionRule_DG.R: R code for a basic outlier detection rule for the DepthGram

7. Sim_MultFunData.R: R function to generate multivariate functional data 
with dependence across dimensions and outliers

8. Simulations_LowDim.R: R code for the low-dimensional simulation study in the paper

9. Simulations.R: R code for the high-dimensional simulation study in the paper

10. Comparison_ComputingTimes.R: R code to compare computing times of DepthGram, MS-plot and FOM (as detailed in the Supplementary Materials)

List of folders needed:

1. msplot_code: Code for ms-plot . To be downloaded from  from  https://www.tandfonline.com/doi/suppl/10.1080/10618600.2018.1473781?scroll=top, supplement.rar

2. Sample_Results_Files: sample output files of the simulation study to be used in Figures_Code.R (to reproduce Figure 6 and similar ones). Available under request (large RData files)
