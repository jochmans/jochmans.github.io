CODES FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner


The data set is data_research.dta, from Ductor, Fafchamps, Goyal and van der Leij (2014, Review of Economics and Statistics) 

The data can be freely downloaded from the stable archive: https://direct.mit.edu/rest/article/96/5/936/58200/Social-Networks-and-Research-Output 


The codes are run in stata (version 19.5) and matlab (version R2024b), including the commands parfor (from Parallel Computing), fminunc (from Optimization) and quantile (from Statistics and Machine Learning)


First run the stata code Code_descriptives.do

To do so, you will need to enter your path on line 10 of the code

All files need to be in the same directory for the codes to run properly

This code produces the descriptive statistics reported in the text


Based on this, process the data for matlab using the matlab code Code_prepare_data


Then, run the matlab code Code_estimation

This produces:

(1) The point estimates for Table 1

(2) The point estimates for Table 2

(3) Figure 1 in the supplement

(4) Figure 2 in the supplement

Time (Intel(R) Core(TM) i9-9900K CPU @ 3.60GHz   3.60 GHz; parallel loop across 8 workers): 3h30


Then, run the matlab code Code_q_diagnostic

This produces the numbers in the last column of Table 1

Time: 1h15


Then, run the matlab code Code_parametric_bootstrap

This produces the standard errors in Tables 1 and 2

Time: approximately 40 hours


In addition, the matlab code Code_simulation produces the numbers in Tables 3 and 4 in the supplement

Time: 6h30


And the matlab code Code_GMM produces the GMM estimates reported in the supplement

Time: 2min:30 (data) and 7min:20 (simulation) 



For questions, contact Stephane Bonhomme at sbonhomme@uchicago.edu

