# RBMR：A Two-Sample Robust Bayesian Mendelian Randomization Method Accounting for Linkage Disequilibrium and Idiosyncratic Pleiotropy
**RBMR** is a Mendelian randomization (MR) method to estimate the causal effect by accounting for the LD structure, systematic pleiotropy and idiosyncratic pleiotropy simultaneously in a unified framework. 
## Installation
To install the development version of **RBMR**, it's easiest to use the 'devtools' package.

Use the following command in R to install the package:
```
#install.packages("devtools")
library(devtools)
install_github("AnqiWang2021/RBMR/RBMR")
```
## Usage
The ['RBMR' Example](https://github.com/AnqiWang2021/RBMR/blob/main/Example/example.R) will provide a good example for two-sample Mendelian randomization analysis using **RBMR** package. This example illustrates the implements of **RBMR** for real data analysis. The
following datasets (‘heart attack myocardial infarction.txt’, ‘exposure c4d.txt’, ‘outcome cardiogram.txt’, ‘all_chr_1000G.bed’,
‘all_chr_1000G.fam’, ‘all_chr_1000G.bim’, ‘fourier_ls-all.bed’) should be prepared. Download here (https://drive.google.com/drive/folders/19iVojq-XbRzzrAwWIGomUW2o-p3lcuLp?usp=sharing).

