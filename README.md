# PvCompAdv R package
### Package containing R scripts and data used in draft paper currently titled *Potential gains of long-distance trade in electricity*.

A preprint of a previous version of the paper can be found [here](https://ui.adsabs.harvard.edu/abs/2022arXiv220501436L/abstract)

Status of the paper and corresponding links to text will be updated as review process proceeds.

### Running the package
The package is meant to be used in combination with a `results_replication.R` script (not provided here until paper is published.

The package can be installed in an R script as follows: 
```R
install.packages('devtools')
library(devtools)
install_github("kawilliges/PvCompAdv")
```

## Instructions for replication of results
Due to technical issues, installation and execution of replication of an R package allowing for streamlined replication of results is in progress. In order to avoid delaying resubmission of the mansucript, and inconveniencing reviewers, we have adpoted a temporary fix to allow for replication of our manuscript's results. Please forgive the inconvenience; we will implement a more streamlined solution as soon as possible. 

In order to replicate the results found in the manuscript, an installed copy of RStudio is required. Additionally, [installation of RTools](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html) is required, to install the package developed for the manuscript from source code.

1. Install RTools, a set of programs required on Windows to build R packages from source, if necessary
  * for installation instructions, see [this installation page](https://ohdsi.github.io/Hades/rSetup.html), subheading "Installing RTools"
2. Open the R script `results_replication.R`
3. Install required libraries by running the R environment setup code found at the beginning of the file (Run the script to line XX)
4. Create a working directory, e.g. `C:\Github\pv_trade_paper_replication` and specify its location in the R script (line XX)
5. Create a subfolder `Data` within the working directory 
6. Download the contents of the `Data` folder of this repository, and place it within the `Data` folder of the working directory
7. You should now be able to execute the remainder of the R script