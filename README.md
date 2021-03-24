# dom_adapt
Code for the NIPS 2018 paper "Domain Adaptation by Using Causal Inference to Predict Invariant Conditional Distributions"


# Installation instructions

The code expects a Clingo solver in the folder ASP (point 1) and a few R packages (point 2):

### 1. Download the new clingo

  * From http://sourceforge.net/projects/potassco/files/clingo and rename it as "clingo"
  * To see whether your clingo installation is working, you can try running:
    ```
    ./ASP/clingo
    ```
### 2. Install the R packages(*)

  * 2.1 Install Bioconductor packages
    * With R *version 3.5 or greater*, use *BiocManager*:
    ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(version = "3.12")
    
    BiocManager::install(c('graph','RBGL','gmp','RcppArmadillo'))
    ```
    * With R *version lesser than 3.5*, install:
    ```
    source('http://bioconductor.org/biocLite.R')
    biocLite(c('graph','RBGL','gmp','RcppArmadillo'))
    ```
  * 2.2 Install the remaining packages:
  
    ```
    install.packages(c('deal','combinat','hash','bnlearn','foreach','doMC','caTools','expm'))
    install.packages(c('pcalg'))
    ```
    Note that:
    * Installing 'pcalg' may require you to install also a few other packages (e.g. robustbase, ggm).
    * For R version 4.0.4 'pcalg' dependencies are managed automatically.

### 3. Start R

  * Navigate to the R/ directory and run:
    ```
    source('load.R')
    loud()
    ```
(*) If you want to keep your global environment unchanged, please, consider using the renv package (https://rstudio.github.io/renv/).
