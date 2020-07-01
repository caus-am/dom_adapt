# dom_adapt
Code for the NIPS 2018 paper "Domain Adaptation by Using Causal Inference to Predict Invariant Conditional Distributions"



# Installation instructions:

The code expects a Clingo solver in the folder ASP (point 1) and a few R packages (point 2):

1. Download the new clingo (http://sourceforge.net/projects/potassco/files/clingo) from the website and rename it as "clingo"
 
You can try running:
```
./ASP/clingo
```
to see whether your clingo installation is working.


2. Install the R packages:

```
source('http://bioconductor.org/biocLite.R')
biocLite(c('graph','RBGL','gmp','RcppArmadillo'))
install.packages(c('deal','pcalg','combinat','hash','bnlearn','foreach','doMC','caTools','expm'))
```

3. Start R in the R/ directory and run:
```
source('load.R')
```
