# dom_adapt
Code for the NIPS 2018 paper "Domain Adaptation by Using Causal Inference to Predict Invariant Conditional Distributions"



# Installation instructions:

1. Check installation for the existence of ./ASP/clingo
If it’s not there download the new clingo (http://sourceforge.net/projects/potassco/files/clingo) from the website and rename it to that file.
 
You can try running:
```
./ASP/clingo
```
to see whether your clingo installation is working.


2. Install R packages:

```
source('http://bioconductor.org/biocLite.R')
biocLite(c('graph','RBGL','gmp','RcppArmadillo'))
install.packages(c('deal','pcalg','combinat','hash','bnlearn','foreach','doMC','caTools','expm'))
```
