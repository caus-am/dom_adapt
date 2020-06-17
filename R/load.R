#These are some of the parameters needed for operation
#

#loads the whole directory to R
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

# since the functions here are not really changed, loading them already here
loud<-function(nproc = 1) { # if -1 it tries to set to all of the available cores.

  library('deal') #used in the bayes independence test
  library('pcalg') #for running pc derivative algorithms
  library('combinat')
  library('graph') #for pcalg output handling
  library('bnlearn') #score-based learning
  library('foreach') # parallel for loops.
  library('doMC')    # parallel backend for the for loops
  #library('expm') # matrix exponential, used in cpag_to_mc
  #library('PRROC') # Plot ROC and PR curves (not a good library, better use Matlab).
  #library('caTools') # compute AUC
  if (nproc == -1){
    if (!system("nproc", ignore.stdout=TRUE, ignore.stderr = TRUE)){
      nproc <- as.numeric(system("nproc", intern=TRUE))
    } else {
      nproc <- 1
    }
  }
  registerDoMC(nproc)   
  sourceDir('./',trace=FALSE)
  sourceDir('./graph_utils/',trace=FALSE)
  sourceDir('./mods/',trace=FALSE)
  sourceDir('./simulator/',trace=FALSE)
  sourceDir('./solvers/',trace=FALSE)
  sourceDir('./tests/',trace=FALSE)
  sourceDir('./utils/',trace=FALSE)
  source('../../simulator/simul_interventions_NIPS.R')
  source('../../simulator/simulateData_NIPS.R')
}


