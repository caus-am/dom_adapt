source("simul_interventions_NIPS.R")
source('simulateData_NIPS.R')

# settings
NSim <- 10
dir.out <- 'data'
# dir.create(file.path('', dir.out), showWarnings = FALSE)

# test
if (0) {
  set.seed(1)
  z <- simulateData.NIPS2017()
  str(z$D[[1]]$data)
  print(z$D[[1]]$data_blinded$experiment)
  print(which(is.nan(z$D[[1]]$data_blinded$v2))) # v2 blinded; corresponds to Y under experiment 1 being blinded
}

# simply a large for loop
for (i in 1:NSim) {
  set.seed(i)
  sim <- simulateData.NIPS2017()
  #write.table(sim$D[[1]]$data, paste0(dir.out, '/sim-', i, '.csv'), row.names=FALSE, sep=",")
  write.csv(sim$D[[1]]$data, paste0('sim-', i, '.csv'), row.names=FALSE, sep=",")
  # write.table(sim$D[[1]]$data_blinded, paste0(dir.out, '/sim-hidden-', i, '.csv'), sep=",") # not necessary
}



