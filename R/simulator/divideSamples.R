divideSamples<-function( n, k ) {
  # Divides n samples to k experiments such that
  # 1) the samples are divided as evenly as possible
  # 2) all samples are used.

  # Start by taking this many experiments from each experiments.
  s<-rep(floor(n/k),k)

  while ( sum(s) < n ) { #Then add one sample to the experiment with minimum
			 #number of exp. until all samples are used.
    i<-which.min(s)
    s[i]<-s[i]+1
  }

  s
}