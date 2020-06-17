gaussCIcontexttest <- function(x, y, S, suffStats) {
  library(ppcor)
  data<-suffStats$data
  # For each combination of values for all context variables, get one index for those values.
  # Allows us to access quickly representative rows for a given combination of context variable values.
  uniqueCVrows<-suffStats$uniqueCVrows
  p<-dim(data)[2]
  # Column indexes for the context variables.
  contextVars<-suffStats$contextVars
  # Regimes: values of the regime column.
  regimes<-suffStats$regimes
  verbose<-suffStats$verbose
  alpha<-suffStats$alpha

  if( verbose ) 
    cat( 'Testing ', x, '_||_', y, '|', S, '\n' )
  if( is.null(S) || length(S)==0 ) {
    pval <- cor.test(data[,x],data[,y])$p.value
  } else {
    # The context variable part of S:
    ScV <- intersect(S,contextVars)
    if( is.null(ScV) || length(ScV)==0 ) {
      #pval <- pcor.test(data[,x],data[,y],data[,S])$p
      data.omit <- na.omit(data[,c(x,y,S)])
      pval <- pcor.test(data.omit[,1],data.omit[,2],data.omit[,3:ncol(data.omit)])$p
    } else {
      # Regime values associated with uniqueCVrows, i.e. unique regimes (not sharing same context values)
      uniqueRegimes<-regimes[uniqueCVrows]
      
      # Subregimes: combinations of context variable values for a given subset of context variables.
      # For example if R=1, I1=0, I2=0; R=2, I1=0, I2=1; R=3, I1=1, I2=2, the subregimes w.r.t. I1 are 0 and 1.
      # This means R=1 and R=2 are indistinguishable w.r.t. I1.
      # Note: optimization of as.matrix(unique(data[,ScV])) considering only unique combinations of contexts.
      subregimes<-as.matrix(unique(data[uniqueCVrows,ScV]))
      pvals<-c()
      # counter for regimes that don't have a constant value for x or y:
      N<-0
      maxN<-dim(subregimes)[1]
      
      # For each subregime, collect the p-value of the test in this subregime in pvals:
      for( R in 1:maxN ) {
        # old and slow
#        inds_R<-1:dim(data)[1]
#        for( i in 1:length(ScV) )
#          inds_R<-intersect(inds_R,which(data[,ScV[i]]==subregimes[R,i]))
        # new and faster?
#        stopifnot(isTRUE(all.equal(inds_R,newinds_R)))

        # After the for loop this will contain the values of the regime variable for this subregime:
        setregimes<-uniqueRegimes
        for( i in 1:length(ScV) )
          setregimes<-intersect(setregimes,which(data[uniqueCVrows,ScV[i]]==subregimes[R,i]))
          
        # All indexes of the rows from datasets with regime values in setregimes:
        inds_R<-which(regimes %in% setregimes)
        if( (max(data[inds_R,x],na.rm=TRUE)>min(data[inds_R,x],na.rm=TRUE)) && (max(data[inds_R,y],na.rm=TRUE)>min(data[inds_R,y],na.rm=TRUE)) )  {
            
          # System variable part of S
          Sremain<-setdiff(S,ScV)
          if( is.null( Sremain ) || length(Sremain) == 0 ) {
            pv<-cor.test(data[inds_R,x],data[inds_R,y])$p.value
          } else {
            #pv<-pcor.test(data[inds_R,x],data[inds_R,y],data[inds_R,Sremain])$p
            data.omit <- na.omit(data[inds_R,c(x,y,Sremain)])
            if (nrow(data.omit) == 0)
              next
            pv <- pcor.test(data.omit[,1],data.omit[,2],data.omit[,3:ncol(data.omit)])$p
          }
          pvals<-c(pvals,pv)
          N<-N+1
          if( verbose )
            cat(N,'/',maxN,':',pvals[N],';',subregimes[R,],'\n')
        }
      }

      # Aggregate the p-values across the different subregimes:
      if (N==0) {
        pval<-NA
        if (verbose)
          cat('empty list->NA\n')
      } else {
        # Fisher's method (should we replace by min over p-values?)
        stat<-(-2.0) * sum(log(pvals))
        # distributed as chi^2 with 2*N d.o.f.
        pval<-pchisq(stat, df=2*N, lower.tail=FALSE)
        if( verbose )
          cat(pvals,': ',stat,'->',pval,'\n')
      }
    }
  }
  if (is.nan(pval))
    pval <- NA
  if (verbose && !is.na(pval) && pval > alpha)
    cat('conclusion: (conditional) independence', x, '_||_', y, '|', S, '\n')
  pval
}
