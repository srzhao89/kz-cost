#   L. Simar, and P. W. Wilson (2019)
# Hypothesis Testing in Nonparametric Models of Production
# using Multiple Sample Splits
#
# Arguments:
#      XOBS   --- np by n matrix of input vectors
#      YOBS   --- nq by n matrix of output vectors
#      NSPLIT  --- Numer of Splits
#      NREP   --- Number of bootstraps
#      ORIENTATION --- equals 1 for input orientation, 2 for output
#                      orientation, or 3 for hyperbolic orientation
#      ISPLIT --- equals "yes" for an uneven split of the sample so that
#                 VRS-DEA and FDH estimates have same order of
#                 estimation error, or "no" for an even split.
# Kneip et al. (2016) used uneven sample splits in their simulations,
# but some care should be taken.  In particular, if the dimensionality
# is too large for a given sample size, uneven splitting may give
# too-few observations to the sub-sample on which VRS-DEA is used to
# yield meaningful estimates.  Daraio et al. (2016) use even splits in
# their separability test to avoid such problems.  The default value 
# is "no".
#
#
#     Shirong Zhao    20 November 2023      DUFE Dalian
#
#########################################################################

test.cost.convexity.multiple.split<-function(x,y,w,NSPLIT=100,NREP=1000,NBCR=100,ISPLIT="no") {
  n=ncol(x)
  # Compute the statistics
  res1=rep(-99,NSPLIT)
  res2=rep(-99,NSPLIT)
  for (iper in 1:NSPLIT) {
    t1=test.cost.convexity.one.split(x,y,w,NBCR=NBCR,ISPLIT=ISPLIT) 
    res1[iper]=t1$that
    res2[iper]=t1$pval
  }
  T.hat=mean(res1)
  K.hat.i=function(pvalue){
    P=ecdf(pvalue) 
    diff=function(x){abs(P(x)-x)}
    out=optimize(diff,interval=c(0,1),maximum=TRUE)
    return(out$objective)
  }
  K.hat=K.hat.i(res2)
  ### Do the bootstrap
  res1.boot=rep(-99,NSPLIT)
  res2.boot=rep(-99,NSPLIT)
  T.hat.boot=rep(-99,NREP)
  K.hat.boot=rep(-99,NREP)
  for (i in 1:NREP) {
    m=sample(n,n,replace=TRUE)
    # Bootstrap the sample
    d1=1/FEAR::dea(x,y,RTS=1,ORIENTATION=1) # d1<=1
    d1.boot=d1[m]
    x.boot=t(t(x)*d1/d1.boot)
    y.boot=y
    #
    for (iper in 1:NSPLIT) {
      t1=test.cost.convexity.one.split(x.boot,y.boot,w,NBCR=NBCR,ISPLIT=ISPLIT) 
      res1.boot[iper]=t1$that
      res2.boot[iper]=t1$pval
    }
    T.hat.boot[i]=mean(res1.boot)
    K.hat.boot[i]=K.hat.i(res2.boot)
    if (i %% 100 ==0) {
    cat("Bootstrap",i,"\n")
    }
  }
  # Calculate the corresponding p-value
  p.value.T=mean(T.hat.boot>=T.hat)
  p.value.K=mean(K.hat.boot>=K.hat)
  # make a list of results to return to calling routine and then quit:
  res=list(tau=c(T.hat,K.hat),pval=c(p.value.T,p.value.K))
  return(res)
}
