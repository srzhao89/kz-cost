#
# Test convexity of the production set versus non-convexity using DEA
# and FDH estimators. Theoretical details are given in:
#   Kneip, A., L. Simar, and P. W. Wilson (2016), Testing Hypotheses
#        in Nonparametric Models of Production, Journal of Business and
#        Economic Statistics 34, 435--456.
# Anyone using this code should cite the above article.  The test
# implemented here uses the randomization algorithm proposed in
#   Daraio, C., L. Simar, and P. W. Wilson (2016), Nonparametric
#        Estimation of Efficiency in the Presence of Environmental
#        Variables, unpublished working paper, Department of Economics,
#        Clemson University, Clemson, South Carolina, USA
# which should also be cited.
#
# Arguments:
#      XOBS   --- np by n matrix of input vectors
#      YOBS   --- nq by n matrix of output vectors
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
# NB: If the test statistic "that" is returned as NA, this could be
#     because all of the FDH estimates in sub-sample 2 are equal to 1.
#     This is an indication of too many dimensions for the amount of 
#     data, i.e., the number of observations.  In such cases, one should
#     explore dimension-reduction possibilities.
#
#     Shirong Zhao    20 November 2023      DUFE Dalian
#
#########################################################################
#
test.cost.convexity.one.split<-function(x,y,w,NBCR=100,ISPLIT="no") {
  np=nrow(x)
  nq=nrow(y)
  kappa=1/(nq+1)
  bc.fac1=1/(2**(2/(nq+2)) - 1)
  bc.fac2=1/(2**kappa - 1)
  icase=ifelse(nq+1<=3,1,2)
  n=ncol(x)
  if (ISPLIT=="yes") {
    foo=dea.fdh.split(n,np=1,nq=nq)
    n1=foo$n1
    n2=foo$n2
  } else {
    n1=floor(n/2)
    n2=n-n1
  }
  n1a=floor(n1/2)
  n1b=n1-n1a
  n2a=floor(n2/2)
  n2b=n2-n2a
  xyw=t(rbind(x,y,w))
  n.row=nrow(xyw)
  xyw=xyw[sample(n.row,n.row,replace=FALSE),]
  rxyw=xyw
  x1=t(rxyw[1:n1,1:np])
  y1=t(rxyw[1:n1,(np+1):(np+nq)])
  w1=t(rxyw[1:n1,(np+nq+1):(np+nq+np)])
  x2=t(rxyw[(n1+1):n,1:np])
  y2=t(rxyw[(n1+1):n,(np+1):(np+nq)])
  w2=t(rxyw[(n1+1):n,(np+nq+1):(np+nq+np)])
  if (np==1) {
    x1=matrix(x1,nrow=1)
    x2=matrix(x2,nrow=1)
    w1=matrix(w1,nrow=1)
    w2=matrix(w2,nrow=1)
  }
  if (nq==1) {
    y1=matrix(y1,nrow=1)
    y2=matrix(y2,nrow=1)
  }
  # compute efficiency estimates:
  d1=1/FEAR::dea.cost.min(XOBS=x1,YOBS=y1,XPRICE=w1,RTS=1)$eff
  d2=1/fdh.cost.min.sz(XOBS=x2,YOBS=y2,XREF=x2,YREF=y2,XPRICE=w2)$eff
  # compute means, variances:
  tm1=mean(d1)
  tm2=mean(d2)
  sig1.2=var(d1)
  sig2.2=var(d2)
  # compute bias corrections via generalized jackknife:
  tbar1=rep(0,n1)
  tbar2=rep(0,n2)
  for (j in 1:NBCR) {
    if (j==1) {
      ind1=c(1:n1)
      ind2=c(1:n2)
    } else {
      ind1=sample(ind1,size=n1)
      ind2=sample(ind2,size=n2)
      x1[,1:n1]=x1[,ind1]
      y1[,1:n1]=y1[,ind1]
      w1[,1:n1]=w1[,ind1]
      x2[,1:n2]=x2[,ind2]
      y2[,1:n2]=y2[,ind2]
      w2[,1:n2]=w2[,ind2]
    }
    d1a=1/FEAR::dea.cost.min(XOBS=matrix(x1[,1:n1a],nrow=np),
                             YOBS=matrix(y1[,1:n1a],nrow=nq),
                             XPRICE=matrix(w1[,1:n1a],nrow=np),
                             RTS=1
                            )$eff
    d1b=1/FEAR::dea.cost.min(XOBS=matrix(x1[,(n1a+1):n1],nrow=np),
                             YOBS=matrix(y1[,(n1a+1):n1],nrow=nq),
                             XPRICE=matrix(w1[,(n1a+1):n1],nrow=np),
                             RTS=1
                            )$eff
    d2a=1/fdh.cost.min.sz(XOBS=matrix(x2[,1:n2a],nrow=np),
                          YOBS=matrix(y2[,1:n2a],nrow=nq),
                          XREF=matrix(x2[,1:n2a],nrow=np),
                          YREF=matrix(y2[,1:n2a],nrow=nq),
                          XPRICE=matrix(w2[,1:n2a],nrow=np)
                          )$eff
    d2b=1/fdh.cost.min.sz(XOBS=matrix(x2[,(n2a+1):n2],nrow=np),
                          YOBS=matrix(y2[,(n2a+1):n2],nrow=nq),
                          XREF=matrix(x2[,(n2a+1):n2],nrow=np),
                          YREF=matrix(y2[,(n2a+1):n2],nrow=nq),
                          XPRICE=matrix(w2[,(n2a+1):n2],nrow=np)
                          )$eff
    #
    tbar1[ind1[1:n1a]]=tbar1[ind1[1:n1a]] +
      d1a - d1[ind1[1:n1a]]
    tbar1[ind1[(n1a+1):n1]]=tbar1[ind1[(n1a+1):n1]] +
      d1b - d1[ind1[(n1a+1):n1]]
    tbar2[ind2[1:n2a]]=tbar2[ind2[1:n2a]] +
      d2a - d2[ind2[1:n2a]]
    tbar2[ind2[(n2a+1):n2]]=tbar2[ind2[(n2a+1):n2]] +
      d2b - d2[ind2[(n2a+1):n2]]
  }
  tbar1=(1/NBCR)*bc.fac1*tbar1
  tbar2=(1/NBCR)*bc.fac2*tbar2
  bc1=mean(tbar1)
  bc2=mean(tbar2)
  # compute test statistic:
  if (icase==1) {
    t1=sig1.2/n1
    t2=sig2.2/n2
    t3=t1+t2
    t4=t3**(-0.5)
    t5=tm2-tm1-(bc2-bc1)
    that=t5*t4
  } else {
    n1k=floor(n1**(2*kappa))
    n2k=floor(n2**(2*kappa))
    tmk1=mean(d1[1:n1k])
    tmk2=mean(d2[1:n2k])
    t1=sig1.2/n1k
    t2=sig2.2/n2k
    t3=t1+t2
    t4=t3**(-0.5)
    t5=tmk2-tmk1-(bc2-bc1)
    that=t5*t4
  }
  #
  # use asy. normality to get p-value:
  pval=pnorm(that,lower.tail=FALSE)
  # make a list of results to return to calling routine and then quit:
  res=list(that=that,pval=pval,muhat1=tm1,muhat2=tm2,
           bc1=bc1,bc2=bc2,v1=sig1.2,v2=sig2.2)
  return(res)
}





