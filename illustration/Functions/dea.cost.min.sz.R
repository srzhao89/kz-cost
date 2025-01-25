dea.cost.min.sz<-function(XOBS=XOBS,YOBS=YOBS,XREF=XREF,YREF=YREF,XPRICE=XPRICE){
  
  n=ncol(XOBS)
  nq=nrow(YOBS)
  n1=ncol(XREF)
  ones=matrix(1,nrow=1,ncol=n)
  ones1=matrix(1,nrow=1,ncol=n1)
  eff=vector(length=n)
  for (i in 1:n) {
    cost=matrix(apply(XOBS*(XPRICE[,i]%*%ones),2,sum),nrow=1,ncol=n)
    costREF=matrix(apply(XREF*(XPRICE[,i]%*%ones1),2,sum),nrow=1,ncol=n1)
    eff[i]=FEAR::dea(XOBS=matrix(cost[1,i],nrow=1,ncol=1),
                      YOBS=matrix(YOBS[,i],nrow=nq,ncol=1),
                      XREF=costREF,
                      YREF=YREF,
                      ORIENTATION=1,
                      RTS=1,
                      METRIC=1
                     )
  }
  
  res=list(eff=eff)
  
}