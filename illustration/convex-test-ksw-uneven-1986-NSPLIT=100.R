# remove all the previous stuff in R 
rm(list = ls(all=TRUE))
t1=proc.time()
require(readxl)
require(digest)
require(FEAR)
require(np)
#
source("./Functions/dea-fdh-split.R")
source("./Functions/fdh.cost.min.sz.R")
source("./Functions/dea.cost.min.sz.R")
source("./Functions/test.cost.convexity.one.split.R")
source("./Functions/test.cost.convexity.multiple.split.R")
#
if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)
######################################################
#######################  BEGIN Import Cleaned Data ##########################
df <- read.table("./Data/kt-data.txt", quote="\"", comment.char="")
######################## End Import Cleaned Data #############################
#

dim(df)
colnames(df)

##################################
res=matrix(NA, nrow=13, ncol=6)
year=1986
i=year-1986+1
df=df[df$V9==year, ]
##################################


W1=df$W1=df$V1
W2=df$W2=df$V2
W3=df$W3=df$V3
X1=df$X1=df$V4
X2=df$X2=df$V5
X3=df$X3=df$V6
Y1=df$Y1=df$V7
ID=df$ID=df$V8
Year=df$Year=df$V9
TC=df$TC=df$V10

for (Year in 1986:1998) {
  
  ii=which(df$Year==Year)
  cat(Year, ": ",length(ii), "\n")
  
}

X=rbind(X1,X2,X3)
Y=rbind(Y1)
W=rbind(W1,W2,W3)

# VRS-DEA, FEAR
cce=1/FEAR::dea.cost.min(XOBS=X,YOBS=Y,XPRICE=W,RTS=1)$eff
res[i,1]=mean(cce)
# FDH, FEAR
ncce=1/fdh.cost.min.sz(XOBS=X,YOBS=Y,XREF=X,YREF=Y,XPRICE=W)$eff
res[i,2]=mean(ncce)

#
s2=test.cost.convexity.multiple.split(x=X,y=Y,w=W,NSPLIT=100,NREP=1000,NBCR=100,ISPLIT="yes")
res[i,3]=s2$tau[1]
res[i,4]=s2$pval[1]
res[i,5]=s2$tau[2]
res[i,6]=s2$pval[2]

print(res)


