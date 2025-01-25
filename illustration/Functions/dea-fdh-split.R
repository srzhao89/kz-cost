#
# This function divides n into n1 nad n2 such that n1<n2 and the order
# of estimation error for FDH estimation on n2 observations is similar
# to the order of estimation error for VRS-DEA estimation on n1
# observations.  Details are given in
#   Kneip, A., L. Simar, and P. W. Wilson (2016), Testing Hypotheses
#        in Nonparametric Models of Production, Journal of Business and
#        Economic Statistics 34, 435--456.
# which should be cited by anyone using this code.
#
# Arguments:
#      n  --- number of observations
#      np --- number of inputs
#      nq --- number of outputs
# A list is returned with the following elements:
#      n1 --- number of observations for use with the VRS-DEA estimator
#      n2 --- number of observations for use with the FDH estimator
#
#     Paul W. Wilson    16 Aug2016      UQ Brisbane
#
#########################################################################
dea.fdh.split <- function(n,np,nq) {
   expon=2*(np+nq)/(np+nq+1)
   a=0
   b=0.7*n
   test=0.5*(b-a)
   tol=10**(-5)
   while(test>tol) {
      c=0.5*(a+b)
      fc=n-c-(c**expon)
      if (fc>0) {
         a=c
	 } else {
	 b=c
	 }
      test=0.5*(b-a)
      }
   n1=floor(c)
   if ((c-n1)>0.5) n1=n1+1
   n2=n-n1
   res=list(n1=n1,n2=n2)
   return(res)
   }
