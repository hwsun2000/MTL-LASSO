
#MTL-LASSO for binary penalized logistic regression
#Title:Robust High Dimensional Logistic Regression through Maximum Trimmed Likelihood  Estimator

#Date:2018-6-17
#Author:Hongwei Sun


#variable list
#x             a numeric matrix containing the predictor variables
#y             a numeric vector containing the binary response variable.
#h             a numeric value giving the percentage of the residuals for which the L1  penalized sum of squares should be minimized (the default is 0.75). 
#ncsteps        a positive integer giving the number of C-steps to perform on all  subsamples inthe first phase of the algorithm (the default is to perform two C-steps)
#formersubset   a numeric value giving the number of subsamples in the first phase of  the algorithm. The default is to first perform nstep C-steps on 500 intitial subsamples,  and then to keep the 10 subsamples with the lowest 
#latersubset    a numeric value giving the number of subsamples in the second phase of  the algorithm.The default is 10 which means keep 10 subsamples with lowest value of the  objective function after performing nstep C-steps on 500 initial samples.

library(glmnet)

lts.logistic<-function(x,y,h=0.95,ncsteps=2,formersubset=600,latersubset=10 ){

  n<-length(y)  # sample size
  p<-dim(x)[2]  # number of variables,
   
  hn<-floor(n*h) 
 indices1=which(y==1) 
  indices0=which(y==0)
  n1=length(indices1)
  n0=length(indices0)

  
  #robustfied point-biserial correlation  
   
  rob_r=rep(0,p)   
  for (j in 1:p)
   {
   rob_r[j]=sqrt(n0*n1/(n*(n-1)))*(median(x[indices1,j])-median(x[indices0,j]))/mad(x [,j])
   }                                 
 
  maxcc=max(abs(rob_r))   # the largest value of lambda


  frac=seq(0.05,1,by=0.05)*maxcc  # lambda sequence for C-step.
      
      center=apply(x,2,median)
     x.c = sweep(x, 2, center)
     mad.c=apply(x,2,mad)
     xsd=sweep(x.c,2,mad.c,"/")

 
  #cstep function

   cstep<-function(indices){
 
     for (i in 1:ncsteps){
      
     #normalize by trimmed mean and standard
     center = colMeans(xsd[indices,])
     x.c = sweep(xsd[indices,], 2, center)
     sd.c=apply(xsd[indices,],2,sd)
     xsd_trim=sweep(x.c,2,sd.c,"/")
 

     LASSO=cv.glmnet(xsd_trim,y[indices],standardize=FALSE,lambda=frac,alpha=1,family="binomial")
     lpredict= predict(LASSO,xsd[1:n,],type="response")  #predictive probablility
     ww=coef(LASSO)       #estimated coefficients
     mybetas=as.matrix(ww)
     res=y*log(lpredict)+(1-y)*log((1-lpredict))   #likelihood function
     o<-order(res,decreasing=TRUE)
     indices<-sort(o[1:hn])    # reserve  hn samples with the largest hn  likelihood  function

     if(sum(y[indices])<=1|sum(y[indices])>=length(y[indices])-1) # if there are less  than 2 samples in either categories of y,then break the loop.
             {break}    
 }    
     crit=sum(res[indices])
     result<-list(indices=indices,crit=crit,mybetas=mybetas)
     return(result)
    }


  # initial subset 
   indices1=which(y==1)
   indices0=which(y==0)
   subset2=matrix(0,6,500)

   for (j in 1:500)
  {
   k1=sample(indices1,3)
   k0=sample(indices0,3)
   subset2[,j]=t(c(k1,k0))
   }    
      
   coef=matrix(0,hn,ncol(subset2))
   output=rep(0,ncol(subset2))    
   xishu=matrix(0,p+1,ncol(subset2))
 
# perform nsctep C-step on formersubset samples
   for (k in 1:ncol(subset2))  
   {
    jieguo<-cstep(subset2[,k])
    output[k]<-jieguo$crit    #likelihood function
    coef[,k]<-jieguo$indices  #subscript of subsamples
    xishu[,k]<-jieguo$mybetas  # coefficient 
    }


   ##remove the subsets that there are less than 2 samples in either categories of y

   series2=rep(0,ncol(subset2))   
   coef2=matrix(0,hn,ncol(subset2))  
   output2=rep(0,ncol(subset2))    
   j=1
   for(k in 1:ncol(subset2)) 
   {
    series2[k]=sum(y[coef[,k]])
    if (series2[k]>1&series2[k]<hn-1)
    {output2[j]=output[k]
     coef2[,j]=coef[,k]
     j=j+1
     }
    }
   output2=output2[-(j:ncol(subset2))]
   coef2=coef2[,-(j:ncol(subset2))]

   ##keep latersubset subsets with largest value of the likelihood function
   betterindices=order(output2,decreasing=TRUE)[1:latersubset]
   a=coef2[,betterindices]

   ##The second C-step function which perform C-step until convergence

   crit=0       
   continueCstep=TRUE
   object=rep(0,50)
   j=0 

   cstep2<-function(indices){
    while (continueCstep){
      previousCrit = crit
  
      center = colMeans(xsd[indices,])
      x.c = sweep(xsd[indices,], 2, center)
      sd.c=apply(xsd[indices,],2,sd)
      xsd_trim=sweep(x.c,2,sd.c,"/")

      LASSO=cv.glmnet(xsd_trim,y[indices],standardize=FALSE,lambda=frac,alpha=1,family="binomial")
      lpredict= predict(LASSO,xsd[1:n,],type="response")
      ww=coef(LASSO)       #the estimated coeficients
      mybetas=as.matrix(ww)
      res=y*log(lpredict)+(1-y)*log((1-lpredict))
      o<-order(res,decreasing=TRUE)
      indices<-sort(o[1:hn])
      crit=sum(res[indices])
      j=j+1
      if(sum(y[indices])<=1|sum(y[indices])>=length(y[indices])-1)
      {break} 
      continueCstep=(abs(previousCrit - crit) > 0.001)
     }
      result<-list(indices=indices, crit=crit,mybetas=mybetas,object=object)
      return(result)
   }
      output3=rep(0,latersubset)
      coef3=matrix(0,hn,latersubset)
      mybetas3=matrix(0,p+1,latersubset)


   ## for latersubset subsets, perform cstep2 function

   for (k in 1:latersubset)  
   {
   jieguo2<-cstep2(a[,k])
   output3[k]<-jieguo2$crit     #likelihood function
   coef3[,k]<-jieguo2$indices   #subscript of subsamples
   mybetas3[,k]<-jieguo2$mybetas # coefficient 
   }

    ##remove the subsets that there are less than 2 samples in either categories of y
    series3=rep(0,latersubset)   
    coef4=matrix(0,hn,latersubset) 
    output4=rep(0,latersubset)    
    j=1
    for(k in 1:latersubset) 
    {
     series3[k]=sum(y[coef3[,k]])
     if (series3[k]>1&series2[k]<hn-1)
     {output4[j]=output3[k]
      coef4[,j]=coef3[,k]
      j=j+1
      }
     }
    output4=output4[-(j:latersubset)]
    coef4=coef4[,-(j:latersubset)]

    #keep the optimal subset with largest likelihood
    lastindices=which.max(output4)
    lastcoef=coef4[,lastindices]
    lastcrit=output4[lastindices]
               
    # perform LASSO on the optimal subset
   center = colMeans(xsd[lastcoef,])
   x.c = sweep(xsd[lastcoef,], 2, center)
   sd.c=apply(xsd[lastcoef,],2,sd)
   xsd_trim=sweep(x.c,2,sd.c,"/")

   LASSO=cv.glmnet(xsd_trim,y[lastcoef],standardize=FALSE,alpha=1,family="binomial")
   raw.betas=coef(LASSO)    
   rawpredict= predict(LASSO,xsd[1:n,],type="response")

   # reweighted step
   perres=(y-rawpredict)/sqrt(rawpredict*(1-rawpredict))
   q <- qnorm(0.9875) 
   ok<- as.integer(abs(perres) <= q)
   rwt.indices<-sort(which(ok>0))

   # perform LASSO on the reweighted subset
   if (sum(y[rwt.indices])>1&sum(y[rwt.indices])<length(y[rwt.indices])-1)  
    {
      center = colMeans(xsd[rwt.indices,])
      x.c = sweep(xsd[rwt.indices,], 2, center)
      sd.c=apply(xsd[rwt.indices,],2,sd)
      xsd_trim=sweep(x.c,2,sd.c,"/")
     rwt.lasso<-cv.glmnet(xsd_trim,y[rwt.indices],standardize=FALSE,alpha=1,family="binomial")
     rwt.betas=coef(rwt.lasso)      
     lpredict= predict(rwt.lasso,newx=xsd[1:n,],type="response")
    }  else
    {  rwt.betas=raw.betas
       lpredict=rawpredict
    }
 
   result<-list(perres=perres,raw.betas=raw.betas,rawpredict=rawpredict,ok=ok,rwt.indices=rwt.indices,rwt.betas=rwt.betas,lpredict=lpredict)
   return(result)
    }





