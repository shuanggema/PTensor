
library(glmnet)

PTensor<-function(X,y,method,threhold,r,lambda,B0){
  
  #####  INPUT #####
  # X: n*p design matrix, n--sample size  p--dimension
  # y: outcome, continuous--n*1 vector  survival--n*2 matrix with survivial time and cersoring indicator
  # method: linear or cox
  # threhold: threhold for the variance of each variable
  # r: rank of CP decomposition
  # lambda: the tuning parameter for Lasso penalty
  # B0: the initial estimation
 
  
  
  
  ### OUTPUT ####
  # fit$beta_final: the estimation of CP decompositon 
  # fit$BIC_result: the value of EBIC
  # tensor_estimation: the final estimation of the coefficients
  
  
  PenaltyMaxIter=30
  PenaltyTolFun=1e-3
  
  
  d=2
  n=dim(X)[1]
  p=dim(X)[2]
  
  
  
  if (method=='cox'){
    cut=-1e+50
  } else {
    cut=1e+50
  }
  
  beta_burnin = B0
  dev0 = Inf
  beta = beta_burnin
  
  
  option.intercept=FALSE
  penalty_factor_s=matrix(1,1,p*r);
  id=seq(from=1,by=p,to=p*r)
  penalty_factor_s[id]=0
  
  
  
  for (iter in 1:PenaltyMaxIter){
    # cyclic update of array regression coefficients
    for (j in 1:d){
      
      B_nd=beta[,,setdiff(1:d,j)]
      
      
      Xj=matrix(0,n,p*r)
      id=rowSums(B_nd!=0)!=0
      B_nd_temp=B_nd[id,]
      X_temp=X[,id]
      if (!is.matrix(X_temp)){
        X_temp=matrix(X_temp,n,1)
      }
      
    
      for (i in 1:n){
        X_dd=X[i,]%*%t(X_temp[i,])
        XX1=X_dd%*%B_nd_temp
        XX=matrix(XX1,1,p*r)
        Xj[i,]=XX
      }
     
    
      
      
      option.penalty_factor=penalty_factor_s
      
      vvv=apply(Xj,2,var)
      vvv2=abs(colMeans(Xj))
      vvv=vvv/(vvv2+1e-100)
      
      
      
      if (sum(abs(vvv)<=threhold)>0){
        
        id=which(abs(vvv)<=threhold)
        id1=which(colSums(Xj)!=0)
        
        betatmp=matrix(0,p*r,1)
        
        temp=matrix(beta[,,j],p*r,1)
        
        
        
        if (method=='ols'){
          betatmp[intersect(id,id1)]=temp[intersect(id,id1)]
          ytmp=y-Xj[,intersect(id,id1)]%*%matrix(betatmp[intersect(id,id1)],length(intersect(id,id1)),1)
        }
        
        if (sum(abs(vvv)>threhold)!=0){
          if (sum(abs(vvv)>threhold)>r){
            option.penalty_factor=option.penalty_factor[abs(vvv)>threhold]
            if (method=='ols'){
              temp=glmnet(Xj[,abs(vvv)>threhold],ytmp,family="gaussian",intercept=option.intercept,penalty.factor=option.penalty_factor,lambda=lambda)
            } else {
              temp=glmnet(Xj[,abs(vvv)>threhold],y,family='cox',penalty.factor=option.penalty_factor,lambda=lambda)
            }
            if (!is.null(temp$beta) && length(temp$lambda)==length(lambda)){
              betatmp[abs(vvv)>threhold]=temp$beta[,length(lambda)]
              beta[,,j] = matrix(betatmp,p,r)
            }
          }
        }
      } else {
     
        if (method=='ols'){
          betatmp=glmnet(Xj,y,family="gaussian",intercept=option.intercept,penalty.factor=option.penalty_factor,lambda=lambda)
          
        } else {
          betatmp=glmnet(Xj,y,family='cox',penalty.factor=option.penalty_factor,lambda=lambda)
        }
        
        
        if (length(betatmp$lambda)==length(lambda)){
          betatmp=betatmp$beta[,length(betatmp$lambda)]
          beta[,,j] = matrix(betatmp,p,r)
        }
        
      }
      
      
    }
    
    
    betatmp=matrix(beta[,,j],p*r,1)
    if (method=='ols'){
      devtmp=sum((y-Xj%*%betatmp)^2)
    } else {
      devtmp=logli_cox(Xj,y,betatmp)
    }
    
    
    #stopping rule
    diffdev = devtmp-dev0
    dev0 = devtmp
    
    
    if ((abs(diffdev)<PenaltyTolFun*(abs(dev0)+1))){
      break
    }
    
    
    
    
    #update scale of array coefficients and standardize
    
    lambda_temp=abs(beta[1,,1]*beta[1,,2])+1e-50
    
    for (j in 1:d) {
      beta[,,j]=beta[,,j]/(matrix(1,p,1)%*%abs(beta[1,,j])+1e-50)
      beta[,,j]=beta[,,j]*(matrix(1,p,1)%*%(lambda_temp)^(1/d))
    }
    
    if (devtmp>1e+50 || is.na(diffdev)){
      dev0=cut
      break
    }
  }
  beta_final = beta
  
  # output the BIC
  cutoff = 1e-8
  if (method=='ols'){
    BIC_result = log(dev0)+ log(n)*max(sum(abs(beta[,,1])>cutoff)+ sum(abs(beta[,,2])>cutoff)-r*r,0)/n
  } else {
    BIC_result = -2*dev0+ log(n)*max(sum(abs(beta[,,1])>cutoff)+ sum(abs(beta[,,2])>cutoff)-r*r,0)/n
  }
  
  
  temp=outer_product(beta_final)
  
  temp=temp+t(temp)
  temp=temp-diag(diag(temp)/2)
  
  tensor_estimation=temp
  
  tensor_estimation[abs(tensor_estimation)<=1e-2]=0
  
  fit=list(beta_final=beta_final,BIC_result=BIC_result,tensor_estimation=tensor_estimation)
  return(fit)
}



logli_cox<-function(X,y,beta){
  yy=y[y[,2]==1,1]
  y_sort=sort(yy)
  ll=matrix(0,length(yy),1)
  linear=X%*%beta
  linear_xx=linear[y[,2]==1]
  for (i in 1:length(yy)){
    ll[i]=sum(linear_xx[yy==y_sort[i]]-log(sum(exp(linear[y[,1]>=y_sort[i]]))))
  }
  logli=sum(ll)
  return(logli)
}


outer_product<-function(beta){
  b1=beta[,,1]
  b2=beta[,,2]
  
  K=dim(b1)[2]
  p=dim(b1)[1]
  
  rr=matrix(0,p,p)
  
  for (i in 1:K){
    rr=b1[,i]%*%t(b2[,i])+rr
  }
  return(rr)
}


