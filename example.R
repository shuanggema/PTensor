library(MASS)


source("PTensor_function.R")


### generate dataset ###

p=1000
n=150

rho=0.3
sigma=matrix(0,p,p)
for (i in 1:p){
  for (j in 1:p){
    sigma[i,j]=rho^(abs(i-j))
  }
  }



  
b=matrix(0,p+1,p+1)
bb=runif(16,0.6,1)
b[1:16,1]=bb
b[1,1:16]=bb
bb2=matrix(runif(36,0.6,1),6,6)
b[1:6,1:6]=(bb2+t(bb2))/2
diag(b)=0
b[1,1]=runif(1,0.6,1)

b_vector=matrix(0,1+p+p+p*(p-1)/2,1);
bb_ll=c(0, seq(from=(p+1),to=1,by=-1))
for (i in 1:(p+1)){
  b_vector[(sum(bb_ll[1:i])+1):sum(bb_ll[1:(i+1)])]=b[i,i:(p+1)]
  }
  
  
  
  
pp= dim(b)
p1=pp[1]
p2=pp[2]


  
X=mvrnorm(n,matrix(0,p,1),sigma)
  

X=cbind(matrix(1,n,1), X)
  
pp_vector=1+p+p+p*(p-1)/2
  
M_vector=matrix(0,n,1+p+p+p*(p-1)/2)
for (i in 1:n){
  each_xx=matrix(X[i,],p+1,1)%*%X[i,]
  for (j in 1:(p+1)){
    M_vector[i,(sum(bb_ll[1:j])+1):sum(bb_ll[1:(j+1)])]=each_xx[j,j:(p+1)]
    }
  }
  
  
sigma_e = 1;  # noise level
e=sigma_e*rnorm(n,0,1)
y=as.numeric(M_vector%*%b_vector+e)
  
  
######### PTensor method  

method='ols'
threhold=0.01

# The rank of the CP decomposition, usually less than ten
r=5

# The tuning parameter of the Lasso penalty
lambda=0.4

# 
X=scale(X)
X[,1]=matrix(1,n,1)


# random initialization of the estimation 
d=2
pp=dim(X)[2]
B0=array(runif(pp*r*d,-1,1),dim=c(pp,r,d)) 


fit=PTensor(X,y,method,threhold,r,lambda,B0)
coef_estimation=fit$tensor_estimation
