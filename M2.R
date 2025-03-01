require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL)
registerDoParallel(5)
ptm=proc.time()
N=2000;d1=25;d2=500;tmax=500;K=5;M=floor(N/K)
logistic=function(x){exp(x)/(1+exp(x))}
func_link=function(x){exp(x)/(1+exp(x))}
beta_pi1_0=0
beta_pi1_S1=1*c(rep(1,d1),rep(0,d1-d1))
beta_pi2_0_0=-3
beta_pi2_0_S1=-c(0,rep(0,d1-1))
beta_pi2_0_S2=-c(0.1,rep(0,d2-1))
beta_pi2_1_0=3
beta_pi2_1_S1=c(0,rep(0,d1-1))
beta_pi2_1_S2=c(0.1,rep(0,d2-1))
beta_g0_0=1
beta_g0_S2=0.3*c(0.99^(0:(d2-1)))
beta_g0_S1=rep(0,d1)
beta_g1_0=-1
beta_g1_S2=-0.3*c(0.99^(0:(d2-1)))
beta_g1_S1=rep(0,d1)
set.seed(123)

binarys=rbinom(N*d1*tmax,1,0.5)
S1_all=matrix(binarys*runif(N*d1*tmax,-1,-0.5)+(1-binarys)*runif(N*d1*tmax,0.5,1),nrow=N*tmax)
#S1_all=matrix(runif(N*d1*tmax,-1,1),nrow=N*tmax)
delta2_all=matrix(runif(N*tmax*d2,-1,1),nrow=N*tmax)

pi1_all=func_link(beta_pi1_0+S1_all%*%beta_pi1_S1)
A1_all=rbinom(N*tmax,1,pi1_all)
hist(pi1_all)
range(pi1_all)
N/mean(1/pi1_all)

W1=matrix(0,nrow=1,ncol=d1)
W0=matrix(0,nrow=1,ncol=d1)
for (j in 1:d1){
  W1[,j]=0.5*0.9^abs(j)
  W0[,j]=0.5*0.9^abs(j)
}

S1_all_11=S1_all%*%t(W1)
S1_all_10=S1_all%*%t(W0)

W2_all=matrix(rep(S1_all_11,d2),nrow=N*tmax,ncol=d2)+delta2_all
W2_all[A1_all==0]=matrix(rep(S1_all_10,d2),nrow=N*tmax,ncol=d2)[A1_all==0]+delta2_all[A1_all==0]
W2_1_all=matrix(rep(S1_all_11,d2),nrow=N*tmax,ncol=d2)+delta2_all
W2_0_all=matrix(rep(S1_all_10,d2),nrow=N*tmax,ncol=d2)+delta2_all
S2_all=W2_all
S2_1_all=W2_1_all
S2_0_all=W2_0_all

gamma_1_0=beta_g1_0
gamma_1_S1=beta_g1_S1+sum(beta_g1_S2)*t(W1)
gamma_0_0=beta_g0_0
gamma_0_S1=beta_g0_S1+sum(beta_g0_S2)*t(W0)

mu2_1_all=beta_g1_0+S1_all%*%beta_g1_S1+W2_1_all%*%beta_g1_S2
mu2_0_all=beta_g0_0+S1_all%*%beta_g0_S1+W2_0_all%*%beta_g0_S2
mu1_1_all=gamma_1_0+S1_all%*%gamma_1_S1
mu1_0_all=gamma_0_0+S1_all%*%gamma_0_S1

theta0=beta_g0_0
theta1=beta_g1_0
theta=theta1-theta0

ww1=beta_pi2_0_0+S1_all%*%beta_pi2_0_S1+S2_0_all%*%beta_pi2_0_S2
pi2_0_all=logistic(ww1)
ww2=beta_pi2_1_0+S1_all%*%beta_pi2_1_S1+S2_1_all%*%beta_pi2_1_S2
pi2_1_all=logistic(ww2)
pi2_all=(1-A1_all)*pi2_0_all+A1_all*pi2_1_all
mean(pi2_1_all)
mean(pi2_0_all)
range(pi1_all)
range(pi2_1_all)
range(pi2_0_all)

A2_all=rbinom(N*tmax,1,pi2_all)
g_all=(A1_all+A2_all==0)*(beta_g0_0+S1_all%*%beta_g0_S1+W2_0_all%*%beta_g0_S2)+(A1_all*A2_all==1)*(beta_g1_0+S1_all%*%beta_g1_S1+W2_1_all%*%beta_g1_S2)
zeta_all=runif(N*tmax,-1,1)
Y_all=g_all+zeta_all
mean(theta1-mu2_1_all)
mean(theta1-mu1_1_all)
mean(theta0-mu2_0_all)
mean(theta0-mu1_0_all)


num_mu1est=4;num_piest=4;num_totest=num_mu1est+1
grid=data.frame(est_mu1=1:num_mu1est,est_pi=1:num_piest)

#pred_theta=c()
#for (t in 1:tmax){
results_par <- foreach (t=1:tmax, .combine=rbind, .packages=c("glmnet")) %dopar% {
  pred_pi1_1=matrix(0,nrow=num_piest,ncol=N)
  pred_pi2_0=matrix(0,nrow=num_piest,ncol=N);pred_pi2_1=matrix(0,nrow=num_piest,ncol=N)
  pred_mu2_1=matrix(0,nrow=num_mu1est,ncol=N);pred_mu1_1=matrix(0,nrow=num_mu1est,ncol=N)
  pred_mu2_0=matrix(0,nrow=num_mu1est,ncol=N);pred_mu1_0=matrix(0,nrow=num_mu1est,ncol=N)
  S1=S1_all[(t-1)*N+1:N,];S2=S2_all[(t-1)*N+1:N,]
  A1=A1_all[(t-1)*N+1:N];A2=A2_all[(t-1)*N+1:N];Y=Y_all[(t-1)*N+1:N]
  pi1=pi1_all[(t-1)*N+1:N]
  pi2_0=pi2_0_all[(t-1)*N+1:N];pi2_1=pi2_1_all[(t-1)*N+1:N]
  mu2_1=mu2_1_all[(t-1)*N+1:N];mu2_0=mu2_0_all[(t-1)*N+1:N]
  mu1_1=mu1_1_all[(t-1)*N+1:N];mu1_0=mu1_0_all[(t-1)*N+1:N]
  #oracle
  pred_pi1_1[num_piest,]=pi1
  pred_pi2_1[num_piest,]=pi2_1;pred_pi2_0[num_piest,]=pi2_0
  pred_mu2_1[num_mu1est,]=mu2_1;pred_mu2_0[num_mu1est,]=mu2_0
  pred_mu1_1[num_mu1est,]=mu1_1;pred_mu1_0[num_mu1est,]=mu1_0
  
  #estimation of pi1, pi2, mu1, mu2
  for (k in 1:K){
    index=(1:M)+(k-1)*M
    index_nk=(1:N)[-index]
    
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index1_nk1=index_nk1[which(A1[index_nk1]==1)];index0_nk1=index_nk1[which(A1[index_nk1]==0)]
    index11_nk1=index_nk1[which(A1[index_nk1]*A2[index_nk1]==1)];index00_nk1=index_nk1[which(A1[index_nk1]+A2[index_nk1]==0)]
    
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index1_nk2=index_nk2[which(A1[index_nk2]==1)];index0_nk2=index_nk2[which(A1[index_nk2]==0)]
    index11_nk2=index_nk2[which(A1[index_nk2]*A2[index_nk2]==1)];index00_nk2=index_nk2[which(A1[index_nk2]+A2[index_nk2]==0)]
    
    index1_nk=index_nk[which(A1[index_nk]==1)];index0_nk=index_nk[which(A1[index_nk]==0)]
    index11_nk=index_nk[which(A1[index_nk]*A2[index_nk]==1)];index00_nk=index_nk[which(A1[index_nk]+A2[index_nk]==0)]
    #pi1
    fit_Log1=cv.glmnet(x=S1[index_nk,],y=A1[index_nk],family="binomial")
    pred_pi1_1[1,index]=predict(fit_Log1,newx=S1[index,],type="response",s="lambda.min") #logistic+lasso
    pred_pi1_1[2,index]=pred_pi1_1[1,index]
    pred_pi1_1[3,index]=pred_pi1_1[1,index]
    #pi2_1
    temp=try({fit_Log2_1=cv.glmnet(x=cbind(S1[index1_nk,],S2[index1_nk,]),y=A2[index1_nk],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_1[1,index]=rep(mean(A2[index1_nk]),length(index))
    } else {
      pred_pi2_1[1,index]=predict(fit_Log2_1,newx=cbind(S1[index,],S2[index,]),type="response",s="lambda.min") #logistic
    }
    pred_pi2_1[2,index]=pred_pi2_1[1,index]
    pred_pi2_1[3,index]=pred_pi2_1[1,index]
    #pi2_0
    temp=try({fit_Log2_0=cv.glmnet(x=cbind(S1[index0_nk,],S2[index0_nk,]),y=A2[index0_nk],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_0[1,index]=rep(mean(A2[index0_nk]),length(index))
    } else {
      pred_pi2_0[1,index]=predict(fit_Log2_0,newx=cbind(S1[index,],S2[index,]),type="response",s="lambda.min") #logistic
    }
    pred_pi2_0[2,index]=pred_pi2_0[1,index]
    pred_pi2_0[3,index]=pred_pi2_0[1,index]
    #mu2_1
    fit_lasso1=cv.glmnet(x=cbind(S1[index11_nk,],S2[index11_nk,]),y=Y[index11_nk],family="gaussian",alpha=1)
    pred_mu2_1[1,index]=predict(fit_lasso1,newx=cbind(S1[index,],S2[index,]),s="lambda.min") #lasso
    pred_mu2_1[2,index]=pred_mu2_1[1,index]
    pred_mu2_1[3,index]=pred_mu2_1[1,index]
    #mu2_0
    fit_lasso0=cv.glmnet(x=cbind(S1[index00_nk,],S2[index00_nk,]),y=Y[index00_nk],family="gaussian",alpha=1)
    pred_mu2_0[1,index]=predict(fit_lasso0,newx=cbind(S1[index,],S2[index,]),s="lambda.min") #lasso
    pred_mu2_0[2,index]=pred_mu2_0[1,index]
    pred_mu2_0[3,index]=pred_mu2_0[1,index]
    #DTL
    #mu1_1
    pred_mu2_1_nocross=predict(fit_lasso1,newx=cbind(S1[index1_nk,],S2[index1_nk,]),s="lambda.min")
    fit_lasso2_nocross=cv.glmnet(x=S1[index1_nk,],y=pred_mu2_1_nocross,family="gaussian",alpha=1)
    pred_mu1_1[2,index]=predict(fit_lasso2_nocross,newx=S1[index,],s="lambda.min")
    #mu1_0
    pred_mu2_0_nocross=predict(fit_lasso0,newx=cbind(S1[index0_nk,],S2[index0_nk,]),s="lambda.min")
    fit_lasso2_nocross=cv.glmnet(x=S1[index0_nk,],y=pred_mu2_0_nocross)
    pred_mu1_0[2,index]=predict(fit_lasso2_nocross,newx=S1[index,],s="lambda.min")
    #D-DRL
    #mu1_1
    fit_lasso1_new1=cv.glmnet(x=cbind(S1[index11_nk1,],S2[index11_nk1,]),y=Y[index11_nk1],family="gaussian",alpha=1)
    temp=try({fit_Log2_1_new1=cv.glmnet(x=cbind(S1[index1_nk1,],S2[index1_nk1,]),y=A2[index1_nk1],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_1_new1=rep(mean(A2[index1_nk1]),length(index1_nk2))
    } else {
      pred_pi2_1_new1=predict(fit_Log2_1_new1,newx=cbind(S1[index1_nk2,],S2[index1_nk2,]),type="response",s="lambda.min")
    }
    pred_mu2_1_new1=predict(fit_lasso1_new1,newx=cbind(S1[index1_nk2,],S2[index1_nk2,]),s="lambda.min")
    Y1_new1=as.numeric(pred_mu2_1_new1+A2[index1_nk2]*(Y[index1_nk2]-pred_mu2_1_new1)/pred_pi2_1_new1)
    fit_lasso2_nocross1_new1=cv.glmnet(x=S1[index1_nk2,],y=Y1_new1,family="gaussian",alpha=1)
    pred_mu1_1_new1=predict(fit_lasso2_nocross1_new1,newx=S1[index,],s="lambda.min")
    
    fit_lasso1_new2=cv.glmnet(x=cbind(S1[index11_nk2,],S2[index11_nk2,]),y=Y[index11_nk2],family="gaussian",alpha=1)
    temp=try({fit_Log2_1_new2=cv.glmnet(x=cbind(S1[index1_nk2,],S2[index1_nk2,]),y=A2[index1_nk2],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_1_new2=rep(mean(A2[index1_nk2]),length(index1_nk1))
    } else {
      pred_pi2_1_new2=predict(fit_Log2_1_new2,newx=cbind(S1[index1_nk1,],S2[index1_nk1,]),type="response",s="lambda.min")
    }
    pred_mu2_1_new2=predict(fit_lasso1_new2,newx=cbind(S1[index1_nk1,],S2[index1_nk1,]),s="lambda.min")
    Y1_new2=as.numeric(pred_mu2_1_new2+A2[index1_nk1]*(Y[index1_nk1]-pred_mu2_1_new2)/pred_pi2_1_new2)
    fit_lasso2_nocross1_new2=cv.glmnet(x=S1[index1_nk1,],y=Y1_new2,family="gaussian",alpha=1)
    pred_mu1_1_new2=predict(fit_lasso2_nocross1_new2,newx=S1[index,],s="lambda.min")
    pred_mu1_1[1,index]=(pred_mu1_1_new1+pred_mu1_1_new2)/2
    
    #mu1_0
    fit_lasso0_new1=cv.glmnet(x=cbind(S1[index00_nk1,],S2[index00_nk1,]),y=Y[index00_nk1],family="gaussian",alpha=1)
    temp=try({fit_Log2_0_new1=cv.glmnet(x=cbind(S1[index0_nk1,],S2[index0_nk1,]),y=A2[index0_nk1],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_0_new1=rep(mean(A2[index0_nk1]),length(index0_nk2))
    } else {
      pred_pi2_0_new1=predict(fit_Log2_0_new1,newx=cbind(S1[index0_nk2,],S2[index0_nk2,]),type="response",s="lambda.min")
    }
    pred_mu2_0_new1=predict(fit_lasso0_new1,newx=cbind(S1[index0_nk2,],S2[index0_nk2,]),s="lambda.min")
    Y0_new1=as.numeric(pred_mu2_0_new1+(1-A2[index0_nk2])*(Y[index0_nk2]-pred_mu2_0_new1)/(1-pred_pi2_0_new1))
    fit_lasso2_nocross0_new1=cv.glmnet(x=S1[index0_nk2,],y=Y0_new1,family="gaussian",alpha=1)
    pred_mu1_0_new1=predict(fit_lasso2_nocross0_new1,newx=S1[index,],s="lambda.min")
    
    fit_lasso0_new2=cv.glmnet(x=cbind(S1[index00_nk2,],S2[index00_nk2,]),y=Y[index00_nk2],family="gaussian",alpha=1)
    temp=try({fit_Log2_0_new2=cv.glmnet(x=cbind(S1[index0_nk2,],S2[index0_nk2,]),y=A2[index0_nk2],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_0_new2=rep(mean(A2[index0_nk2]),length(index0_nk1))
    } else {
      pred_pi2_0_new2=predict(fit_Log2_0_new2,newx=cbind(S1[index0_nk1,],S2[index0_nk1,]),type="response",s="lambda.min")
    }
    pred_mu2_0_new2=predict(fit_lasso0_new2,newx=cbind(S1[index0_nk1,],S2[index0_nk1,]),s="lambda.min")
    Y0_new2=as.numeric(pred_mu2_0_new2+(1-A2[index0_nk1])*(Y[index0_nk1]-pred_mu2_0_new2)/(1-pred_pi2_0_new2))
    fit_lasso2_nocross0_new2=cv.glmnet(x=S1[index0_nk1,],y=Y0_new2,family="gaussian",alpha=1)
    pred_mu1_0_new2=predict(fit_lasso2_nocross0_new2,newx=S1[index,],s="lambda.min")
    pred_mu1_0[1,index]=(pred_mu1_0_new1+pred_mu1_0_new2)/2
    #D-DRL'
    #mu1_1
    fit_lasso1_new=cv.glmnet(x=cbind(S1[index11_nk,],S2[index11_nk,]),y=Y[index11_nk],family="gaussian",alpha=1)
    temp=try({fit_Log2_1_new=cv.glmnet(x=cbind(S1[index1_nk,],S2[index1_nk,]),y=A2[index1_nk],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_1_new=rep(mean(A2[index1_nk]),length(index1_nk))
    } else {
      pred_pi2_1_new=predict(fit_Log2_1_new,newx=cbind(S1[index1_nk,],S2[index1_nk,]),type="response",s="lambda.min")
    }
    pred_mu2_1_new=predict(fit_lasso1_new,newx=cbind(S1[index1_nk,],S2[index1_nk,]),s="lambda.min")
    Y1_new=as.numeric(pred_mu2_1_new+A2[index1_nk]*(Y[index1_nk]-pred_mu2_1_new)/pred_pi2_1_new)
    fit_lasso2_nocross1_new=cv.glmnet(x=S1[index1_nk,],y=Y1_new,family="gaussian",alpha=1)
    pred_mu1_1[3,index]=predict(fit_lasso2_nocross1_new,newx=S1[index,],s="lambda.min")
    
    #mu1_0
    fit_lasso0_new=cv.glmnet(x=cbind(S1[index00_nk,],S2[index00_nk,]),y=Y[index00_nk],family="gaussian",alpha=1)
    temp=try({fit_Log2_0_new=cv.glmnet(x=cbind(S1[index0_nk,],S2[index0_nk,]),y=A2[index0_nk],family="binomial")},TRUE)
    if ('try-error' %in% class(temp)){
      pred_pi2_0_new=rep(mean(A2[index0_nk]),length(index0_nk))
    } else {
      pred_pi2_0_new=predict(fit_Log2_0_new,newx=cbind(S1[index0_nk,],S2[index0_nk,]),type="response",s="lambda.min")
    }
    pred_mu2_0_new=predict(fit_lasso0_new,newx=cbind(S1[index0_nk,],S2[index0_nk,]),s="lambda.min")
    Y0_new=as.numeric(pred_mu2_0_new+(1-A2[index0_nk])*(Y[index0_nk]-pred_mu2_0_new)/(1-pred_pi2_0_new))
    fit_lasso2_nocross0_new=cv.glmnet(x=S1[index0_nk,],y=Y0_new,family="gaussian",alpha=1)
    pred_mu1_0[3,index]=predict(fit_lasso2_nocross0_new,newx=S1[index,],s="lambda.min")
  }
  
  #MSE check
  mse_mu2_0=rowMeans((pred_mu2_0[,which(A1==0)]-matrix(rep(mu2_0[which(A1==0)],each=num_mu1est),nrow=num_mu1est))^2)
  mse_mu1_0=rowMeans((pred_mu1_0[,which(A1==0)]-matrix(rep(mu1_0[which(A1==0)],each=num_mu1est),nrow=num_mu1est))^2)
  mse_mu2_1=rowMeans((pred_mu2_1[,which(A1==1)]-matrix(rep(mu2_1[which(A1==1)],each=num_mu1est),nrow=num_mu1est))^2)
  mse_mu1_1=rowMeans((pred_mu1_1[,which(A1==1)]-matrix(rep(mu1_1[which(A1==1)],each=num_mu1est),nrow=num_mu1est))^2)
  mse_pi1=rowMeans((pred_pi1_1-matrix(rep(pi1,each=num_piest),ncol=N))^2)
  mse_pi2_0=rowMeans((pred_pi2_0[,which(A1==0)]-matrix(rep(pi2_0[which(A1==0)],each=num_piest),nrow=num_piest))^2)
  mse_pi2_1=rowMeans((pred_pi2_1[,which(A1==1)]-matrix(rep(pi2_1[which(A1==1)],each=num_piest),nrow=num_piest))^2)
  
  #estimation of theta
  pred_theta_t=c();pred_sig2_t=c()
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    gamma2_1=A1*A2/(pred_pi1_1[x[2],]*pred_pi2_1[x[2],]);gamma1_1=A1/pred_pi1_1[x[2],]
    psi_1=gamma2_1*Y-(gamma1_1-1)*pred_mu1_1[x[1],]-(gamma2_1-gamma1_1)*pred_mu2_1[x[1],]
    gamma2_0=(1-A1)*(1-A2)/((1-pred_pi1_1[x[2],])*(1-pred_pi2_0[x[2],]));gamma1_0=(1-A1)/(1-pred_pi1_1[x[2],])
    psi_0=gamma2_0*Y-(gamma1_0-1)*pred_mu1_0[x[1],]-(gamma2_0-gamma1_0)*pred_mu2_0[x[1],]
    pred_theta_t=c(pred_theta_t,mean(psi_1-psi_0))
  }
  #asymptotic variance estimator
  for (j in 1:nrow(grid)){
    x=as.numeric(grid[j,])
    gamma2_1=A1*A2/(pred_pi1_1[x[2],]*pred_pi2_1[x[2],]);gamma1_1=A1/pred_pi1_1[x[2],]
    psi_1=gamma2_1*Y-(gamma1_1-1)*pred_mu1_1[x[1],]-(gamma2_1-gamma1_1)*pred_mu2_1[x[1],]
    gamma2_0=(1-A1)*(1-A2)/((1-pred_pi1_1[x[2],])*(1-pred_pi2_0[x[2],]));gamma1_0=(1-A1)/(1-pred_pi1_1[x[2],])
    psi_0=gamma2_0*Y-(gamma1_0-1)*pred_mu1_0[x[1],]-(gamma2_0-gamma1_0)*pred_mu2_0[x[1],]
    psi=psi_1-psi_0-pred_theta_t[j]
    pred_sig2_t=c(pred_sig2_t,mean(psi^2))
  }
  pred_theta_t=c(mean(Y[which(A1*A2==1)])-mean(Y[which(A1+A2==0)]),pred_theta_t)
  pred_sig2_t=c(var(Y[which(A1*A2==1)])*sum(A1*A2==1)/N+var(Y[which(A1+A2==0)])*sum(A1+A2==0)/N,pred_sig2_t)
  c(pred_theta=pred_theta_t,pred_sig2=pred_sig2_t,mse_mu2_0,mse_mu1_0,mse_mu2_1,mse_mu1_1,mse_pi1,mse_pi2_0,mse_pi2_1)
  #pred_theta=cbind(pred_theta,c(pred_theta=pred_theta_t,pred_sig2=pred_sig2_t,mse_mu2_0,mse_mu1_0,mse_mu2_1,mse_mu1_1,mse_pi1,mse_pi2_0,mse_pi2_1))
}
pred_theta=results_par[,1:num_totest]
pred_sig2=results_par[,(num_totest+1):(2*num_totest)];
mse_mu2_0=apply(results_par[,(2*num_totest+1):(2*num_totest+num_mu1est)],2,median)
mse_mu1_0=apply(results_par[,(2*num_totest+num_mu1est+1):(2*num_totest+num_mu1est+num_mu1est)],2,median)
mse_mu2_1=apply(results_par[,(2*num_totest+num_mu1est+num_mu1est+1):(2*num_totest+2*num_mu1est+num_mu1est)],2,median)
mse_mu1_1=apply(results_par[,(2*num_totest+2*num_mu1est+num_mu1est+1):(2*num_totest+2*num_mu1est+2*num_mu1est)],2,median)
mse_pi1=apply(results_par[,(2*num_totest+2*num_mu1est+2*num_mu1est+1):(2*num_totest+2*num_mu1est+2*num_mu1est+num_piest)],2,median)
mse_pi2_0=apply(results_par[,(2*num_totest+2*num_mu1est+2*num_mu1est+num_piest+1):(2*num_totest+2*num_mu1est+2*num_mu1est+2*num_piest)],2,median)
mse_pi2_1=apply(results_par[,(2*num_totest+2*num_mu1est+2*num_mu1est+2*num_piest+1):(2*num_totest+2*num_mu1est+2*num_mu1est+3*num_piest)],2,median)
bias=apply(pred_theta-theta,2,median)
RMSE=apply((pred_theta-theta)^2,2,median)^0.5
AL=apply(2*sqrt(pred_sig2/N)*qnorm(1-0.05/2),2,median)
AC=apply(abs(pred_theta-theta)<=sqrt(pred_sig2/N)*qnorm(1-0.05/2),2,mean)
ESD=apply(pred_theta,2,function(x){1.4826*median(abs(x-median(x)))})
ASD=apply(pred_sig2^0.5,2,median)/sqrt(N)
est_pi=factor(grid$est_pi,levels=1:num_piest,labels=c("log-DDRL","log-DTL","log-DDRL'","oracle"))
est_mu1=factor(grid$est_mu1,levels=1:num_mu1est,labels=c("DDRL","DTL","DDRL'","oracle"))
est=c("empdiff",apply(t(1:nrow(grid)),2,function(x){paste(est_pi[x],est_mu1[x],sep=" + ")}))
table_result=function(x){
  paste(est[x],"&",format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(AL[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(AC[x],digits=3),nsmall=3,scientific=FALSE),"&",
        format(round(ESD[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(ASD[x],digits=3),nsmall=3,scientific=FALSE),sep="")
}
result=apply(t(1:num_totest),2,table_result)
result=list(result=c("estimator&Bias&RMSE&AL&AC&ESD&ASD",result),est_pi=levels(est_pi),mse_pi1=paste((mse_pi1)),
            mse_pi2_0=paste((mse_pi2_0)),mse_pi2_1=paste((mse_pi2_1)),est_mu1=levels(est_mu1),
            mse_mu2_0=paste((mse_mu2_0)),mse_mu1_0=paste((mse_mu1_0)),
            mse_mu2_1=paste((mse_mu2_1)),mse_mu1_1=paste((mse_mu1_1)))
result
nrow(pred_theta)
proc.time()-ptm
