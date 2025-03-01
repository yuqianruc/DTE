require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL);library(ggplot2)
ptm=proc.time()
registerDoParallel(4)
library(readr)
library(haven)
dat <- read_sas("~/Downloads/bencost.sas7bdat")
dat[,34]
data <- read_csv("~/Downloads/data.csv")
dim(data)
#View(data)
tmp=cbind(dat[,34],data[,c(2:2340)])
tmp=tmp[complete.cases(tmp),]
tmp=as.data.frame(tmp)
dim(tmp)
#data1=tmp[tmp[,2] !=-1 & tmp[,2] !=0 & tmp[,2] !=2 & tmp[,3] !=-1,];dim(data1);set.seed(123)
data1=tmp[tmp[,2] !=-1 & tmp[,3] !=-1,];dim(data1);set.seed(123)

s1_all=c()
for (i in 5:913){
  s1_all=cbind(s1_all,data1[,i])
}
S1_all=as.matrix(s1_all);dim(S1_all)
A1_all=data1[,2]
s2_all=c()
for (i in 914:2340){
  s2_all=cbind(s2_all,data1[,i])
}
S2_all=as.matrix(s2_all);dim(S2_all)
S2_1_all=S2_all[A1_all==3,];dim(S2_1_all)
S2_0_all=S2_all[A1_all==2,];dim(S2_0_all)
A2_all=data1[,3]
Y_all=data1[,1]


indicator1 <- function(data) I(data == 2) + 1 - 1
indicator2 <- function(data) I(data == 1) + 1 - 1
indicator3 <- function(data) I(data == 1) + 1 - 1

K=5;N=length(A1_all);M=floor(N/K)
S1=S1_all;S2=S2_all
A1_3=indicator1(A1_all)
A2_3=indicator1(A2_all)
A1_1=indicator2(A1_all)
A2_1=indicator2(A2_all)
Y=log(Y_all+1)

mean(A1_3*A2_3)*8570
mean(A1_1*A2_1)*8570

num_mu1est=1;num_piest=1;num_totest=num_mu1est+1
pred_pi1_3=matrix(0,nrow=num_piest,ncol=N)
pred_pi2_3=matrix(0,nrow=num_piest,ncol=N)
pred_mu2_3=matrix(0,nrow=num_mu1est,ncol=N)
pred_mu1_3_DTL=matrix(0,nrow=num_mu1est,ncol=N)
pred_mu1_3_DDRL=matrix(0,nrow=num_mu1est,ncol=N)

for (k in 1:K){
  if (k==K) {
    index=(((k-1)*M+1):N)
  }
  if (k<K) {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  
  index_nk1=index_nk[1:floor(length(index_nk)/2)]
  index3_nk1=index_nk1[which(A1_3[index_nk1]==1)]
  index33_nk1=index_nk1[which(A1_3[index_nk1]*A2_3[index_nk1]==1)]
  
  index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  index3_nk2=index_nk2[which(A1_3[index_nk2]==1)]
  index33_nk2=index_nk2[which(A1_3[index_nk2]*A2_3[index_nk2]==1)]
  
  index3_nk=index_nk[which(A1_3[index_nk]==1)]
  index33_nk=index_nk[which(A1_3[index_nk]*A2_3[index_nk]==1)]
  
  #pi1
  fit_Log1_3=cv.glmnet(x=S1[index_nk,],y=A1_3[index_nk],family="binomial")
  pred_pi1_3[1,index]=predict(fit_Log1_3,newx=S1[index,],type="response",s="lambda.min") #logistic+lasso
  #pi2_1
  fit_Log2_3=cv.glmnet(x=cbind(S1[index3_nk,],S2[index3_nk,]),y=A2_3[index3_nk],family="binomial")
  pred_pi2_3[1,index]=predict(fit_Log2_3,newx=cbind(S1[index,],S2[index,]),type="response",s="lambda.min") #logistic
  #mu2_1
  fit_lasso1_3=cv.glmnet(x=cbind(S1[index33_nk,],S2[index33_nk,]),y=Y[index33_nk],family="gaussian",alpha=1)
  pred_mu2_3[1,index]=predict(fit_lasso1_3,newx=cbind(S1[index,],S2[index,]),s="lambda.min") #lasso
  #mu1_1
  pred_mu2_3_nocross=predict(fit_lasso1_3,newx=cbind(S1[index3_nk,],S2[index3_nk,]),s="lambda.min")
  fit_lasso2_3_nocross=cv.glmnet(x=S1[index3_nk,],y=pred_mu2_3_nocross,family="gaussian",alpha=1)
  pred_mu1_3_DTL[1,index]=predict(fit_lasso2_3_nocross,newx=S1[index,],s="lambda.min")
  
  
  
  #mu1_1
  fit_lasso1_3_new1=cv.glmnet(x=cbind(S1[index33_nk1,],S2[index33_nk1,]),y=Y[index33_nk1],family="gaussian",alpha=1)
  fit_Log2_3_new1=cv.glmnet(x=cbind(S1[index3_nk1,],S2[index3_nk1,]),y=A2_3[index3_nk1],family="binomial")
  pred_mu2_3_new1=predict(fit_lasso1_3_new1,newx=cbind(S1[index3_nk2,],S2[index3_nk2,]),s="lambda.min")
  pred_pi2_3_new1=predict(fit_Log2_3_new1,newx=cbind(S1[index3_nk2,],S2[index3_nk2,]),type="response",s="lambda.min")
  Y1_new1=as.numeric(pred_mu2_3_new1+A2_3[index3_nk2]*(Y[index3_nk2]-pred_mu2_3_new1)/pred_pi2_3_new1)
  fit_lasso2_nocross3_new1=cv.glmnet(x=S1[index3_nk2,],y=Y1_new1,family="gaussian",alpha=1)
  pred_mu1_3_new1=predict(fit_lasso2_nocross3_new1,newx=S1[index,],s="lambda.min")
  
  fit_lasso1_3_new2=cv.glmnet(x=cbind(S1[index33_nk2,],S2[index33_nk2,]),y=Y[index33_nk2],family="gaussian",alpha=1)
  fit_Log2_3_new2=cv.glmnet(x=cbind(S1[index3_nk2,],S2[index3_nk2,]),y=A2_3[index3_nk2],family="binomial")
  pred_mu2_3_new2=predict(fit_lasso1_3_new2,newx=cbind(S1[index3_nk1,],S2[index3_nk1,]),s="lambda.min")
  pred_pi2_3_new2=predict(fit_Log2_3_new2,newx=cbind(S1[index3_nk1,],S2[index3_nk1,]),type="response",s="lambda.min")
  Y1_new2=as.numeric(pred_mu2_3_new2+A2_3[index3_nk1]*(Y[index3_nk1]-pred_mu2_3_new2)/pred_pi2_3_new2)
  fit_lasso2_nocross3_new2=cv.glmnet(x=S1[index3_nk1,],y=Y1_new2,family="gaussian",alpha=1)
  pred_mu1_3_new2=predict(fit_lasso2_nocross3_new2,newx=S1[index,],s="lambda.min")
  
  pred_mu1_3_DDRL[1,index]=(pred_mu1_3_new1+pred_mu1_3_new2)/2
}

pred_pi1_1=matrix(0,nrow=num_piest,ncol=N)
pred_pi2_1=matrix(0,nrow=num_piest,ncol=N)
pred_mu2_1=matrix(0,nrow=num_mu1est,ncol=N)
pred_mu1_1_DTL=matrix(0,nrow=num_mu1est,ncol=N)
pred_mu1_1_DDRL=matrix(0,nrow=num_mu1est,ncol=N)


for (k in 1:K){
  if (k==K) {
    index=(((k-1)*M+1):N)
  }
  if (k<K) {
    index=(1:M)+(k-1)*M
  }
  index_nk=(1:N)[-index]
  
  index_nk1=index_nk[1:floor(length(index_nk)/2)]
  index1_nk1=index_nk1[which(A1_1[index_nk1]==1)]
  index11_nk1=index_nk1[which(A1_1[index_nk1]*A2_1[index_nk1]==1)]
  
  index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  index1_nk2=index_nk2[which(A1_1[index_nk2]==1)]
  index11_nk2=index_nk2[which(A1_1[index_nk2]*A2_1[index_nk2]==1)]
  
  index1_nk=index_nk[which(A1_1[index_nk]==1)]
  index11_nk=index_nk[which(A1_1[index_nk]*A2_1[index_nk]==1)]
  #pi1
  fit_Log1_1=cv.glmnet(x=S1[index_nk,],y=A1_1[index_nk],family="binomial")
  pred_pi1_1[1,index]=predict(fit_Log1_1,newx=S1[index,],type="response",s="lambda.min") #logistic+lasso
  #pi2_1
  fit_Log2_1=cv.glmnet(x=cbind(S1[index1_nk,],S2[index1_nk,]),y=A2_1[index1_nk],family="binomial")
  pred_pi2_1[1,index]=predict(fit_Log2_1,newx=cbind(S1[index,],S2[index,]),type="response",s="lambda.min") #logistic
  #mu2_1
  fit_lasso1_1=cv.glmnet(x=cbind(S1[index11_nk,],S2[index11_nk,]),y=Y[index11_nk],family="gaussian",alpha=1)
  pred_mu2_1[1,index]=predict(fit_lasso1_1,newx=cbind(S1[index,],S2[index,]),s="lambda.min") #lasso
  #mu1_1
  pred_mu2_1_nocross=predict(fit_lasso1_1,newx=cbind(S1[index1_nk,],S2[index1_nk,]),s="lambda.min")
  fit_lasso2_1_nocross=cv.glmnet(x=S1[index1_nk,],y=pred_mu2_1_nocross,family="gaussian",alpha=1)
  pred_mu1_1_DTL[1,index]=predict(fit_lasso2_1_nocross,newx=S1[index,],s="lambda.min")
  
  
  #mu1_1
  fit_lasso1_1_new1=cv.glmnet(x=cbind(S1[index11_nk1,],S2[index11_nk1,]),y=Y[index11_nk1],family="gaussian",alpha=1)
  fit_Log2_1_new1=cv.glmnet(x=cbind(S1[index1_nk1,],S2[index1_nk1,]),y=A2_1[index1_nk1],family="binomial")
  pred_mu2_1_new1=predict(fit_lasso1_1_new1,newx=cbind(S1[index1_nk2,],S2[index1_nk2,]),s="lambda.min")
  pred_pi2_1_new1=predict(fit_Log2_1_new1,newx=cbind(S1[index1_nk2,],S2[index1_nk2,]),type="response",s="lambda.min")
  Y1_new1=as.numeric(pred_mu2_1_new1+A2_1[index1_nk2]*(Y[index1_nk2]-pred_mu2_1_new1)/pred_pi2_1_new1)
  fit_lasso2_nocross1_new1=cv.glmnet(x=S1[index1_nk2,],y=Y1_new1,family="gaussian",alpha=1)
  pred_mu1_1_new1=predict(fit_lasso2_nocross1_new1,newx=S1[index,],s="lambda.min")
  
  fit_lasso1_1_new2=cv.glmnet(x=cbind(S1[index11_nk2,],S2[index11_nk2,]),y=Y[index11_nk2],family="gaussian",alpha=1)
  fit_Log2_1_new2=cv.glmnet(x=cbind(S1[index1_nk2,],S2[index1_nk2,]),y=A2_1[index1_nk2],family="binomial")
  pred_mu2_1_new2=predict(fit_lasso1_1_new2,newx=cbind(S1[index1_nk1,],S2[index1_nk1,]),s="lambda.min")
  pred_pi2_1_new2=predict(fit_Log2_1_new2,newx=cbind(S1[index1_nk1,],S2[index1_nk1,]),type="response",s="lambda.min")
  Y1_new2=as.numeric(pred_mu2_1_new2+A2_1[index1_nk1]*(Y[index1_nk1]-pred_mu2_1_new2)/pred_pi2_1_new2)
  fit_lasso2_nocross1_new2=cv.glmnet(x=S1[index1_nk1,],y=Y1_new2,family="gaussian",alpha=1)
  pred_mu1_1_new2=predict(fit_lasso2_nocross1_new2,newx=S1[index,],s="lambda.min")
  
  pred_mu1_1_DDRL[1,index]=(pred_mu1_1_new1+pred_mu1_1_new2)/2
}

#DTL trimmed
psi3_DTL=c()
psi1_DTL=c()
for (i in 1:length(pred_pi1_3[1,])) {
  if (pred_pi1_3[,i]*pred_pi2_3[,i]>0 & pred_pi1_1[,i]*pred_pi2_1[,i]>0) {
    gamma2_3=A1_3[i]*A2_3[i]/(pred_pi1_3[,i]*pred_pi2_3[,i])
    gamma1_3=A1_3[i]/pred_pi1_3[,i]
    psi_3_DTL=gamma2_3*Y[i]-(gamma1_3-1)*pred_mu1_3_DTL[,i]-(gamma2_3-gamma1_3)*pred_mu2_3[,i]
    psi3_DTL=c(psi3_DTL,psi_3_DTL)
    gamma2_1=A1_1[i]*A2_1[i]/(pred_pi1_1[,i]*pred_pi2_1[,i])
    gamma1_1=A1_1[i]/pred_pi1_1[,i]
    psi_1_DTL=gamma2_1*Y[i]-(gamma1_1-1)*pred_mu1_1_DTL[,i]-(gamma2_1-gamma1_1)*pred_mu2_1[,i]
    psi1_DTL=c(psi1_DTL,psi_1_DTL)
  }
}
length(psi3_DTL); length(psi1_DTL)
psi31_DTL=psi3_DTL-psi1_DTL
pred_theta_DTL=mean(psi31_DTL)
ob_DTL=length(psi31_DTL)
trimmed_DTL=length(A2_all)-ob_DTL
pred_theta0_DTL=mean(psi1_DTL)


pred_sig2_DTL=c()
for (i in 1:length(pred_pi1_1[1,])) {
  if (pred_pi1_3[,i]*pred_pi2_3[,i]>0 & pred_pi1_1[,i]*pred_pi2_1[,i]>0) {
    gamma2_3=A1_3[i]*A2_3[i]/(pred_pi1_3[,i]*pred_pi2_3[,i])
    gamma1_3=A1_3[i]/pred_pi1_3[,i]
    psi_3_DTL=gamma2_3*Y[i]-(gamma1_3-1)*pred_mu1_3_DTL[,i]-(gamma2_3-gamma1_3)*pred_mu2_3[,i]
    gamma2_1=A1_1[i]*A2_1[i]/(pred_pi1_1[,i]*pred_pi2_1[,i])
    gamma1_1=A1_1[i]/pred_pi1_1[,i]
    psi_1_DTL=gamma2_1*Y[i]-(gamma1_1-1)*pred_mu1_1_DTL[,i]-(gamma2_1-gamma1_1)*pred_mu2_1[,i]
    psi=psi_3_DTL-psi_1_DTL-pred_theta_DTL
    pred_sig2_DTL=c(pred_sig2_DTL,psi^2)
  }
}
length(pred_sig2_DTL)

SE_DTL=sqrt(mean(pred_sig2_DTL))/sqrt(ob_DTL)
p_value_DTL=2*pnorm(-abs(pred_theta_DTL)/SE_DTL)
observations=length(A2_all)


#DDRL trimmed
psi3_DDRL=c()
psi1_DDRL=c()
for (i in 1:length(pred_pi1_3[1,])) {
  if (pred_pi1_3[,i]*pred_pi2_3[,i]>0 & pred_pi1_1[,i]*pred_pi2_1[,i]>0) {
    gamma2_3=A1_3[i]*A2_3[i]/(pred_pi1_3[,i]*pred_pi2_3[,i])
    gamma1_3=A1_3[i]/pred_pi1_3[,i]
    psi_3_DDRL=gamma2_3*Y[i]-(gamma1_3-1)*pred_mu1_3_DDRL[,i]-(gamma2_3-gamma1_3)*pred_mu2_3[,i]
    psi3_DDRL=c(psi3_DDRL,psi_3_DDRL)
    gamma2_1=A1_1[i]*A2_1[i]/(pred_pi1_1[,i]*pred_pi2_1[,i])
    gamma1_1=A1_1[i]/pred_pi1_1[,i]
    psi_1_DDRL=gamma2_1*Y[i]-(gamma1_1-1)*pred_mu1_1_DDRL[,i]-(gamma2_1-gamma1_1)*pred_mu2_1[,i]
    psi1_DDRL=c(psi1_DDRL,psi_1_DDRL)
  }
}
length(psi3_DDRL); length(psi1_DDRL)
psi31_DDRL=psi3_DDRL-psi1_DDRL
pred_theta_DDRL=mean(psi31_DDRL)
ob_DDRL=length(psi31_DDRL)
trimmed_DDRL=length(A2_all)-ob_DDRL
pred_theta0_DDRL=mean(psi1_DDRL)


pred_sig2_DDRL=c()
for (i in 1:length(pred_pi1_1[1,])) {
  if (pred_pi1_3[,i]*pred_pi2_3[,i]>0 & pred_pi1_1[,i]*pred_pi2_1[,i]>0) {
    gamma2_3=A1_3[i]*A2_3[i]/(pred_pi1_3[,i]*pred_pi2_3[,i])
    gamma1_3=A1_3[i]/pred_pi1_3[,i]
    psi_3_DDRL=gamma2_3*Y[i]-(gamma1_3-1)*pred_mu1_3_DDRL[,i]-(gamma2_3-gamma1_3)*pred_mu2_3[,i]
    gamma2_1=A1_1[i]*A2_1[i]/(pred_pi1_1[,i]*pred_pi2_1[,i])
    gamma1_1=A1_1[i]/pred_pi1_1[,i]
    psi_1_DDRL=gamma2_1*Y[i]-(gamma1_1-1)*pred_mu1_1_DDRL[,i]-(gamma2_1-gamma1_1)*pred_mu2_1[,i]
    psi=psi_3_DDRL-psi_1_DDRL-pred_theta_DDRL
    pred_sig2_DDRL=c(pred_sig2_DDRL,psi^2)
  }
}
length(pred_sig2_DDRL)

SE_DDRL=sqrt(mean(pred_sig2_DDRL))/sqrt(ob_DDRL)
p_value_DDRL=2*pnorm(-abs(pred_theta_DDRL)/SE_DDRL)
observations=length(A2_all)





pred_theta0_DTL;pred_theta_DTL;SE_DTL;p_value_DTL;observations;trimmed_DTL
pred_theta0_DDRL;pred_theta_DDRL;SE_DDRL;p_value_DDRL;observations;trimmed_DDRL


max_3=max(c(length(as.numeric(pred_pi1_3)[A1_3==1]), length(as.numeric(pred_pi1_3)[A1_3==0])))
df3=data.frame(col1=c(as.numeric(pred_pi1_3)[A1_3==1],rep(NA, max_3 - length(as.numeric(pred_pi1_3)[A1_3==1]))),
               col2=c(as.numeric(pred_pi1_3)[A1_3==0],rep(NA, max_3 - length(as.numeric(pred_pi1_3)[A1_3==0]))))
ggplot(df3, aes(x=x)) + 
  geom_histogram(aes(x=df3[,1], y = ..density..), colour = 1, fill = "white") +
  geom_histogram(aes(x=df3[,2], y = -..density..), colour = 1, fill = "white") +
  xlab("Propensity Score")+
  ylab("Density")

max_1=max(c(length(as.numeric(pred_pi1_1)[A1_1==1]), length(as.numeric(pred_pi1_1)[A1_1==0])))
df1=data.frame(col1=c(as.numeric(pred_pi1_1)[A1_1==1],rep(NA, max_1 - length(as.numeric(pred_pi1_1)[A1_1==1]))),
               col2=c(as.numeric(pred_pi1_1)[A1_1==0],rep(NA, max_1 - length(as.numeric(pred_pi1_1)[A1_1==0]))))
ggplot(df1, aes(x =x)) + 
  geom_histogram(aes(x=df1[,1], y = ..density..), colour = 1, fill = "white") +
  geom_histogram(aes(x=df1[,2], y = -..density..), colour = 1, fill = "white") +
  xlab("Propensity Score")+
  ylab("Density")

max_33=max(c(length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1]), length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0])))
df33=data.frame(col1=c(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1],rep(NA, max_33 - length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1]))),
                col2=c(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0],rep(NA, max_33 - length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0]))))
ggplot(df33, aes(x=x)) + 
  geom_histogram(aes(x=df33[,1], y = ..density..), colour = 1, fill = "white") +
  geom_histogram(aes(x=df33[,2], y = -..density..), colour = 1, fill = "white") +
  xlab("Propensity Score")+
  ylab("Density")

max_11=max(c(length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1]), length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0])))
df11=data.frame(col1=c(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1],rep(NA, max_11 - length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1]))),
                col2=c(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0],rep(NA, max_11 - length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0]))))
ggplot(df11, aes(x=x)) + 
  geom_histogram(aes(x=df11[,1], y = ..density..), colour = 1, fill = "white") +
  geom_histogram(aes(x=df11[,2], y = -..density..), colour = 1, fill = "white") +
  xlab("Propensity Score")+
  ylab("Density")









max_3=max(c(length(as.numeric(pred_pi1_3)[A1_3==1]), length(as.numeric(pred_pi1_3)[A1_3==0])))
df3=data.frame(col1=c(as.numeric(pred_pi1_3)[A1_3==1],rep(NA, max_3 - length(as.numeric(pred_pi1_3)[A1_3==1]))),
               col2=c(as.numeric(pred_pi1_3)[A1_3==0],rep(NA, max_3 - length(as.numeric(pred_pi1_3)[A1_3==0]))))
ggplot(df3, aes(x=x)) + 
  geom_histogram(aes(x=df3[,1], y = ..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.4, y=5, label="treated"),colour = 1)+
  geom_histogram(aes(x=df3[,2], y = -..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.4, y=-5, label="controls"),colour = 1)+
  xlab("Propensity Score")+
  ylab("Density")

max_1=max(c(length(as.numeric(pred_pi1_1)[A1_1==1]), length(as.numeric(pred_pi1_1)[A1_1==0])))
df1=data.frame(col1=c(as.numeric(pred_pi1_1)[A1_1==1],rep(NA, max_1 - length(as.numeric(pred_pi1_1)[A1_1==1]))),
               col2=c(as.numeric(pred_pi1_1)[A1_1==0],rep(NA, max_1 - length(as.numeric(pred_pi1_1)[A1_1==0]))))
ggplot(df1, aes(x =x)) + 
  geom_histogram(aes(x=df1[,1], y = ..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.5, y=5, label="treated"),colour = 1)+
  geom_histogram(aes(x=df1[,2], y = -..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.5, y=-4, label="controls"),colour = 1)+
  xlab("Propensity Score")+
  ylab("Density")

max_33=max(c(length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1]), length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0])))
df33=data.frame(col1=c(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1],rep(NA, max_33 - length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==1]))),
                col2=c(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0],rep(NA, max_33 - length(as.numeric(pred_pi2_3)[A1_3==1 & A2_3==0]))))
ggplot(df33, aes(x=x)) + 
  geom_histogram(aes(x=df33[,1], y = ..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.8, y=2.5, label="treated"),colour = 1)+
  geom_histogram(aes(x=df33[,2], y = -..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.8, y=-2, label="controls"),colour = 1)+
  xlab("Propensity Score")+
  ylab("Density")

max_11=max(c(length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1]), length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0])))
df11=data.frame(col1=c(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1],rep(NA, max_11 - length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==1]))),
                col2=c(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0],rep(NA, max_11 - length(as.numeric(pred_pi2_1)[A1_1==1 & A2_1==0]))))
ggplot(df11, aes(x=x)) + 
  geom_histogram(aes(x=df11[,1], y = ..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.75, y=8, label="treated"),colour = 1)+
  geom_histogram(aes(x=df11[,2], y = -..density..), colour = 1, fill = "white") +
  geom_label( aes(x=0.75, y=-10, label="controls"),colour = 1)+
  xlab("Propensity Score")+
  ylab("Density")

