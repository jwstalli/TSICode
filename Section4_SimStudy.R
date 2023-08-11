#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)


##################################
#Effect generation information
##################################

m_main<-6
m_2FI<-5
m_quad<-3
eff_main_offset<-2.5
eff_main_rate<-1
eff_int_offset<-2.5
eff_int_rate<-1
eff_quad_offset<-2.5
eff_quad_rate<-1

err_sigma<-1

true.heredity<-"strong"
model<-"FullQuad"



##################################
#Generate effects and error terms
##################################


num_sims<-1000
k<-8
n<-21

if(model=="ME"){
  p=k
}
if(model=="2FI"){
  p=k+choose(k,2)
}
if(model=="FullQuad"){
  p=2*k+choose(k,2)
}

set.seed(1234)



active<-matrix(0,nrow=num_sims,ncol=k)
for(i in 1:num_sims){
  active[i,sample(1:k,m_main)]<-1
}
errors<-matrix(rnorm(num_sims*n,0,err_sigma),nrow=num_sims,ncol=n)

betas<-matrix(0,nrow=num_sims,ncol=p)
model_terms<-matrix(0,nrow=num_sims,ncol=p)

for(i in 1:num_sims){
  MeanModel<-ModelSupportAndEffects(active=active[i,],
                                    model=model,
                                    heredity = true.heredity,
                                    m_2FI=m_2FI,
                                    m_quad=m_quad,
                                    eff_main_offset = eff_main_offset,
                                    eff_main_rate = eff_main_rate,
                                    eff_int_offset = eff_int_offset ,
                                    eff_int_rate = eff_int_rate,
                                    eff_quad_offset = eff_quad_offset,
                                    eff_quad_rate = eff_quad_rate)  
  betas[i,]=MeanModel$b
  model_terms[i,]=MeanModel$model_terms
}



if(k==6){
  D<-matrix(c(0,1,1,1,1,1,
              0,-1,-1,-1,-1,-1,
              1,0,1,1,-1,1,
              -1,0,-1,-1,1,-1,
              1,-1,0,1,1,-1,
              -1,1,0,-1,-1,1,
              1,-1,-1,0,1,1,
              -1,1,1,0,-1,-1,
              1,1,-1,-1,0,1,
              -1,-1,1,1,0,-1,
              1,-1,1,-1,-1,0,
              -1,1,-1,1,1,0,
              1,1,-1,1,-1,-1,
              -1,-1,1,-1,1,1,
              1,1,1,-1,1,-1,
              -1,-1,-1,1,-1,1,
              0,0,0,0,0,0),nrow=17,ncol=6,byrow=TRUE)
}else{
  D<-matrix(c(0,1,1,1,1,1,1,1,
              0,-1,-1,-1,-1,-1,-1,-1,
              1,0,-1,-1,-1,-1,1,1,
              -1,0,1,1,1,1,-1,-1,
              1,-1,0,-1,1,1,-1,-1,
              -1,1,0,1,-1,-1,1,1,
              1,-1,-1,0,1,1,1,1,
              -1,1,1,0,-1,-1,-1,-1,
              1,-1,1,1,0,-1,-1,1,
              -1,1,-1,-1,0,1,1,-1,
              1,-1,1,1,-1,0,1,-1,
              -1,1,-1,-1,1,0,-1,1,
              1,1,-1,1,-1,1,0,-1,
              -1,-1,1,-1,1,-1,0,1,
              1,1,-1,1,1,-1,-1,0,
              -1,-1,1,-1,-1,1,1,0,
              1,1,1,-1,-1,1,-1,1,
              -1,-1,-1,1,1,-1,1,-1,
              1,1,1,-1,1,-1,1,-1,
              -1,-1,-1,1,-1,1,-1,1,
              0,0,0,0,0,0,0,0),nrow=21,ncol=8,byrow=TRUE)
}










##################################
#Create Model Matrix (sans intercept)
##################################

if(model=="ME"){
  X=D
}



if(model=="2FI"){
  X=D
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      X<-cbind(X,D[,i]*D[,j])
    }
  }
}


if(model=="FullQuad"){
  X=D
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      X<-cbind(X,D[,i]*D[,j])
    }
  }
  
  X<-cbind(X,D^2)
  
}


alphaME=0.05
pool="ME"
heredity="strong"
maxpower=FALSE
if(maxpower==TRUE){
  alpha2nd=0.05
}else{
  alpha2nd=0.20
}



Results<-matrix(NA,nrow=3,ncol=4)



###################################
#Guided subsets maxmod = c/2
###################################

GS1.all.model.terms<-matrix(NA,nrow=num_sims,ncol=ncol(X))

ME.TPR<-rep(NA,num_sims)
ME.FPR<-rep(NA,num_sims)
TPR<-rep(NA,num_sims)
FPR<-rep(NA,num_sims)
MOEs.mean<-rep(NA,num_sims)
MSE<-rep(NA,num_sims)
MSE.DF<-rep(NA,num_sims)


for(iter in 1:num_sims){
  y<-10+X%*%betas[iter,]+errors[iter,]
  
  
  GS1.Analysis<-GuidedSubsetsAnalysis(y=y,D=D,model=model,heredity=heredity,alphaME=alphaME,alpha2nd=alpha2nd,pool=pool,
                                     maxpower=maxpower,maxmod="0.5")  
  
  GS1.all.model.terms[iter,]<-GS1.Analysis$est_model_terms
  
  
  #True Model Information
  active_ME_terms<-which(model_terms[iter,1:k]==1)
  n.active.ME<-length(active_ME_terms)
  inert_ME_terms<-which(model_terms[iter,1:k]==0)
  n.inert.ME<-length(inert_ME_terms)
  
  if(n.active.ME>0){
    ME.TPR[iter]<-length(intersect(which(GS1.Analysis$est_model_terms[1:k]==1),active_ME_terms))/n.active.ME
  }else{
    ME.TPR[iter]=1
  }
  
  if(n.inert.ME>0){
    ME.FPR[iter]<-length(intersect(which(GS1.Analysis$est_model_terms[1:k]==1),inert_ME_terms))/n.inert.ME
  }else{
    ME.FPR[iter]<-0
  }  
  
  
  active_model_terms<-which(model_terms[iter,]==1)
  n.active<-length(active_model_terms)
  inert_model_terms<-which(model_terms[iter,]==0)
  n.inert<-length(inert_model_terms)
  
  if(n.active>0){
    TPR[iter]<-length(intersect(which(GS1.Analysis$est_model_terms==1),active_model_terms))/n.active
  }else{
    TPR[iter]=1
  }
  
  if(n.inert>0){
    FPR[iter]<-length(intersect(which(GS1.Analysis$est_model_terms==1),inert_model_terms))/n.inert
  }else{
    FPR[iter]<-0
  }
  
  
  MOEs<-apply(GS1.Analysis$CIs,1,diff)/2  
  MOEs.mean[iter]<-mean(MOEs)
  
  MSE[iter]<-GS1.Analysis$PureMSE
  MSE.DF[iter]<-GS1.Analysis$PureDF
  
  print(iter)
  #if(iter%%20==0) print(iter)
}



c(mean(ME.TPR),mean(ME.FPR),mean(TPR),mean(FPR),mean(MOEs.mean),mean(MSE),mean(MSE.DF))



#2FIs
indx<-which(model_terms[,c(1:choose(k,2))+k]==1)
if(length(indx)==0){
  print(0)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS1.all.model.terms[,c(1:choose(k,2))+k])==-1))/(choose(k,2)*num_sims)
}else{
  #TP
  TP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS1.all.model.terms[,c(1:choose(k,2))+k])[indx]==0))/length(indx)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS1.all.model.terms[,c(1:choose(k,2))+k])[-indx]==-1))/(choose(k,2)*num_sims-length(indx))
}

Results[1,1:2]<-c(TP,FP)


#Quadratic Effects
indx<-which(model_terms[,(ncol(X)-(k-1)):ncol(X)]==1)
#TP
TP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-GS1.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[indx]==0))/length(indx)
#FP
FP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-GS1.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[-indx]==-1))/(k*num_sims-length(indx))

Results[1,3:4]<-c(TP,FP)



###################################
#Guided subsets maxmod = c-1
###################################

GS2.all.model.terms<-matrix(NA,nrow=num_sims,ncol=ncol(X))

ME.TPR<-rep(NA,num_sims)
ME.FPR<-rep(NA,num_sims)
TPR<-rep(NA,num_sims)
FPR<-rep(NA,num_sims)
MOEs.mean<-rep(NA,num_sims)
MSE<-rep(NA,num_sims)
MSE.DF<-rep(NA,num_sims)


for(iter in 1:num_sims){
  y<-10+X%*%betas[iter,]+errors[iter,]
  
  
  GS2.Analysis<-GuidedSubsetsAnalysis(y=y,D=D,model=model,heredity=heredity,alphaME=alphaME,alpha2nd=alpha2nd,pool=pool,
                                  maxpower=maxpower,maxmod="r-1")  
  
  GS2.all.model.terms[iter,]<-GS2.Analysis$est_model_terms
  
  
  #True Model Information
  active_ME_terms<-which(model_terms[iter,1:k]==1)
  n.active.ME<-length(active_ME_terms)
  inert_ME_terms<-which(model_terms[iter,1:k]==0)
  n.inert.ME<-length(inert_ME_terms)
  
  if(n.active.ME>0){
    ME.TPR[iter]<-length(intersect(which(GS2.Analysis$est_model_terms[1:k]==1),active_ME_terms))/n.active.ME
  }else{
    ME.TPR[iter]=1
  }
  
  if(n.inert.ME>0){
    ME.FPR[iter]<-length(intersect(which(GS2.Analysis$est_model_terms[1:k]==1),inert_ME_terms))/n.inert.ME
  }else{
    ME.FPR[iter]<-0
  }  
  
  
  active_model_terms<-which(model_terms[iter,]==1)
  n.active<-length(active_model_terms)
  inert_model_terms<-which(model_terms[iter,]==0)
  n.inert<-length(inert_model_terms)
  
  if(n.active>0){
    TPR[iter]<-length(intersect(which(GS2.Analysis$est_model_terms==1),active_model_terms))/n.active
  }else{
    TPR[iter]=1
  }
  
  if(n.inert>0){
    FPR[iter]<-length(intersect(which(GS2.Analysis$est_model_terms==1),inert_model_terms))/n.inert
  }else{
    FPR[iter]<-0
  }
  
  
  MOEs<-apply(GS2.Analysis$CIs,1,diff)/2  
  MOEs.mean[iter]<-mean(MOEs)
  
  MSE[iter]<-GS2.Analysis$PureMSE
  MSE.DF[iter]<-GS2.Analysis$PureDF
  
  print(iter)
  #if(iter%%20==0) print(iter)
}



c(mean(ME.TPR),mean(ME.FPR),mean(TPR),mean(FPR),mean(MOEs.mean),mean(MSE),mean(MSE.DF))



#2FIs
indx<-which(model_terms[,c(1:choose(k,2))+k]==1)
if(length(indx)==0){
  print(0)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS2.all.model.terms[,c(1:choose(k,2))+k])==-1))/(choose(k,2)*num_sims)
}else{
  #TP
  TP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS2.all.model.terms[,c(1:choose(k,2))+k])[indx]==0))/length(indx)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-GS2.all.model.terms[,c(1:choose(k,2))+k])[-indx]==-1))/(choose(k,2)*num_sims-length(indx))
}

Results[2,1:2]<-c(TP,FP)


#Quadratic Effects
indx<-which(model_terms[,(ncol(X)-(k-1)):ncol(X)]==1)
#TP
TP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-GS2.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[indx]==0))/length(indx)
#FP
FP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-GS2.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[-indx]==-1))/(k*num_sims-length(indx))

Results[2,3:4]<-c(TP,FP)




############################################
# All subsets
############################################

AS.all.model.terms<-matrix(NA,nrow=num_sims,ncol=ncol(X))

ME.TPR<-rep(NA,num_sims)
ME.FPR<-rep(NA,num_sims)
TPR<-rep(NA,num_sims)
FPR<-rep(NA,num_sims)
MOEs.mean<-rep(NA,num_sims)
MSE<-rep(NA,num_sims)
MSE.DF<-rep(NA,num_sims)


for(iter in 1:num_sims){
  y<-10+X%*%betas[iter,]+errors[iter,]
  
  
  AS.Analysis<-AllSubsetsAnalysis(y=y,D=D,model=model,heredity=heredity,alphaME=alphaME,pool=pool)  
  
  AS.all.model.terms[iter,]<-AS.Analysis$est_model_terms
  
  active_ME_terms<-which(model_terms[iter,1:k]==1)
  n.active.ME<-length(active_ME_terms)
  inert_ME_terms<-which(model_terms[iter,1:k]==0)
  n.inert.ME<-length(inert_ME_terms)
  
  if(n.active.ME>0){
    ME.TPR[iter]<-length(intersect(which(AS.Analysis$est_model_terms[1:k]==1),active_ME_terms))/n.active.ME
  }else{
    ME.TPR[iter]=1
  }
  
  if(n.inert.ME>0){
    ME.FPR[iter]<-length(intersect(which(AS.Analysis$est_model_terms[1:k]==1),inert_ME_terms))/n.inert.ME
  }else{
    ME.FPR[iter]<-0
  }  
  
  
  active_model_terms<-which(model_terms[iter,]==1)
  n.active<-length(active_model_terms)
  inert_model_terms<-which(model_terms[iter,]==0)
  n.inert<-length(inert_model_terms)
  
  if(n.active>0){
    TPR[iter]<-length(intersect(which(AS.Analysis$est_model_terms==1),active_model_terms))/n.active
  }else{
    TPR[iter]=1
  }
  
  if(n.inert>0){
    FPR[iter]<-length(intersect(which(AS.Analysis$est_model_terms==1),inert_model_terms))/n.inert
  }else{
    FPR[iter]<-0
  }
  
  
  MOEs<-apply(AS.Analysis$CIs,1,diff)/2  
  MOEs.mean[iter]<-mean(MOEs)
  
  MSE[iter]<-AS.Analysis$PureMSE
  MSE.DF[iter]<-AS.Analysis$PureDF
  
  print(iter)
  #if(iter%%25==0) print(iter)
}



c(mean(ME.TPR),mean(ME.FPR),mean(TPR),mean(FPR),mean(MOEs.mean),mean(MSE),mean(MSE.DF))



#2FIs
indx<-which(model_terms[,c(1:choose(k,2))+k]==1)
if(length(indx)==0){
  print(0)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-AS.all.model.terms[,c(1:choose(k,2))+k])==-1))/(choose(k,2)*num_sims)
}else{
  #TP
  TP=length(which( (model_terms[,c(1:choose(k,2))+k]-AS.all.model.terms[,c(1:choose(k,2))+k])[indx]==0))/length(indx)
  #FP
  FP=length(which( (model_terms[,c(1:choose(k,2))+k]-AS.all.model.terms[,c(1:choose(k,2))+k])[-indx]==-1))/(choose(k,2)*num_sims-length(indx))
}

Results[3,1:2]<-c(TP,FP)


#Quadratic Effects
indx<-which(model_terms[,(ncol(X)-(k-1)):ncol(X)]==1)
#TP
TP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-AS.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[indx]==0))/length(indx)
#FP
FP=length(which( (model_terms[,(ncol(X)-(k-1)):ncol(X)]-AS.all.model.terms[,(ncol(X)-(k-1)):ncol(X)])[-indx]==-1))/(k*num_sims-length(indx))

Results[3,3:4]<-c(TP,FP)




#Results from 1000 sims had to be hardcoded

k6.n17.ME4.2FI2.Q1<-matrix(c(0.9780, 0.01400000, 0.840, 0.0258,
                             0.9630, 0.02300000, 0.835, 0.0412,
                             0.9495, 0.05607692, 0.928, 0.1098),nrow=3,ncol=4,byrow=TRUE)

k6.n17.ME4.2FI4.Q3<-matrix(c(0.19150, 0.01545455, 0.08666667, 0.03833333,
                             0.59350, 0.08227273, 0.51600000, 0.20700000,
                             0.62275, 0.09000000, 0.57000000, 0.22800000),nrow=3,ncol=4,byrow=TRUE)


k8.n21.ME6.2FI3.Q2<-matrix(c(0.7776667, 0.03348, 0.466, 0.0570000,
                             0.7560000, 0.04116, 0.453, 0.0675000,
                             0.7160000, 0.07136, 0.522, 0.1153333),nrow=3,ncol=4,byrow=TRUE)


k8.n21.ME6.2FI5.Q3<-matrix(c(0.3836, 0.07682609, 0.206, 0.1140,
                             0.3932, 0.09182609, 0.257, 0.1422,
                             0.4090, 0.11495652, 0.299, 0.1684),nrow=3,ncol=4,byrow=TRUE)






par(mfrow=c(1,2))

plot(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[1,1], k6.n17.ME4.2FI4.Q3[1,1], k8.n21.ME6.2FI3.Q2[1,1], k8.n21.ME6.2FI5.Q3[1,1]),type="l",lwd=4,ylim=c(0,1),
     ylab="Proportion",xlab="Cases",main="TPR/FPR 2FIs",axes=FALSE)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[2,1], k6.n17.ME4.2FI4.Q3[2,1], k8.n21.ME6.2FI3.Q2[2,1], k8.n21.ME6.2FI5.Q3[2,1]),lwd=3,col=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[3,1], k6.n17.ME4.2FI4.Q3[3,1], k8.n21.ME6.2FI3.Q2[3,1], k8.n21.ME6.2FI5.Q3[3,1]),lwd=3,col=3)
box()
axis(1,at=c(1:4))
axis(2)

lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[1,2], k6.n17.ME4.2FI4.Q3[1,2], k8.n21.ME6.2FI3.Q2[1,2], k8.n21.ME6.2FI5.Q3[1,2]),lwd=4,lty=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[2,2], k6.n17.ME4.2FI4.Q3[2,2], k8.n21.ME6.2FI3.Q2[2,2], k8.n21.ME6.2FI5.Q3[2,2]),lwd=3,col=2,lty=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[3,2], k6.n17.ME4.2FI4.Q3[3,2], k8.n21.ME6.2FI3.Q2[3,2], k8.n21.ME6.2FI5.Q3[3,2]),lwd=3,col=3,lty=2)



plot(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[1,3], k6.n17.ME4.2FI4.Q3[1,3], k8.n21.ME6.2FI3.Q2[1,3], k8.n21.ME6.2FI5.Q3[1,3]),type="l",lwd=4,ylim=c(0,1),
     ylab="Proportion",xlab="Cases",main="TPR/FPR Quads",axes=FALSE)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[2,3], k6.n17.ME4.2FI4.Q3[2,3], k8.n21.ME6.2FI3.Q2[2,3], k8.n21.ME6.2FI5.Q3[2,3]),lwd=3,col=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[3,3], k6.n17.ME4.2FI4.Q3[3,3], k8.n21.ME6.2FI3.Q2[3,3], k8.n21.ME6.2FI5.Q3[3,3]),lwd=3,col=3)
box()
axis(1,at=c(1:4))
axis(2)

lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[1,4], k6.n17.ME4.2FI4.Q3[1,4], k8.n21.ME6.2FI3.Q2[1,4], k8.n21.ME6.2FI5.Q3[1,4]),lwd=4,lty=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[2,4], k6.n17.ME4.2FI4.Q3[2,4], k8.n21.ME6.2FI3.Q2[2,4], k8.n21.ME6.2FI5.Q3[2,4]),lwd=3,col=2,lty=2)
lines(x=1:4,y=c(k6.n17.ME4.2FI2.Q1[3,4], k6.n17.ME4.2FI4.Q3[3,4], k8.n21.ME6.2FI3.Q2[3,4], k8.n21.ME6.2FI5.Q3[3,4]),lwd=3,col=3,lty=2)

