#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)

###################
#Simulation Study
###################
#Only choose one of the D objects below to run

#To run ADSD load this design
D<-matrix(c(0,1,1,1,1,1,1,
            0,-1,-1,-1,-1,-1,-1,
            1,0,1,-1,1,1,1,
            -1,0,-1,1,-1,-1,-1,
            1,-1,0,1,-1,1,1,
            -1,1,0,-1,1,-1,-1,
            1,1,-1,0,1,-1,1,
            -1,-1,1,0,-1,1,-1,
            1,-1,1,-1,0,1,-1,
            -1,1,-1,1,0,-1,1,
            1,-1,-1,1,-1,0,1,
            -1,1,1,-1,1,0,-1,
            1,-1,-1,-1,1,-1,0,
            -1,1,1,1,-1,1,0,
            1,1,-1,-1,-1,1,-1,
            -1,-1,1,1,1,-1,1,
            1,1,1,-1,-1,-1,1,
            -1,-1,-1,1,1,1,-1,
            1,1,1,1,-1,-1,-1,
            -1,-1,-1,-1,1,1,1,
            1,-1,1,1,1,-1,-1,
            -1,1,-1,-1,-1,1,1,
            1,1,-1,1,1,1,-1,
            -1,-1,1,-1,-1,-1,1),nrow=24,ncol=7,byrow=TRUE)

#To run optimal constrained LOF design, load this design

D<-matrix(c(1,-1,0,-1,1,0,-1,
            1,-1,0,1,-1,1,1,
            1,0,0,-1,1,1,1,
            1,1,1,0,1,-1,-1,
            -1,0,1,1,1,0,1,
            1,1,-1,0,0,-1,1,
            1,1,0,-1,-1,1,-1,
            0,-1,0,-1,1,-1,0,
            0,0,0,0,0,0,0,
            1,0,-1,-1,-1,0,-1,
            -1,1,0,1,-1,0,1,
            1,-1,-1,1,1,-1,1,
            1,-1,-1,1,0,1,1,
            0,1,0,-1,0,0,-1,
            0,1,0,1,-1,1,0,
            -1,1,1,-1,-1,1,-1,
            -1,-1,0,1,1,-1,1,
            -1,0,0,1,-1,-1,-1,
            -1,1,1,-1,0,-1,-1,
            -1,1,-1,0,1,1,-1,
            -1,-1,-1,0,-1,1,1,
            -1,-1,1,0,0,1,-1,
            -1,1,0,-1,1,-1,-1,
            1,-1,1,0,-1,-1,1),nrow=24,ncol=7,byrow=TRUE)


model<-"FullQuad"

k<-7
n<-24


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


##################################
#Effect generation information
##################################

m_main<-5
m_2FI<-4
m_quad<-3
eff_main_offset<-1.5
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


num_sims<-100


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








alphaME=0.05
pool="ME"
heredity="strong"
maxpower=FALSE
if(maxpower==TRUE){
  alpha2nd=0.05
}else{
  alpha2nd=0.20
}



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





#TPRF
mean(ME.TPR)
#FPRF
mean(ME.FPR)
#Fhat = F
mean(apply(AS.all.model.terms[,c(1:k)]==model_terms[,c(1:k)],1,sum)==6)
#TP2FI
mean(apply(model_terms[,c(1:choose(k,2)+k)]*AS.all.model.terms[,c(1:choose(k,2))+k],1,sum)/m_2FI)
#FP2FI
mean(apply((1-model_terms[,c(1:choose(k,2)+k)])*AS.all.model.terms[,c(1:choose(k,2))+k],1,sum)/(choose(k,2)-m_2FI))
#TPQ
mean(apply(model_terms[,-c(1:(choose(k,2)+k))]*AS.all.model.terms[,-c(1:(choose(k,2)+k))],1,sum)/m_quad)
#FPQ
mean(apply((1-model_terms[,-c(1:(choose(k,2)+k))])*AS.all.model.terms[,-c(1:(choose(k,2)+k))],1,sum)/(k-m_quad))
# Ahat = A
mean(apply(AS.all.model.terms==model_terms,1,sum)==ncol(X))
# Model Size
mean(apply(AS.all.model.terms,1,sum))

