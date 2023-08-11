#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)

#Simulation Study

Reactor<-matrix(c(1,   -1,   -1,    -1,    1,
                  -1,   -1,   -1,   -1,   -1,
                  -1,   -1,   -1,   1,    1,
                  1,   -1,    1,    1,   -1,
                  -1,    1,   -1,   -1,    1,
                  -1,    1,    1,   1,   -1,
                  1,    1,    1,    -1,   -1,
                  1,    1,    1,    1,    1,
                  -1,   -1,    1,   -1,    1,
                  1,    1,   -1,    1,   -1,
                  1,   -1,   -1,    -1,    1,
                  -1,    1,    1,   1,   -1),nrow=12,ncol=5,byrow=TRUE)

y<-c(63,61,44,60,70,95,61,82,59,93)

Measure(Reactor,model="ME2FI",alpha=0.10,tau2=20,l=0,TSI.pure=TRUE)$Primary/5

Model.ERSS(D=Reactor,model="ME2FI",n.so=1,tau2=1)


#Plugging in potential values using projection
# for run 11 we can try values 55, 56, 59
# for run 12 we can try values 93, 94, 98


rep.y<-expand.grid(c(55,56,59),c(93,94,98))
pvals<-matrix(NA,nrow=9,ncol=5)
est_model_terms<-matrix(NA,nrow=9,ncol=5+10)
sig2s<-rep(NA,9)
est.vals<-pvals

design.std.errors<-sqrt(diag(solve(crossprod(cbind(1,Reactor)))))[-1]

for(i in 1:9){
  
  y<-c(63,61,44,60,70,95,61,82,59,93)
  y<-c(y,rep.y[i,1],rep.y[i,2])
  
  sig2<-1/2*(var(c(y[1],y[11]))+var(c(y[6],y[12])))
  sig2s[i]<-sig2
  
  est<-coef(lm(y~Reactor))[-1]
  est.vals[i,]<-est
  tvals<-est/(sqrt(sig2)*design.std.errors)
  
  pvals[i,]<-2*pt(abs(tvals),2,lower.tail=FALSE)
  
  Analysis<-AllSubsetsAnalysis(y=y,D=Reactor,model="2FI",heredity="strong",alphaME=0.10,pool="ME")
  
  est_model_terms[i,]<-Analysis$est_model_terms
  
}


#TPR_F
mean(apply(pvals[,c(2,4,5)]<0.1,1,function(x){sum(x)/3}))
#FPR_F
mean(apply(pvals[,c(1,3)]<0.1,1,function(x){sum(x)/2}))
#Fhat = F
mean(apply(pvals<0.1,1,function(x){sum(x== c(0,1,0,1,1))})==5)
#TPR_2FI
mean(apply(est_model_terms[,c(11,15)],1,function(x){sum(x)/2}))
#FPR_2FI
mean(apply(est_model_terms[,-c(1:5,11,15)],1,function(x){sum(x)/8}))
#Ahat = A
mean(apply(est_model_terms,1,function(x){
  sum(x==c(0,1,0,1,1,0,0,0,0,0,1,0,0,0,1))
})==15)
#|Ahat|
mean(apply(est_model_terms,1,sum))



###################
#Sim Reactor
###################

nsim<-100
pvals<-matrix(NA,nrow=nsim,ncol=5)
est_model_terms<-matrix(NA,nrow=nsim,ncol=5+10)


set.seed(1234)
for(i in 1:nsim){
  Ftrue<-sample(1:5)
  Dmix<-Reactor[,Ftrue]
  
    
  y<-65.5 + Dmix[,c(2,4,5)]%*%c(9.75,5.375,-3.125) + 6.625*Dmix[,2]*Dmix[,4]-5.5*Dmix[,4]*Dmix[,5]+rnorm(12,0,sd=3.331)
  
  
  sig2<-1/2*(var(c(y[1],y[11]))+var(c(y[6],y[12])))
  
  est<-coef(lm(y~Dmix))[-1]
  tvals<-est/(sqrt(sig2)*design.std.errors[Ftrue])
  
  pvals[i,]<-2*pt(abs(tvals),2,lower.tail=FALSE)
  
  Analysis<-AllSubsetsAnalysis(y=y,D=Dmix,model="2FI",heredity="strong",alphaME=0.10,pool="ME")
  
  est_model_terms[i,]<-Analysis$est_model_terms
  
}


#TPR_F
mean(apply(pvals[,c(2,4,5)]<0.1,1,function(x){sum(x)/3}))
#FPR_F
mean(apply(pvals[,c(1,3)]<0.1,1,function(x){sum(x)/2}))
#Fhat = F
mean(apply(pvals<0.1,1,function(x){sum(x== c(0,1,0,1,1))})==5)
#TPR_2FI
mean(apply(est_model_terms[,c(11,15)],1,function(x){sum(x)/2}))
#FPR_2FI
mean(apply(est_model_terms[,-c(1:5,11,15)],1,function(x){sum(x)/8}))
#Ahat = A
mean(apply(est_model_terms,1,function(x){
  sum(x==c(0,1,0,1,1,0,0,0,0,0,1,0,0,0,1))
})==15)
#|Ahat|
mean(apply(est_model_terms,1,sum))




###################
#Sim Reactor 2
###################

nsim<-100
pvals<-matrix(NA,nrow=nsim,ncol=5)
est_model_terms<-matrix(NA,nrow=nsim,ncol=5+10)


set.seed(1234)
for(i in 1:nsim){
  Ftrue<-sample(1:5)
  Dmix<-Reactor[,Ftrue]
  
  y<-65.5 + Dmix[,c(2,4,5)]%*%c(9.75,5.375,-3.125) + 6.625*Dmix[,2]*Dmix[,4]-5.5*Dmix[,4]*Dmix[,5]+10*Dmix[,2]*Dmix[,5]+rnorm(12,0,sd=3.331)
  
  
  sig2<-1/2*(var(c(y[1],y[11]))+var(c(y[6],y[12])))
  
  est<-coef(lm(y~Dmix))[-1]
  tvals<-est/(sqrt(sig2)*design.std.errors[Ftrue])
  
  pvals[i,]<-2*pt(abs(tvals),2,lower.tail=FALSE)
  
  Analysis<-AllSubsetsAnalysis(y=y,D=Dmix,model="2FI",heredity="strong",alphaME=0.10,pool="ME")
  
  est_model_terms[i,]<-Analysis$est_model_terms
  
}


#TPR_F
mean(apply(pvals[,c(2,4,5)]<0.1,1,function(x){sum(x)/3}))
#FPR_F
mean(apply(pvals[,c(1,3)]<0.1,1,function(x){sum(x)/2}))
#Fhat = F
mean(apply(pvals<0.1,1,function(x){sum(x== c(0,1,0,1,1))})==5)
#TPR_2FI
mean(apply(est_model_terms[,c(11,15)],1,function(x){sum(x)/2}))
#FPR_2FI
mean(apply(est_model_terms[,-c(1:5,11,15)],1,function(x){sum(x)/8}))
#Ahat = A for scenario 2
mean(apply(est_model_terms,1,function(x){
  sum(x==c(0,1,0,1,1,0,0,0,0,0,1,1,0,0,1))
})==15)
#|Ahat|
mean(apply(est_model_terms,1,sum))


