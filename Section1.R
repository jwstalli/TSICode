#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)


#NRFFD
Reactor.NRFFD<-matrix(c(-1,1,-1,1,1,
                   -1,1,1,1,-1,
                   -1,-1,1,-1,1,
                   1,1,1,1,1,
                   1,1,1,-1,-1,
                   1,-1,-1,1,-1,
                   -1,1,-1,-1,1,
                   1,-1,1,1,1,
                   1,1,-1,-1,-1,
                   -1,-1,1,-1,-1,
                   1,-1,-1,-1,1,
                   -1,-1,-1,1,-1),nrow=12,ncol=5,byrow=TRUE)

y.NRFFD<-c(78,95,59,82,61,61,70,42,61,53,63,69)

DesignEval(Reactor.NRFFD,model="ME2FI")

summary(lm(y.NRFFD~Reactor.NRFFD))



#Bayes D
Reactor.BayesD<-matrix(c(-1,-1,1,1,-1,
                   -1,-1,-1,-1,-1,
                   1,-1,-1,1,-1,
                   -1,1,-1,-1,1,
                   1,-1,-1,-1,1,
                   1,-1,1,1,1,
                   -1,1,-1,1,-1,
                   -1,-1,-1,1,1,
                   -1,-1,1,-1,1,
                   -1,1,1,1,1,
                   1,1,1,-1,-1,
                   1,1,-1,1,1),nrow=12,ncol=5,byrow=TRUE)

y.BayesD<-c(66,61,61,70,63,42,94,44,59,81,61,77)

DesignEval(Reactor.BayesD,model="ME2FI")

summary(lm(y.BayesD~Reactor.BayesD))


#EDMA

Reactor.EDMA<-matrix(c(1,-1,-1,1,-1,
                   -1,1,1,-1,1,
                   1,1,1,1,-1,
                   -1,-1,-1,-1,1,
                   1,    1,    1,   -1,   -1,
                   -1,-1,-1,1,1,
                   1,-1,-1,-1,1,
                   -1,1,1,1,-1,
                   1,1,-1,1,1,
                   -1,-1,1,-1,-1,
                   -1,1,-1,-1,-1,
                   1,-1,1,1,1),nrow=12,ncol=5,byrow=TRUE)

y.EDMA<-c(61,67,98,56,61,44,63,95,77,53,63,42)

DesignEval(Reactor.EDMA,model="ME2FI")

summary(lm(y.EDMA~Reactor.EDMA))






#ECI optimal design

Reactor.ECI<-matrix(c( 1,   -1,   -1,   -1,    1,
                  -1,   -1,   -1,   -1,   -1,
                  -1,   -1,   -1,    1,    1,
                   1,   -1,    1,    1,   -1,
                  -1,    1,   -1,   -1,    1,
                  -1,    1,    1,    1,   -1,
                   1,    1,    1,   -1,   -1,
                   1,    1,    1,    1,    1,
                  -1,   -1,    1,   -1,    1,
                   1,    1,   -1,    1,   -1,
                   1,   -1,   -1,   -1,    1,
                  -1,    1,    1,    1,   -1),nrow=12,ncol=5,byrow=TRUE)

y.ECI<-c(63,61,44,60,70,95,61,82,59,93) #Only partial responses


#Plugging in potential values using projection
# for run 11 we can try values 55, 56, 59
# for run 12 we can try values 93, 94, 98

DesignEval(Reactor.ECI,model="ME2FI")

rep.y<-expand.grid(c(55,56,59),c(93,94,98))
pvals<-matrix(NA,nrow=9,ncol=5)
est_model_terms<-matrix(NA,nrow=9,ncol=5+10)
sig2s<-rep(NA,9)
est.vals<-pvals

design.std.errors<-sqrt(diag(solve(crossprod(cbind(1,Reactor)))))[-1]

for(i in 1:9){
  
  y<-c(y.ECI,rep.y[i,1],rep.y[i,2])
  
  #Pure Error Estimate
  sig2<-1/2*(var(c(y[1],y[11]))+var(c(y[6],y[12])))
  sig2s[i]<-sig2
  
  est<-coef(lm(y~Reactor))[-1]
  est.vals[i,]<-est
  tvals<-est/(sqrt(sig2)*design.std.errors)
  
  pvals[i,]<-2*pt(abs(tvals),2,lower.tail=FALSE)
  
  Analysis<-GuidedSubsetsAnalysis(y=y,D=Reactor,model="2FI",heredity="strong",alphaME=0.10,alpha2nd=0.20,pool="all",
                                  maxpower=FALSE)
  
  est_model_terms[i,]<-Analysis$est_model_terms
  
}


design.std.errors
colMeans(est.vals)

boxplot(pvals,xlab="Factor",ylab="P-value",lwd=2)
abline(h=0.10,lwd=4,col=2)
