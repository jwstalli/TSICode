#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)


##################################
thisn=24
thisk=7
model="FullQuad"  #ME, ME2FI, FullQuad
n.init<-2000
nrep=0
l=1
alpha=0.05
TSI.pure=TRUE
if(l>0){
  TSI.pure=FALSE
}

tau2=20

#Optional
StoreAll<-TRUE
SN<-1


CoorCand<-vector(mode="list",length=thisk)
for(j in 1:thisk){
  #CoorCand[[j]]<-c(-1,1)
  CoorCand[[j]]<-c(-1,0,1)
  #CoorCand[[j]]<-seq(-1,1,by=0.01)
}



#Gather all designs that satisfy above and choose among those the one with smallest potential SEs

set.seed(1234) #This is just one example, results from all examples saved
Try<-TSI.Construct.r.l(thisn,
                       thisk,
                       model,
                       tau2,
                       n.init,
                       nrep,
                       l,
                       alpha,
                       CoorCand,
                       TSI.pure,
                       StoreAll,
                       SN)

#Combine Results
filepath="/Volumes/jwstalli/Documents/Research/Publications/2023/Two-Stage Inference Screening/Section 5/n24k7SimStudy/n24.k7.FullQuad.Rdata"
load(file=filepath)


minSN<-0 #normally 0
maxSN<-5  #normally 1
alpha=0.05
TSI.pure=TRUE
tau2<-20

AllDesigns<-c(AllTry$Try.1.1$Store.D,
              AllTry$Try.0.1$Store.D,
              AllTry$Try.2.0$Store.D)
Primary.Measures<-unlist(lapply(AllDesigns,function(x){
  Measure(x,model="FullQuad",alpha=alpha,l=0,TSI.pure=TSI.pure,tau2=tau2)$Primary
}))


Candidates<-which(Primary.Measures/thisk <= maxSN & Primary.Measures/thisk >= minSN)
length(Candidates)


hist(Primary.Measures[Candidates]/thisk,breaks=100,xlab="Expected CI MOE",main="")



minSN<-0 #normally 0
maxSN<-1  #normally 1
Candidates<-which(Primary.Measures/thisk <= maxSN & Primary.Measures/thisk >= minSN)
length(Candidates)


n.so<-6
ERSS.measures<-rep(NA,length(Candidates))
p<-choose(thisk,2)+thisk
model.samp<-combinations(1:p,n.so)  
if(nrow(model.samp)>5000){
  model.samp<-model.samp[sample(1:nrow(model.samp),5000),]
}
for(i in 1:length(Candidates)){
  ERSS.measures[i]<-Model.ERSS(D=AllDesigns[[Candidates[i]]],model,n.so,tau2=1,model.samp=model.samp)
  if(i%%50==0)print(round(i/length(Candidates),3))
}



AugDSD<-matrix(c(0,1,1,1,1,1,1,
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

Primary.AugDSD<-Measure(D=AugDSD,model="FullQuad",alpha=alpha,tau2=tau2,l=0,TSI.pure=TSI.pure)$Primary
Primary.Measures<-c(Primary.Measures[Candidates],Primary.AugDSD)/thisk
ERSS.measures<-c(ERSS.measures,Model.ERSS(AugDSD,model,n.so,tau2=1,model.samp=model.samp))




Utopia<-c(min(Primary.Measures),max(ERSS.measures))

distances<-abs(Primary.Measures-Utopia[1])+abs(ERSS.measures-Utopia[2])

plot(x=Primary.Measures,y=ERSS.measures,pch=16,col="grey",xlab="Expected CI MOE",ylab="Min Expected LOF")


points(x=Utopia[1],y=Utopia[2],pch=16,col=2,cex=3)

Opt.val<-max(round(ERSS.measures[-length(ERSS.measures)],12))
Opt.all<-which(round(ERSS.measures,12)==Opt.val)
Opt<-Opt.all[which.min(Primary.Measures[Opt.all])]

points(x=Primary.Measures[Opt],y=ERSS.measures[Opt],pch=15,col=5,cex=3)

points(x=Primary.Measures[length(Candidates)+1],
       y=ERSS.measures[length(Candidates)+1],pch=18,col=4,cex=3)


Dbest<-AllDesigns[[Candidates[Opt]]]

Dbest

DesignEval(Dbest,model)
DesignEval(AugDSD,model)


##########################################
# Extra: Leonard Edwards Design Evaluation
##########################################

LeonardD<-matrix(c(-1,-1,-1,1,-1,1,1,
                   -1,-1,1,0,-1,-1,0,
                   -1,0,-1,-1,1,0,0,
                   -1,0,-1,-1,1,0,0,
                   -1,0,0,0,0,0,-1,
                   -1,1,-1,1,1,1,-1,
                   -1,1,1,-1,-1,1,1,
                   -1,1,1,1,1,-1,1,
                   0,-1,0,-1,0,-1,1,
                   0,0,0,0,1,1,0,
                   0,0,0,0,1,1,0,
                   0,0,1,1,-1,0,-1,
                   0,0,1,1,-1,0,-1,
                   0,1,-1,0,0,0,0,
                   1,-1,-1,-1,-1,1,-1,
                   1,-1,-1,1,1,-1,-1,
                   1,-1,1,0,1,0,1,
                   1,-1,1,0,1,0,1,
                   1,0,-1,0,-1,-1,1,
                   1,0,1,1,0,1,0,
                   1,0,1,1,0,1,0,
                   1,1,0,1,-1,0,0,
                   1,1,0,1,-1,0,0,
                   1,1,1,-1,1,-1,-1),nrow=24,ncol=7,byrow=TRUE)
DesignEval(LeonardD,model)


