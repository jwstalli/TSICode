#Source necessary R functions 
path="" #File path to "Function Library.R"
source(path)


##################################
thisn=12
thisk=5
model="ME2FI"  #ME, ME2FI, FullQuad
n.init<-2000
nrep=1
l=0
alpha=0.10
TSI.pure=TRUE
if(l>0){
  TSI.pure=FALSE
}

tau2=20

#Optional
StoreAll<-TRUE
SN<-1.25


CoorCand<-vector(mode="list",length=thisk)
for(j in 1:thisk){
  CoorCand[[j]]<-c(-1,1)
}



#Gather all designs that satisfy above and choose among those the one with smallest potential SEs


set.seed(1234)
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


hist(Try$Store.Measure.Primary/thisk,breaks=100)
DesignEval(D=Try$Dbest,model)




minSN<-0 #normally 0
maxSN<-1  #normally 1
Candidates<-which(Try$Store.Measure.Primary/thisk <= maxSN & Try$Store.Measure.Primary/thisk >= minSN)
length(Candidates)


n.so<-2
ERSS.measures<-rep(NA,length(Candidates))
#p<-thisk+choose(thisk,2)
p<-choose(thisk,2)
model.samp<-combinations(1:p,n.so)  
if(nrow(model.samp)>5000){
  model.samp<-model.samp[sample(1:nrow(model.samp),5000),]
}
for(i in 1:length(Candidates)){
  ERSS.measures[i]<-Model.ERSS(D=Try$Store.D[[Candidates[i]]],model,n.so,tau2=1,model.samp=model.samp)
  print(i)
}

#Optimal LOF value
ERSS.measures[Candidates[Opt]]

#Best design
Dbest<-Try$Store.D[[Candidates[Opt]]]

#Display design and properties
Dbest
DesignEval(Dbest,model)

