library(Matrix)
library(MASS)
library(arrangements)
library(lmSubsets) #for fast all subsets



#############################################################
#############################################################
#############################################################

# Preliminary Functions for Design Construction and Analysis

#############################################################
#############################################################
#############################################################




####################
#Projection Matrix
####################
Proj<-function(X){
  P<-X%*%ginv(crossprod(X))%*%t(X)
  return(P)
}



##############################################################
#Transform D to X for given model
##############################################################

#Function to generate X matrix
# thisD = Design matrix
# model = "ME", "ME2FI", "FullQuad"

DtoX <- function(thisD,model="ME2FI"){
  
  thisn <- nrow(thisD)
  thisk <- ncol(thisD)
  
  if(model=="ME"){
    
    X<-cbind(1,thisD)
    
  }else if(model=="ME2FI"){
    
    X.2FI<-matrix(NA,nrow=thisn,ncol=choose(thisk,2))
    p<-1
    for(i in 1:(thisk-1)){
      for(ii in (i+1):thisk){
        X.2FI[,p]<-thisD[,i]*thisD[,ii]
        p        <-p+1
      }
    }
    
    X<-cbind(1,thisD,X.2FI)
    
  }else if(model=="FullQuad"){
    X.2FI<-matrix(NA,nrow=thisn,ncol=choose(thisk,2))
    p<-1
    for(i in 1:(thisk-1)){
      for(ii in (i+1):thisk){
        X.2FI[,p]<-thisD[,i]*thisD[,ii]
        p        <-p+1
      }
    }
    
    X.Q<-thisD^2
    
    X<-cbind(1,thisD,X.2FI,X.Q)
    
    
  }else print(paste("Model ",model," Not Recognized.")); stop;  
  
  return(X)
}











#################################################
### Evaluate properties of design and print out
#################################################

DesignEval<-function(D,model){
  n=nrow(D)
  k=ncol(D)
  
  X<-DtoX(D,model)
  
  X1<-X[,1:(k+1)]
  X2<-X[,-c(1:(k+1))]
  
  MEvars<-diag(solve(crossprod(X1)))[-1]
  
  print("ME model design Standard Errors")
  print(sqrt(MEvars))
  
  Alias<-(solve(crossprod(X1))%*%t(X1)%*%X2)[-1,]
  
  print("Aliasing")
  print(summary(abs(c(Alias))))
  
  d<-n-qr(crossprod(DtoX(D,model)))$rank
  PureDF<-sum(unlist(lapply(FindReps(D),length))-1)
  LoFDF<-d-PureDF
  print("Variance Estimation")
  print(paste("d= ",d,"   PureDF= ",PureDF,"    LOFDF= ",LoFDF,sep=""))
  
  
}







##############################################################
#Initial Design Generates Random Design According to CoorCand
##############################################################
#thisn    <- n total runs
#thisk    <- k factors
#nrep     <- number of required replicated rows
#CoorCand <- list of candidate values

GenInit<-function(thisn,thisk,nrep=0,CoorCand){
  if(nrep==0){
    Du<-matrix(NA,nrow=thisn,ncol=thisk)
    for(j in 1:thisk){
      Du[,j]<-sample(CoorCand[[j]],thisn,replace=TRUE)
    }
    D<-Du
    Dr<-NULL
    Dr.indx=NULL
    reps=NULL
  }else{
    thisn.u<-thisn-nrep
    Du<-matrix(NA,nrow=thisn.u,ncol=thisk)
    for(j in 1:thisk){
      Du[,j]<-sample(CoorCand[[j]],thisn.u,replace=TRUE)
    }
    reps<-sort(sample(1:thisn.u,nrep,replace=TRUE))
    Dr<-Du[reps,]
    D<-rbind(Du,Dr)
    Dr.indx=c(1:nrep)+thisn.u
  }
  return(list(D=D,Du=Du,Dr=Dr,Dr.indx=Dr.indx,whichrep=reps))
  
  #Output
  #D         = matrix with Du and Dr combined
  #Du        = matrix of unrestricted design points
  #Dr        = matrix of replicated rows corresponding to rows in Du
  #Dr.indx   = which rows in D correspond to Dr
  #whichreps = vector connecting rows in Dr.indx and Du
  
}










##############################################################
#Find replicated rows for pure error
##############################################################

FindReps<-function(thisD){
  runs<-which(duplicated(thisD))
  D.unique<-thisD[-runs,]
  D.rep<-thisD[runs,]
  if(length(runs)==0){
    RepVec<-lapply(1:nrow(thisD),function(x){x})
  }else{
    if(length(runs)==1){
      D.rep<-t(as.matrix(D.rep))
    }
    RepVec<-vector(mode="list",length=dim(D.unique)[1])
    for(i in 1:length(RepVec)){
      check<-apply(D.rep,1,function(x){sum(abs(D.unique[i,]-x))})==0
      ireps<-c(i,which(check)+length(RepVec))
      names(ireps)<-NULL
      RepVec[[i]]<-ireps
    }
  }
  
  #D<-rbind(D.unique,D.rep)
  
  return(RepVec=RepVec)
  
  #Output: list with elements = # unique rows of D and their replicated rows
  #Note: Could include more rows than what is in Dr
  
}



















####################################################
####################################################
####################################################

# Analysis Code

####################################################
####################################################
####################################################



SequentialSS<-function(y,D,X,A.ME,PureMSE,PureDF,model="FullQuad",heredity="none"){
  
  n<-nrow(D)
  k<-ncol(D)
  
  
  #Active factor labels
  ActiveLabels<-which(A.ME==1)
  m.ME<-sum(A.ME)
  
  ################################################################################
  #Calculate SS for insig main effects and pool with earlier MSE
  ################################################################################
  SSE.full<-anova(lm(y~D))$"Sum Sq"[2]
  SSE.red<-anova(lm(y~D[,ActiveLabels]))$"Sum Sq"[2]
  
  PureSSE<-PureMSE*PureDF
  PureSSE<-PureSSE+(SSE.red-SSE.full)
  PureDF<-PureDF+(k-m.ME)
  
  PureMSE<-PureSSE/PureDF
  
  
  X1<-cbind(1,D)
  P1<-Proj(X1)
  
  
  if(heredity=="none"){
    
    if(model=="2FI"){
      allowed_terms=c(1:choose(k,2))+k
      excluded_terms=c()
    }
    
    if(model=="FullQuad"){
      allowed_terms=c(1:(choose(k,2)+k)+k)
      excluded_terms=c()
    }
    
  }
  
  
  
  if(heredity=="strong"){
    
    #2FI model, Need more than 1 active factor for possible interaction
    if(model=="2FI"){
      
      
      #Interaction labels with 3 columns: columns 1 and 2 give int pair, 
      #column 3 = index for "terms" object, 
      
      IntLabels<-cbind(combinations(1:k,2),
                       c(1:choose(k,2)+k)
      )
      
      
      if(sum(A.ME)==1){  
        
        allowed_2FIs_terms<-c()
        excluded_2FIs_terms<-c(1:choose(k,2))+k
        
      }
      
      if(sum(A.ME)>1){
        
        allowed_2FIs<-data.frame(combinations(ActiveLabels,2))
        
        #Terms index of overlap of allowed 2FIs and all possible 2FIs
        allowed_2FIs_terms<-c(merge(allowed_2FIs, data.frame(IntLabels))[,3])
        excluded_2FIs_terms<-setdiff((1:choose(k,2))+k,allowed_2FIs_terms)
        
      }
      
      allowed_terms=allowed_2FIs_terms
      excluded_terms=excluded_2FIs_terms
      
      
      
    }
    
    
    if(model=="FullQuad"){
      
      
      if(sum(A.ME)==1){  #Only add quadratic effect for one sig main effect
        
        allowed_terms<-ActiveLabels+k+choose(k,2)
        excluded_terms<-c( c(1:(choose(k,2)+k))+k )[-c(allowed_terms-k)]
        
      }
      
      if(sum(A.ME)>1){
        
        #Interaction labels with 3 columns: columns 1 and 2 give int pair, 
        #column 3 = index for "terms" object, 
        
        IntLabels<-cbind(combinations(1:k,2),
                         c(1:choose(k,2)+k)
        )
        
        #Quadratic labels just 1 column, same order as main effects
        QuadLabels<-c(1:k)+k+choose(k,2)
        
        
        allowed_2FIs<-data.frame(combinations(ActiveLabels,2))
        
        #Terms index of overlap of allowed 2FIs and all possible 2FIs
        allowed_2FIs_terms<-c(merge(allowed_2FIs, data.frame(IntLabels))[,3])
        excluded_2FIs_terms<-setdiff((1:choose(k,2))+k,allowed_2FIs_terms)
        
        allowed_Quad_terms<-QuadLabels[ActiveLabels]
        excluded_Quad_terms<-QuadLabels[which(A.ME==0)]
        
        
        allowed_terms<-c(allowed_2FIs_terms,allowed_Quad_terms)
        excluded_terms<-c(excluded_2FIs_terms,excluded_Quad_terms)
        
        
      }
      
      
    }
    
  } #End heredity=strong
  
  
  
  
  if(heredity=="weak"){
    
    #2FI model, Need more than 1 active factor for possible interaction
    if(model=="2FI"){
      
      
      #Interaction labels with 3 columns: columns 1 and 2 give int pair, 
      #column 3 = index for "terms" object, 
      
      IntLabels<-cbind(combinations(1:k,2),
                       c(1:choose(k,2)+k)
      )
      
      allowed_2FIs<-apply(IntLabels[,-3],1,function(x){sum(ActiveLabels %in% x)>0})
      allowed_2FIs_terms <-IntLabels[ which(allowed_2FIs==1) ,3]
      excluded_2FIs_terms<-IntLabels[ which(allowed_2FIs==0) ,3]
      
      
      allowed_terms<-allowed_2FIs_terms
      excluded_terms<-excluded_2FIs_terms
      
    }
    
    
    if(model=="FullQuad"){
      
      #Interaction labels with 3 columns: columns 1 and 2 give int pair, 
      #column 3 = index for "terms" object, 
      
      IntLabels<-cbind(combinations(1:k,2),
                       c(1:choose(k,2)+k)
      )
      
      #Quadratic labels just 1 column, same order as main effects
      QuadLabels<-c(1:k)+k+choose(k,2)
      
      
      allowed_2FIs<-apply(IntLabels[,-3],1,function(x){sum(ActiveLabels %in% x)>0})
      allowed_2FIs_terms <-IntLabels[ which(allowed_2FIs==1) ,3]
      excluded_2FIs_terms<-IntLabels[ which(allowed_2FIs==0) ,3]
      
      
      allowed_Quad_terms<-QuadLabels[ActiveLabels]
      excluded_Quad_terms<-QuadLabels[which(A.ME==0)]
      
      
      allowed_terms<-c(allowed_2FIs_terms,allowed_Quad_terms)
      excluded_terms<-c(excluded_2FIs_terms,excluded_Quad_terms)
      
      
    }
    
    
    
  }#End heredity=weak
  
  
  
  if(is.null(allowed_terms)){
    X2.fit<-matrix(0,nrow=n,ncol=1)
  }else{
    X2.fit<-X[,allowed_terms]
  }
  
  if(is.null(excluded_terms)){
    X2.exclude<-matrix(0,nrow=n,ncol=1)
  }else{
    X2.exclude<-X[,excluded_terms]  
  }
  
  
  Xfull<-cbind(1,X)
  PX<-Proj(Xfull)
  
  X.fit<-cbind(1,D,X2.fit)
  PX.fit<-Proj(X.fit)
  
  X2adj.fit<-(PX.fit-P1)%*%X2.fit
  X2adj.fit[which(abs(X2adj.fit)<1e-13)]=0
  if(ncol(X2adj.fit)==0){
    r2.fit=0
  }else{
    r2.fit<-qr(X2adj.fit)$rank
  }
  yadj.fit<-(PX-P1)%*%y
  
  
  
  X2adj.null<-(PX-PX.fit)%*%X2.exclude
  X2adj.null[which(abs(X2adj.null)<1e-13)]=0
  if(ncol(X2adj.null)==0){
    r2.null=0
  }else{
    r2.null<-qr(X2adj.null)$rank
  }
  yadj.null<-(PX-PX.fit)%*%y
  
  
  return(list(PureMSE=PureMSE,PureDF=PureDF,
              allowed_terms=allowed_terms,
              excluded_terms=excluded_terms,
              X2adj.fit=X2adj.fit,
              r2.fit=r2.fit,
              yadj.fit=yadj.fit,
              X2adj.null=X2adj.null,
              r2.null=r2.null,
              yadj.null=yadj.null))
  
}





#pool = none, ME, all
GuidedSubsetsAnalysis<-function(y,D,model="FullQuad",heredity="none",alphaME=0.05,alpha2nd=0.2,
                                pool="none",maxpower=TRUE,maxmod="0.5"){
  
  
  n<-nrow(D)
  k<-ncol(D)
  
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
  
  
  
  
  #############################
  #Generate error estimate (we require either pure error of LOF)
  #############################
  
  #Grab pure+FF MSE from full model analysis in lm
  fullfit<-lm(y~X)
  PureMSE<-summary(fullfit)$sigma^2
  PureDF<-summary(fullfit)$df[2]
  
  
  #############################
  #Main effects only analysis
  #############################
  
  #Use lm to fit main effect model and extract estimates
  MEfit<-lm(y~D)
  MEest<-MEfit$coefficients[-1]
  
  ME.design.standard.errors<-sqrt(diag(solve(crossprod(cbind(1,D))))[-1])
  
  #Calculate confidence intervals
  tcritME<-abs(qt(alphaME/2,PureDF))
  lowCI<-MEest-sqrt(PureMSE)*ME.design.standard.errors*tcritME
  UppCI<-MEest+sqrt(PureMSE)*ME.design.standard.errors*tcritME
  CIs<-cbind(lowCI,UppCI)
  
  
  #Check for containing 0 or not
  A.ME<-1-apply(CIs,1,function(x){x[1] < 0 & 0 < x[2]})
  m.ME<-sum(A.ME)
  
  
  #############################
  #Stage 2 Analysis
  #############################        
  
  
  if(m.ME > 0){
    
    
    #Calculate sequential sum-of-squares based on heredity
    Stage2<-SequentialSS(y=y,D=D,X=X,A.ME=A.ME,PureMSE=PureMSE,PureDF=PureDF,model=model,heredity=heredity)
    
    #If we allow main effect pooling
    if(pool=="ME" | pool=="all"){
      
      PureMSE<-Stage2$PureMSE
      PureDF<-Stage2$PureDF
      
    }
    
    #Further pooling of null terms, if allowed
    if(pool=="all"){
      
      PureMSE<-(PureMSE*PureDF+c(crossprod(Stage2$yadj.null)))/(PureDF+Stage2$r2.null)
      PureDF<-PureDF+Stage2$r2.null  
      
    }
    
    
    
    
    
    #Calculate F-value for largest submodel
    if(Stage2$r2.fit>0){
      Global.F<-c((crossprod(Stage2$yadj.fit)/Stage2$r2.fit)/PureMSE)
      Global.p<-pf(Global.F,Stage2$r2.fit,PureDF,lower.tail=FALSE)  
    }else{
      Global.F<-0
      Global.p<-1
    }
    
    
    
    
    if(Global.p < alpha2nd){
      
      if(Stage2$r2.fit==1){
        best.mod<-c(1:ncol(Stage2$X2adj.fit))
      }else{
        if(maxmod=="0.5"){
          max.size<-min(Stage2$r2.fit,floor((Stage2$r2.fit+Stage2$r2.null)/2))
        }else{
          max.size<-min(Stage2$r2.fit,floor((Stage2$r2.fit+Stage2$r2.null)-1))
        }
        
        mod.size<-1
        LoFcheck<-0
        
        n.terms<-ncol(Stage2$X2adj.fit)
        
        
        #######################
        ###### OLD STUFF #####
        #######################
        
        while(mod.size <= max.size & LoFcheck==0){
          
          models<-combinations(1:n.terms,mod.size)
          
          LoFdf<-Stage2$r2.fit-mod.size
          LoF.pvals<-apply(models,1,function(x){
            fit<-lm(Stage2$yadj.fit~Stage2$X2adj.fit[,x])
            Fval<-(anova(fit)$"Sum Sq"[2]/LoFdf)/PureMSE
            pval<-pf(Fval,LoFdf,PureDF,lower.tail=FALSE)
          })
          
          
          best.mod.indx<-which.max(LoF.pvals)
          
          #If no LOF then compile all effects from models where LOF > alpha2nd
          if(LoF.pvals[best.mod.indx] > alpha2nd){
            
            LoFcheck<-1
            #  if(maxpower==TRUE){
            #    best.mod<-sort(unique(c(models[ which(LoF.pvals > alpha2nd), ])))  
            #  }else{
            best.mod<-sort((models[best.mod.indx, ]))
            #}
            
          }
          
          mod.size<-mod.size+1
        }
        
        #Still more important terms, just accept the rest and say inconclusive
        if(mod.size>max.size & LoFcheck==0){
          best.mod<-sort((models[best.mod.indx, ]))
        }
        #if(maxpower==TRUE){
        #    best.mod<-c(1:ncol(Stage2$X2adj.fit))
        #    print(paste(m.ME,"All Accept"))
        #}else{
        
        #}
        #}
        
      }
      
      #Use allowed_terms from Stage2 object to know what terms selected
      allowed_terms<-Stage2$allowed_terms
      
      est_model_terms<-rep(0,ncol(X))
      est_model_terms[1:k]<-A.ME
      est_model_terms[allowed_terms[best.mod]]<-1
      
      
    }else{
      est_model_terms<-rep(0,ncol(X))
      est_model_terms[1:k]<-A.ME
    } #End guided subset check
    
    
  }else{  #End check for at least one sig main effect
    
    Stage2=NULL
    est_model_terms<-rep(0,ncol(X))
    
  } 
  
  
  return(list(est_model_terms=est_model_terms,
              PureMSE=PureMSE,
              PureDF=PureDF,
              CIs=CIs,
              Stage2=Stage2))
  
}







#pool = none, ME, all
AllSubsetsAnalysis<-function(y,D,model="FullQuad",heredity="none",alphaME=0.05,alpha2nd=0.20,pool="none"){
  
  
  n<-nrow(D)
  k<-ncol(D)
  
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
  
  
  
  
  #############################
  #Generate error estimate (we require either pure error of LOF)
  #############################
  
  #Grab pure+FF MSE from full model analysis in lm
  fullfit<-lm(y~X)
  PureMSE<-summary(fullfit)$sigma^2
  PureDF<-summary(fullfit)$df[2]
  
  
  #############################
  #Main effects only analysis
  #############################
  
  #Use lm to fit main effect model and extract estimates
  MEfit<-lm(y~D)
  MEest<-MEfit$coefficients[-1]
  
  ME.design.standard.errors<-sqrt(diag(solve(crossprod(cbind(1,D))))[-1])
  
  #Calculate confidence intervals
  tcritME<-abs(qt(alphaME/2,PureDF))
  lowCI<-MEest-sqrt(PureMSE)*ME.design.standard.errors*tcritME
  UppCI<-MEest+sqrt(PureMSE)*ME.design.standard.errors*tcritME
  CIs<-cbind(lowCI,UppCI)
  
  
  #Check for containing 0 or not
  A.ME<-1-apply(CIs,1,function(x){x[1] < 0 & 0 < x[2]})
  m.ME<-sum(A.ME)
  
  
  #############################
  #Stage 2 Analysis
  #############################        
  
  
  if(m.ME > 0){
    
    
    #Calculate sequential sum-of-squares based on heredity
    Stage2<-SequentialSS(y=y,D=D,X=X,A.ME=A.ME,PureMSE=PureMSE,PureDF=PureDF,model=model,heredity=heredity)
    
    #If we allow main effect pooling
    if(pool=="ME" | pool=="all"){
      
      PureMSE<-Stage2$PureMSE
      PureDF<-Stage2$PureDF
      
    }
    
    #Further pooling of null terms, if allowed
    if(pool=="all"){
      
      PureMSE<-(PureMSE*PureDF+c(crossprod(Stage2$yadj.null)))/(PureDF+Stage2$r2.null)
      PureDF<-PureDF+Stage2$r2.null  
      
    }
    
    
    
    
    
    #Calculate F-value for largest submodel
    if(Stage2$r2.fit>0){
      Global.F<-c((crossprod(Stage2$yadj.fit)/Stage2$r2.fit)/PureMSE)
      Global.p<-pf(Global.F,Stage2$r2.fit,PureDF,lower.tail=FALSE)  
    }else{
      Global.F<-0
      Global.p<-1
    }
    
    
    
    
    if(Global.p < alpha2nd){
      
      if(Stage2$r2.fit==1){
        best.mod<-c(1:ncol(Stage2$X2adj.fit))
      }else{
        max.size<-min(Stage2$r2.fit,floor(Stage2$r2.fit+Stage2$r2.null-1))
        
        
        if(nrow(Stage2$X2adj.fit)<=ncol(Stage2$X2adj.fit)){
          fix.x<-rbind(Stage2$X2adj.fit,Stage2$X2adj.fit)
          fix.y<-c(Stage2$yadj.fit,Stage2$yadj.fit)
          fit<-lmSubsets_fit(x=fix.x,y=fix.y,
                             nbest=50000,intercept=FALSE)  
          fit$submodel$RSS<-fit$submodel$RSS/2
        }else{
          fit<-lmSubsets_fit(x=Stage2$X2adj.fit,y=Stage2$yadj.fit,
                             nbest=50000,intercept=FALSE,nmax=max.size)  
        }
        fit.subset<-which(is.na(fit$submodel$RSS) | fit$submodel$SIZE>max.size)
        fit.RSS<-fit$submodel[-fit.subset,]
        fit.RSS$BIC<-NA
        fit.RSS$good<-0
        fit.submodels<-fit$subset[-fit.subset,]
        
        fit.RSS$BIC<-fit.RSS$RSS/PureMSE+log(Stage2$r2.fit+Stage2$r2.null)*fit.RSS$SIZE
        
        best.mod<-which(as.matrix(fit.submodels[which.min(fit.RSS$BIC),]))
        
      }
      
      
      #Use allowed_terms from Stage2 object to know what terms selected
      allowed_terms<-Stage2$allowed_terms
      
      
      est_model_terms<-rep(0,ncol(X))
      est_model_terms[1:k]<-A.ME
      est_model_terms[allowed_terms[best.mod]]<-1
      
      
    }else{
      est_model_terms<-rep(0,ncol(X))
      est_model_terms[1:k]<-A.ME
    } #End guided subset check
    
    
  }else{  #End check for at least one sig main effect
    
    Stage2=NULL
    est_model_terms<-rep(0,ncol(X))
    
  } 
  
  
  
  return(list(est_model_terms=est_model_terms,
              PureMSE=PureMSE,
              PureDF=PureDF,
              CIs=CIs,
              Stage2=Stage2))
  
}























#############################################################
#############################################################
#############################################################

# Construction Code

#############################################################
#############################################################
#############################################################



########################################
# Calculate multiplier
########################################

#TSI.pure means we only do pure DF and so do the expected standard deviation formula
cmult<-function(X1,X2,alpha,tau2,l=0,TSI.pure=FALSE,A){
  
  D<-X1[,-1] #Assumes X1 is always intercept + ME's
  X<-cbind(X1,X2)
  n<-nrow(X)
  
  
  #########################
  ## Find d from null rank
  #########################
  d<-n-qr(crossprod(X))$rank
  
  
  ###############################
  ## Find pure error of pure LOF
  ###############################
  PureDF<-sum(unlist(lapply(FindReps(D),length))-1)
  LoFDF<-d-PureDF
  
  
  ##############################
  ## Find contaminated FF, if necessary
  ##############################
  tr.bias=0
  l.more=l-LoFDF #num of contaminated fake factors
  if(l.more > 0 & !TSI.pure){
    
    C21<-crossprod(X2 - tcrossprod(X1,t(A)))
    C21[which(abs(C21)<1e-13)]<-0
    
    rank.X21<-qr(C21)$rank
    C21.evals<-Re(eigen(C21)$values)
    
    tr.bias<-sum(C21.evals[c(rank.X21:(rank.X21-l.more+1))])
    
    d<-d+l.more
  }
  
  ###############################
  ## Calculate t value
  ###############################
  if(d==0){d=0.1} #Just in case! Will cause multiplier to blow up
  tval<-qt(alpha/2,d,lower.tail=FALSE)
  
  if(TSI.pure){
    cm<-tval*sqrt(2/d)*gamma((d+1)/2)/gamma(d/2)
  }else{
    cm<-tval*sqrt(1+tau2/d*tr.bias)
  }
  
  
  return(list(cm=cm,d=d,PureDF=PureDF,LoFDF=LoFDF,Contaminate=l.more))
  
}





Measure<-function(D,model,alpha,tau2,l,TSI.pure){
  
  thisn<-nrow(D)
  thisk<-ncol(D)
  
  X<-DtoX(D,model=model)
  
  X1<-X[,c(1:(thisk+1))]
  X1tX1inv<-chol2inv(chol(crossprod(X1)))
  X2<-X[,-c(1:(thisk+1))]
  
  
  A<-X1tX1inv%*%crossprod(X1,X2)
  Sqrt.Alias<-sqrt(apply(A[-1,]^2,1,sum)) #removes intercept
  
  RepVec<-FindReps(D)
  PureDF<-sum(unlist(lapply(RepVec,length))-1)
  
  rank.X<-qr(crossprod(X))$rank
  d<-thisn-rank.X
  
  cm<-cmult(X1,X2,alpha,tau2,l,TSI.pure,A)$cm
  
  Primary<-sqrt(2/pi)*sqrt(tau2)*sum(Sqrt.Alias)+
    cm*(  sum(sqrt(diag(X1tX1inv)[-1])))    
  
  return(list(Primary=Primary))
  
}




########################################
# Update matrix V with row exchange (works for both main effects and updating beta2 standard errors)
########################################

TSI.Update.V<-function(V,ri,xi,xtilde){
  
  F1<-cbind(xtilde,-1*xi)
  F2<-cbind(xtilde,xi)
  
  v.xi<-c(t(xi)%*%V%*%xi)
  v.xtilde<-c(t(xtilde)%*%V%*%xtilde)
  v.xi.xtilde<-c(t(xi)%*%V%*%xtilde)
  
  Det<-(1+ri*v.xtilde)*(1-ri*v.xi)+ri^2*v.xi.xtilde^2
  
  if(Det<1e-16) print("Singular Matrix")
  
  InvMat<-ri/Det*matrix(c(1-ri*v.xi, ri*v.xi.xtilde,-ri*v.xi.xtilde, 1+ri*v.xtilde),nrow=2,ncol=2,byrow=TRUE)  
  
  Vtilde<-V-crossprod(V,F1)%*%InvMat%*%crossprod(F2,V)
  
  return(list(Vtilde=Vtilde,Det=Det))
}




#Xr has columns of new x's
TSI.Update.V.reps<-function(Vu,Xr){
  
  F1<-Xr
  F2<-t(F1)
  
  nrep<-ncol(Xr)
  
  InvMat<-chol2inv(chol(diag(nrep)+F2%*%Vu%*%F1))
  
  Vtilde<-Vu-crossprod(Vu,F1)%*%InvMat%*%crossprod(F1,Vu)
  
  return(Vtilde)
}





#V here is the previous main effects model covariance matrix, A is the previous alias matrix
# rows and col say which coordinate or coordinates we are changing
TSI.CoorUpdate<-function(D,X,V,A,rows,col,model,alpha,tau2,l,CoorCand,Dr.indx,TSI.pure=FALSE){
  
  n<-nrow(D)
  k<-ncol(D)
  
  nexch<-length(CoorCand[[col]])
  k<-length(CoorCand)
  X1<-X[,1:(k+1)]
  X2<-X[,-c(1:(k+1))]
  
  
  
  ##############################
  # Candidate Designs
  ##############################
  Dcand<-lapply(CoorCand[[col]],
                function(x){
                  Dnew<-D
                  Dnew[rows,col]=x
                  return(Dnew)
                }
  )
  
  
  ##############################
  # Candidate X's
  ##############################
  Xcand<-lapply(Dcand,function(D){DtoX(D,model=model)})
  X1cand<-lapply(Xcand,function(X){X[,1:(k+1)]})
  X2cand<-lapply(Xcand,function(X){X[,-c(1:(k+1))]})
  
  
  ##############################
  # Candidate V1tilde
  ##############################
  
  xi=as.matrix(X1[rows[1],],nrow=ncol(X1),ncol=1)
  Xtilde=lapply(X1cand,function(x){as.matrix(x[rows[1],],nrow=ncol(X1),ncol=length(rows))})
  
  V1tilde<-lapply(Xtilde,function(xtilde){
    Update<-TSI.Update.V(V,
                         ri=length(rows),
                         xi=xi,
                         xtilde=xtilde)
    if(Update$Det<1e-16){
      Update<-FALSE
    }else{
      Update<-Update$Vtilde
    }
    return(Update)
  }
  )
  
  #Update if singular design
  CheckSingular<-unlist(lapply(V1tilde,is.matrix))
  drop.cand<-which(CheckSingular==FALSE)
  
  if(length(drop.cand)>0){
    modCoorCand<-CoorCand[[col]][-drop.cand]
    
    nexch<-length(modCoorCand)
    ##############################
    # Candidate Designs
    ##############################
    Dcand<-lapply(modCoorCand,
                  function(x){
                    Dnew<-D
                    Dnew[rows,col]=x
                    return(Dnew)
                  }
    )
    
    
    ##############################
    # Candidate X's
    ##############################
    Xcand<-lapply(Dcand,DtoX,model=model)
    X1cand<-lapply(Xcand,function(X){X[,1:(k+1)]})
    X2cand<-lapply(Xcand,function(X){X[,-c(1:(k+1))]})
    
    V1tilde<-V1tilde[-drop.cand]
  }
  
  
  ##############################
  # Candidate Alias Matrices
  ##############################
  
  Acand<-lapply(1:nexch,function(cand){
    if(sum(V1tilde[[cand]])==0){
      return(Inf)
    }else{
      A<-crossprod(V1tilde[[cand]],crossprod(X1cand[[cand]],X2cand[[cand]]))
      return(A) 
    }
  }
  )
  
  
  
  ##############################
  # Multiplier
  ##############################
  
  cmultiplier<-unlist(lapply(1:nexch,function(cand){
    if(sum(V1tilde[[cand]])==0){
      return(FALSE)
    }else{
      cmult(X1cand[[cand]],
            X2cand[[cand]],
            alpha=alpha,
            tau2=tau2,l=l,TSI.pure=TSI.pure,A=Acand[[cand]])$cm
    }
  }))
  
  
  ##############################
  # Final Measures
  ##############################
  
  Measures<-unlist(lapply(1:nexch,function(cand){
    sqrt(tau2*2/pi)*sum(sqrt(apply(Acand[[cand]][-1,]^2,1,sum)))+
      cmultiplier[cand]*(sum(sqrt(diag(V1tilde[[cand]]))[-1]))
  }
  )
  )
  
  BestExch<-which.min(Measures)
  
  ########
  #Output#
  ########
  D=Dcand[[BestExch]]
  X=Xcand[[BestExch]]
  V=V1tilde[[BestExch]]
  A=Acand[[BestExch]]
  
  
  return(list(D=D,X=X,V=V,A=A,Primary=Measures[BestExch]))
  
}



########################################
# Stage 2: Row Updates
########################################
TSI.RowUpdate<-function(Du,Xu,Vu,rows,model,alpha,tau2,l,TSI.pure=FALSE){
  
  k<-ncol(Du)
  X1u<-Xu[,1:(k+1)]
  X2u<-Xu[,-c(1:(k+1))]
  
  if(length(rows)==1){
    Dr<-t(as.matrix(Du[rows,]))
  }else{
    Dr<-as.matrix(Du[rows,])
  }
  
  
  ##############################
  # Candidate Design
  ##############################
  Dcand<-rbind(Du,Dr)
  Dr.indx<-c((nrow(Du)+1):nrow(Dcand))
  
  
  ##############################
  # Candidate X's
  ##############################
  Xcand<-DtoX(Dcand,model=model)
  X1cand<-Xcand[,1:(k+1)]
  X2cand<-Xcand[,-c(1:(k+1))]
  
  ##############################
  # V1tilde
  ##############################
  
  #Xr must be column vectors
  if(length(rows)==1){
    Xr<-as.matrix(X1cand[Dr.indx,])
  }else{
    Xr=t(as.matrix(X1cand[Dr.indx,]))  
  }
  
  V1tilde<-TSI.Update.V.reps(Vu,Xr)
  
  ##############################
  # Alias
  ##############################
  A<-crossprod(V1tilde,crossprod(X1cand,X2cand))
  
  
  ##############################
  # Multiplier
  ##############################
  
  cmultiplier<-cmult(X1cand,
                     X2cand,
                     alpha=alpha,
                     tau2=tau2,l=l,TSI.pure=TSI.pure,A=A)$cm
  
  
  ##############################
  # Final Measure
  ##############################
  
  Measure<-sqrt(tau2*2/pi)*sum(sqrt(apply(A[-1,]^2,1,sum)))+
    cmultiplier*(sum(sqrt(diag(V1tilde))[-1]))
  
  
  ########
  #Output#
  ########
  D=Dcand
  X=Xcand
  V=V1tilde
  A=A
  
  return(list(D=D,X=X,V=V,A=A,Primary=Measure))
  
  
}











TSI.Construct.r.l<-function(thisn,
                            thisk,
                            model,
                            tau2,
                            n.init=10,
                            nrep=0,
                            l=0,
                            alpha=0.01,
                            CoorCand,
                            TSI.pure=FALSE,
                            StoreAll=FALSE,
                            SN=1){
  
  #Store Designs
  Store.D<-vector(mode="list")
  #Store Measures
  Store.Measure.Primary<-c()
  
  
  for(i in 1:n.init){
    print(paste("Design",i))
    
    GoodStart<-0
    
    while(GoodStart==0){
      
      Dgen<-GenInit(thisn=thisn,thisk=thisk,nrep=nrep,CoorCand=CoorCand)
      
      D.init<-as.matrix(Dgen$D)
      X<-DtoX(D.init,model)
      X1<-X[,1:(thisk+1)]
      if(qr(X1)$rank==(1+thisk)) GoodStart=1
      
    }
    
    
    #For debugging
    #Dgen<-GenInit(thisn=thisn,thisk=thisk,nrep=nrep,CoorCand=CoorCand)
    #Dhold<-Dgen
    #D.init<-as.matrix(Dgen$D)
    
    Du.init<-as.matrix(Dgen$Du)
    Dr.init<-Dgen$Dr
    if(!is.null(Dr.init)) Dr.init<-matrix(c(Dgen$Dr),nrow=nrep,ncol=thisk)
    Dr.indx<-Dgen$Dr.indx
    whichrep<-Dgen$whichrep
    
    
    #Initialize
    X<-DtoX(D.init,model)
    X1<-X[,1:(thisk+1)]
    X2<-X[,-c(1:(thisk+1))]
    
    V1<-chol2inv(chol(crossprod(X1)))
    A<-crossprod(V1,crossprod(X1,X2))
    
    BestMeasure<-Measure(D=D.init,model=model,alpha=alpha,tau2=tau2,l=l,TSI.pure=TSI.pure)$Primary
    print(paste("Initial Value: ", round(BestMeasure,4)))
    
    
    
    stop=0
    while(stop==0){
      
      MeasureVal=BestMeasure
      
      
      
      for(ii in 1:thisk){
        for(j in 1:nrow(Du.init)){
          rows<-j
          if(rows %in% whichrep){
            rows<-c(rows,(thisn-nrep)+which(rows==whichrep))  
          }
          Exch<-TSI.CoorUpdate(D=D.init,X=X,V=V1,A=A,rows=rows,col=ii,model=model,
                               alpha=alpha,tau2=tau2,l=l,CoorCand=CoorCand,TSI.pure=TSI.pure)
          
          D.init<-Exch$D
          Du.init<-D.init[1:nrow(Du.init),]
          Dr.init<-D.init[Dr.indx,]
          X<-Exch$X
          V1=Exch$V
          A=Exch$A
          MeasureVal=Exch$Primary
          
          #For debugging
          #j=j+1
          #j
          #print(paste(ii,j,MeasureVal))
        }
      }
      
      X1<-X[,1:(thisk+1)]
      X2<-X[,-c(1:(thisk+1))]
      
      
      
      ###################################
      #Block row exchange for replicates
      ###################################
      
      
      #Need library "arrangements"
      if(nrep>0){
        Xu=X[-Dr.indx,]
        Vu=chol2inv(chol(crossprod(Xu[,1:(thisk+1)])))
        
        Combos<-combinations(1:nrow(Du.init),k=nrep,replace=TRUE)
        if(dim(Combos)[1]>5000){
          Combos<-Combos[sample(1:dim(Combos)[1],5000),]
          Combos<-rbind(Combos,whichrep)
        }
        Updates<-apply(Combos,1,function(rows){
          TSI.RowUpdate(Du=Du.init,
                        Xu=X[-Dr.indx,],
                        Vu=Vu,
                        rows=rows,
                        model,alpha,tau2,l,TSI.pure)$Primary
        }
        )
        whichrep<-sort(Combos[which.min(Updates),])
        
        Dr.init<-matrix(Du.init[whichrep,],nrow=nrep,ncol=thisk)
        D.init<-rbind(Du.init,Dr.init)
        X<-DtoX(D.init,model)
        X1<-X[,1:(thisk+1)]
        X2<-X[,-c(1:(thisk+1))]
        V1<-chol2inv(chol(crossprod(X1)))
        A<-crossprod(V1,crossprod(X1,X2))
        
        MeasureVal<-Measure(D=D.init,model=model,alpha=alpha,tau2=tau2,l=l,TSI.pure=TSI.pure)$Primary
        
      }
      
      
      if( (MeasureVal/thisk) < SN & StoreAll==TRUE){
        Store.D<-append(Store.D,list(D.init))
        Store.Measure.Primary<-c(Store.Measure.Primary,MeasureVal)
      }
      
      
      #Check Convergence
      delta<-BestMeasure-MeasureVal
      if(abs(delta)<0.00001) stop=1
      
      BestMeasure=MeasureVal
      
      print(paste("Initial Value: ", round(BestMeasure,4)))
    }
    
    #if(StoreAll==FALSE){
    Store.D<-append(Store.D,list(D.init))
    Store.Measure.Primary<-c(Store.Measure.Primary,MeasureVal)
    #Store.Measure.Potential<-c(Store.Measure.Potential,MeasureVal$Potential)
    #}
    
    
    print(paste("Final Measure: ", round(BestMeasure,4)))
    
  }
  
  Dbest<-Store.D[[which.min(Store.Measure.Primary)]]
  BestMeasure.Primary<-min(Store.Measure.Primary)
  #BestMeasure.Potential<-min(Store.Measure.Potential)
  
  return(list(Dbest=Dbest,
              BestMeasure.Primary=BestMeasure.Primary,
              Store.D=Store.D,
              Store.Measure.Primary=Store.Measure.Primary))
  
}



#################################################
### Calculate LOF criterion
#################################################

Model.ERSS<-function(D,model,n.so,tau2,model.samp=NA){
  
  n<-nrow(D)
  k<-ncol(D)
  
  if(model=="ME2FI"){
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
  
  #Adjust for main effects
  X1<-cbind(1,D)
  P1<-Proj(X1)
  
  #Full model projection
  PX<-Proj(cbind(1,X))
  
  #Adjusted X2
  X2<-(PX-P1)%*%X[,-c(1:k)]
  X2[which(abs(X2)<1e-13)]=0 #Computational adjustment
  max.r<-qr(X2)$rank
  p<-ncol(X2)
  
  #X2tX2<-tcrossprod(X2)
  
  if(max.r < n.so){
    print("max.r < n.so")
    return(0)
  }else{
    
    if(length(model.samp)==1){
      model.samp<-combinations(1:p,n.so)  
      if(nrow(model.samp)>5000){
        model.samp<-model.samp[sample(1:nrow(model.samp),5000),]
      }
    }
    
    
    Lambda<-apply(model.samp,1,function(x){
      Zstar<-X2[,x]
      PZstar<-Proj(Zstar)
      
      PZstar[which(abs(PZstar)<1e-13)]=0 #Computational adjustment
      
      P.LoF<-PX-P1-PZstar
      #lambda<-tau2*sum(diag(crossprod(X2tX2,P.LoF)))
      #lambda<-sum(sort(diag(crossprod(X2tX2,P.LoF))))
      lambda<-sum(sort(round(diag(t(X2)%*%P.LoF%*%X2),2)[-x])[1:(max.r-n.so)])
      return(lambda)
    })
    
    #return(quantile(Lambda,0.01))
    return(min(Lambda))
    
  }
  
  
}

























#############################################################
#############################################################
#############################################################

# Simulation Model Generation Code

#############################################################
#############################################################
#############################################################


#Generate 2nd-order terms from given active set
# k = number of factors
# active = k length vector of 0's and 1's (1 = active)
# model = "FullQuad", "2FI", "ME"
# heredity = "strong", "weak"
# p_int, p_quad = prob of interaction and quadratic
# eff_main_offset, eff_int_offset, eff_quad_offset = offset for exponential distribution
# eff_main_rate, eff_int_rate, eff_quad_rate = rate for exponential (mean = 1/rate)

ModelSupportAndEffects<-function(active,model="FullQuad",heredity="strong",
                                 m_2FI=1,
                                 m_quad=1,
                                 eff_main_offset=1.5,eff_int_offset=1,eff_quad_offset=1,
                                 eff_main_rate=0.5,eff_int_rate=1,eff_quad_rate=1){
  
  k<-length(active)
  
  if(model=="ME"){
    
    model_terms<-active
    #Generate k effects but 0/1 model terms keeps active/inactive
    b<-model_terms*(rexp(k,rate=eff_main_rate)+eff_main_offset)
    #Random signs
    b<-b*(2*rbinom(k,1,0.5)-1)
    
  }else{
    
    ######################
    #Initialize model_terms
    ######################
    
    if(model=="2FI"){
      model_terms<-c(active,rep(0,choose(k,2)))
      b<-rep(0,length(model_terms))
      
      #Generate main effects but 0/1 model terms keeps active/inactive
      b[1:k]<-active*(rexp(k,rate=eff_main_rate)+eff_main_offset)
      #Random signs
      b[1:k]<-b[1:k]*(2*rbinom(k,1,0.5)-1)
      
    }
    
    
    
    if(model=="FullQuad"){
      model_terms<-c(active,
                     rep(0,choose(k,2)),
                     rep(0,k))
      b<-rep(0,length(model_terms))
      
      #Generate main effects but 0/1 model terms keeps active/inactive
      b[1:k]<-active*(rexp(k,rate=eff_main_rate)+eff_main_offset)
      #Random signs
      b[1:k]<-b[1:k]*(2*rbinom(k,1,0.5)-1)
    }
    
    
    #Interaction labels with 3 columns: columns 1 and 2 give int pair, 
    #column 3 = index for "terms" object, 
    
    IntLabels<-cbind(combinations(1:k,2),
                     c(1:choose(k,2)+k)
    )
    
    #Quadratic labels just 1 column, same order as main effects
    QuadLabels<-c(1:k)+k+choose(k,2)
    
    
    ######################
    #Generate active second order terms
    ######################
    
    
    if(heredity=="strong" & sum(active)>0){
      
      #Active factor labels
      ActiveLabels<-which(active==1)
      
      #2FI model, Need more than 1 active factor for possible interaction
      if(model=="2FI" & sum(active)>1){
        
        allowed_2FIs<-data.frame(combinations(ActiveLabels,2))
        
        #Terms index of overlap of allowed 2FIs and all possible 2FIs
        allowed_2FIs_terms<-c(merge(allowed_2FIs, data.frame(IntLabels))[,3])
        
        #When we had a p_int term
        #model_terms[allowed_2FIs_terms]<-rbinom(length(allowed_2FIs_terms),
        #                                  1,p_int)
        chosen_2FIs_terms<-sort(sample(allowed_2FIs_terms,m_2FI))
        model_terms[chosen_2FIs_terms]<-1
        
        
        #Generate interactions but 0/1 model terms keeps active/inactive
        b[IntLabels[,3]]<-model_terms[IntLabels[,3]]*(rexp(choose(k,2),
                                                           rate=eff_int_rate)+eff_int_offset)
        #Random signs
        b[IntLabels[,3]]<-b[IntLabels[,3]]*(2*rbinom(choose(k,2),1,0.5)-1)
        
        
      }
      
      
      if(model=="FullQuad"){
        
        if(sum(active)>1){ #Need more than 1 active factor for possible interaction
          allowed_2FIs<-data.frame(combinations(ActiveLabels,2))
          
          #Terms index of overlap of allowed 2FIs and all possible 2FIs
          allowed_2FIs_terms<-c(merge(allowed_2FIs, data.frame(IntLabels))[,3])
          
          #model_terms[allowed_2FIs_terms]<-rbinom(length(allowed_2FIs_terms),
          #                                        1,p_int)
          
          if(length(allowed_2FIs_terms)==1){
            chosen_2FIs_terms=allowed_2FIs_terms
          }else{
            chosen_2FIs_terms<-sort(sample(allowed_2FIs_terms,m_2FI))  
          }
          
          model_terms[chosen_2FIs_terms]<-1
          
          
          #Generate interactions but 0/1 model terms keeps active/inactive
          b[IntLabels[,3]]<-model_terms[IntLabels[,3]]*(rexp(choose(k,2),
                                                             rate=eff_int_rate)+eff_int_offset)
          #Random signs
          b[IntLabels[,3]]<-b[IntLabels[,3]]*(2*rbinom(choose(k,2),1,0.5)-1)
          
        }
        
        
        #model_terms[ActiveLabels+(k+choose(k,2))]<-rbinom(length(ActiveLabels),
        #                                            1,p_quad)
        
        
        chosen_quad_terms<-sort(sample(ActiveLabels+(k+choose(k,2)),m_quad))
        model_terms[chosen_quad_terms]<-1
        
        
        #Generate interactions but 0/1 model terms keeps active/inactive
        b[QuadLabels]<-model_terms[QuadLabels]*(rexp(k,rate=eff_quad_rate)+eff_quad_offset)
        #Random signs
        b[QuadLabels]<-b[QuadLabels]*(2*rbinom(k,1,0.5)-1)
        
        
        
      } #End strong heredity, full quad
      
      
    } #End strong heredity
    
    
    
    
    
    
    if(heredity=="weak"){
      
      #Active factor labels
      ActiveLabels<-which(active==1)
      
      #############
      #2FI model
      #############
      if(model=="2FI"){
        
        allowed_2FIs<-apply(IntLabels[,-3],1,function(x){sum(ActiveLabels %in% x)>0})
        
        allowed_2FIs_terms<-IntLabels[which(allowed_2FIs),3]
        
        #model_terms[allowed_2FIs_terms]<-rbinom(length(allowed_2FIs_terms),
        #                                        1,p_int)
        
        if(length(allowed_2FIs_terms)==1){
          chosen_2FIs_terms=allowed_2FIs_terms
        }else{
          chosen_2FIs_terms<-sort(sample(allowed_2FIs_terms,m_2FI))  
        }
        model_terms[chosen_2FIs_terms]<-1
        
        
        #Generate interactions but 0/1 model terms keeps active/inactive
        b[IntLabels[,3]]<-model_terms[IntLabels[,3]]*(rexp(choose(k,2),
                                                           rate=eff_int_rate)+eff_int_offset)
        #Random signs
        b[IntLabels[,3]]<-b[IntLabels[,3]]*(2*rbinom(choose(k,2),1,0.5)-1)
        
      }
      
      
      if(model=="FullQuad"){
        
        ################
        #Interactions
        ################
        allowed_2FIs<-apply(IntLabels[,-3],1,function(x){sum(ActiveLabels %in% x)>0})
        
        allowed_2FIs_terms<-IntLabels[which(allowed_2FIs),3]
        
        #model_terms[allowed_2FIs_terms]<-rbinom(length(allowed_2FIs_terms),
        #                                        1,p_int)
        
        if(length(allowed_2FIs_terms)==1){
          chosen_2FIs_terms=allowed_2FIs_terms
        }else{
          chosen_2FIs_terms<-sort(sample(allowed_2FIs_terms,m_2FI))  
        }
        model_terms[chosen_2FIs_terms]<-1
        
        
        #Generate interactions but 0/1 model terms keeps active/inactive
        b[IntLabels[,3]]<-model_terms[IntLabels[,3]]*(rexp(choose(k,2),
                                                           rate=eff_int_rate)+eff_int_offset)
        #Random signs
        b[IntLabels[,3]]<-b[IntLabels[,3]]*(2*rbinom(choose(k,2),1,0.5)-1)
        
        
        #Only allow quadratic effects with significant linear effects
        
        chosen_quad_terms<-sort(sample(ActiveLabels+(k+choose(k,2)),m_quad))
        model_terms[chosen_quad_terms]<-1
        
        #model_terms[ActiveLabels+(k+choose(k,2))]<-rbinom(length(ActiveLabels),
        #                                                  1,p_quad)
        
        #Generate interactions but 0/1 model terms keeps active/inactive
        b[QuadLabels]<-model_terms[QuadLabels]*(rexp(k,rate=eff_quad_rate)+eff_quad_offset)
        #Random signs
        b[QuadLabels]<-b[QuadLabels]*(2*rbinom(k,1,0.5)-1)
        
        
      } #End weak heredity, full quad
      
      
      
      
    } #End weak heredity
    
    
    
  }
  #End else for ME model
  
  return(list(model_terms=model_terms,b=b))
  
}
