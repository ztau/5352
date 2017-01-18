OPTIMISATION1_universal_NEW_Adiag_covminmax<-function(Data,  X_cov, Phi_cov, phi_cov_places, k,OHC, P,Q,procedure=c("SEMP"), it, envir, Adiag, covmin=NA, covmax=NA,...){
    
    # procedures:  "S", "SE", "SEM", "SEMP", "NM", "LBFGSB", "DE", "PSO"
    
    
    library("nnet")
    
    source("hmtd_mu_cov.R")       
    source("hmtd_sigma_nocov.R")
    source("forward_procedure_cov.R")
    source("backward_procedure_cov.R")
    source("hmtd_sigma_arch_nocov.R")
    source("epsilon_cov.R")
    source("gamma1.R")
    source("tgamma.R")
    
    source("hmtd_mu_NOcov.R")       
    source("forward_procedure_NOcov.R")
    source("backward_procedure_NOcov.R")
    source("epsilon_NOcov.R")
    
    source("LL_after_hidden_level_reestimation_cov.R")
    source("LL_after_hidden_level_reestimation_NOcov.R")
    
    
    source("SEM_v2.R")
    source("SEMP_v2.R")
    
    source("Viterbi.R")
    
    r=max(P,Q)
    seqlengths<-array(0,c(nrow(Data),1))
    for(i in 1:nrow(Data)){
        seqlengths[i]=min(min(which(is.na(Data[i,]),arr.ind=TRUE)-1),ncol(Data))
    }
    Data=array(Data[c(which(seqlengths>=r+OHC+2)),],c(length(c(which(seqlengths>=r+OHC+2))),ncol(Data))) #,drop=FALSE               # min length = r+OHC+2 (voir fp())
    seqlengths=seqlengths[which(seqlengths>=r+OHC+2)]
    Data=matrix(as.numeric(Data),c(nrow(Data),ncol(Data)))
    
    n=ncol(Data)
    spec=0
    
    kOHC<-k^OHC
    if(OHC==0){kOHC=k}
    
    
    
    maxP=max(P)
    maxQ=max(Q)
    
    
    Phi=matrix(cbind(matrix(rnorm(k,mean(Data,na.rm=TRUE),sqrt(mean(apply(Data,1,function(x) var(x,na.rm=TRUE))))/4),1,k),matrix(runif(k*maxP,-0.2,0.2),1,k*maxP) ),c(k,maxP+1))
    Theta=matrix(abs(rnorm(k*(maxQ+1),sqrt(mean(apply(Data,1,function(x) var(x,na.rm=TRUE)))),sqrt(mean(apply(Data,1,function(x) var(x,na.rm=TRUE))))/4)), c(k,maxQ+1))
    
    
    r<-max(P,Q)
    
    if(OHC==0){
        Pi<-array(0,dim=c(1,k,1))     # initialise Pi - distribution of the first OHC hidden states
        Pi[1,,1]=array(1/k, 1,k)
    } else{
        Pi<-array(0,dim=c(k^(OHC-1),k,OHC))     # initialise Pi - distribution of the first OHC hidden states
        for (h in 1:OHC){
            Pi[1:(k^(h-1)),,h]=array(1/k, c((k^(h-1)),k))
        }
    }
    
    #if(OHC==1){Pi=Pi[,,1]}
    # Pi<-array(0,dim=c(k^(OHC-1),k,OHC))
    if(OHC==0){
        A<-array(0,dim=c(1,k))
        RA<-array(0,dim=c(1,k))
        RA<-0.1+matrix(runif(1*k), ncol=k)  
        RA[1,]=RA[1,]/sum(RA)
        A=RA
    } else{
        A<-array(0,dim=c(kOHC,kOHC))
        RA<-array(0,dim=c(kOHC,k))
        if(Adiag==1){
            RA=diag(1,k,k)
            A=RA
        } else{
            RA<-0.1+matrix(runif(kOHC*k), ncol=k)   # 0.1 pour éliminer les probabilités trop petites ?
            for(i in 1:kOHC){
                RA[i,]=RA[i,]/sum(RA[i,][1:length(RA[i,])])
                for(b in 1:kOHC){
                    for (co in 1:k){
                        A[b,(co-1)*k^(OHC-1)+ceiling(b/(k))]=RA[b,co]
                    }
                }
            }
        }
    }
    
    
    
    
    if(length(X_cov)<2){   # are there any visible level covariates?
        
        
        LL_seq<-array(0,c(nrow(Data),1))
        for(seq in 1:nrow(Data)){
            
            
            
            answers_fp<-forward_procedure_NOcov(Data[seq,], seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi, P,Theta,Q,spec)
            
            
            LL_seq[seq,]<-answers_fp[[3]]
            
        }
        
        LLAlpha=sum(LL_seq)
        
        # 3 # MODIFICATION DES PARAMETRES ----------------------------------------------------------------------------------------------
        
        lengthTheta=length(Theta)
        lengthPhi=length(Phi) 
        
        # create a vector of all visible parameters:
        VECTPAR=c(Theta, Phi)   
        lengthVECTPAR=length(VECTPAR)
        
        
        
        # Limits of the parameters:
        #####    
        
        Theta_Min= 0.5*min(apply(Data,1,function(x) var(x[!is.na(x)]))) #0.5*var of the data
        if(Theta_Min==0)Theta_Min=0.0001
        Theta_Max= 1.5*max(apply(Data,1,function(x) var(x[!is.na(x)]))) #1.5*var of the data
        Phi1_Min= min(Data[!is.na(Data)])#*(min(Data)/(2*max(Data)))   # min of the data   
        Phi1_Max= 1.5*max(Data[!is.na(Data)])#*((2*max(Data))/min(Data))  # max(Data)*1.5
        
        Phi_Min=-2
        Phi_Max=2
        VECT_limits_min=c(rep(Theta_Min,k),rep(Phi1_Min,k),rep(Phi_Min, k*(ncol(Phi)-1)) )                        
        VECT_limits_max=c(rep(Theta_Max,k),rep(Phi1_Max,k),rep(Phi_Max, k*(ncol(Phi)-1)) )                         
        
        # Function to optimise (logL)
        
        
        LL_10<-function(VECTPAR){
            if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again in LL_10...") }
            Phi=matrix(VECTPAR[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
            if(any(is.na(Phi))==TRUE){stop("Phi has NA again in LL_10...") }
            Phi_cov=array(0,c(k,nrow(X_cov)))
            Phi_cov[phi_cov_places] = VECTPAR[(lengthTheta+lengthPhi+1):length(VECTPAR)]    #,k,ncol(Phi_cov))
            Theta=matrix(VECTPAR[1:(lengthTheta)],k,max(Q)+1)
            Theta[Theta<c(Theta_Min)]=Theta_Min
            VECTPAR[which(Theta<c(Theta_Min))]=Theta_Min
            output=LL_after_hidden_level_reestimation_NOcov(Data, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi  ,P,maxP,Theta,Q,maxQ=maxQ,spec=0,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
            A<<-output[[2]]
            Pi<<-output[[3]]
            return=output[[1]] # output[[2]] output[[3]]
            
        }
        
        LL_11<-function(VECTPAR){
            if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again in LL_11...") }
            Phi=matrix(VECTPAR[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
            if(any(is.na(Phi))==TRUE){stop("Phi has NA again in LL_11...") }
            Phi_cov=array(0,c(k,nrow(X_cov)))
            Phi_cov[phi_cov_places] = VECTPAR[(lengthTheta+lengthPhi+1):length(VECTPAR)]    #,k,ncol(Phi_cov))
            Theta=matrix(VECTPAR[1:(lengthTheta)],k,max(Q)+1)
            Theta[Theta<c(Theta_Min)]=Theta_Min
            VECTPAR[which(Theta<=c(Theta_Min))]=Theta_Min
            output=LL_after_hidden_level_reestimation_NOcov(Data, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi  ,P,maxP,Theta,Q,maxQ=maxQ,spec=0,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
            A<<-output[[2]]
            Pi<<-output[[3]]  
            return=-output[[1]]
            
        }
        
        
        if (procedure=="SEM_v2"){
            SEM_v2(LL_10,VECTPAR,lengthPhi, lengthTheta,k,0, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 211          
        }
        else if (procedure=="SEMP_v2"){
            SEMP_v2(LL_10,VECTPAR,lengthPhi, lengthTheta,k,0, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 211          
        }
        else if (procedure=="SEMP"){
            nouvelle_appr_consecutive_modification_par_importance_ordre_v2_maxiterations(LL_10,VECTPAR, lengthPhi, lengthTheta, k, 0, 0.3, 0.005, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 221          
        }
        else if (procedure=="NM"){
            rslt <- optim(VECTPAR,LL_11, control=list(maxit=it))# Nelder-Mead
            envir$param=rslt$par
            envir$LL=-rslt$value
            envir$fct_calls=rslt$counts[1]
        }
        else if (procedure=="L-BFGS-B"){
            
            rslt <- optim(VECTPAR,LL_11,method ="L-BFGS-B", control=list(maxit=it),lower= VECT_limits_min ,upper= VECT_limits_max)#
            envir$param=rslt$par
            envir$LL=-rslt$value
            envir$fct_calls=rslt$counts[1]
        }
        else if (procedure=="DE"){
            library(DEoptim)
            rslt <- DEoptim(LL_11,lower= VECT_limits_min ,upper= VECT_limits_max, control=list(itermax=it/(10*length(VECT_limits_min)))) # (10*length(lower)) is the default number of population members
            envir$param=rslt$optim$bestmem
            envir$LL=-rslt$optim$bestval
            envir$fct_calls=-rslt$optim$nfeval
        }
        else if (procedure=="GA"){
            library(GA)
            rslt <- ga(type = "real-valued",fitness = function(VECTPAR) LL_10(VECTPAR),min = VECT_limits_min , max = VECT_limits_max, maxiter = it/50) # 50 function calls in one iteration
            envir$param=rslt@solution
            envir$LL=rslt@fitnessValue
            envir$fct_calls=rslt@iter*rslt@popSize
        }
        else {
            library("pso")
            rslt <- psoptim(rep(NA,length(VECTPAR)),function(x) LL_10(x),lower= VECT_limits_min ,upper= VECT_limits_max,control=list(fnscale=-1,maxf=it))#,maxit=round(it/15)) ) # 15 function calls in one iteration
            envir$param=rslt[[1]]
            envir$LL=rslt[[2]]
            envir$fct_calls=rslt[[3]]
        }
        
        
        # Compute A and Pi with the optimal VECTPAR
        Phi=matrix(envir$param[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
        Theta=matrix(envir$param[1:(lengthTheta)],k,max(Q)+1)
        Theta[Theta<c(Theta_Min)]=Theta_Min
        envir$param[which(Theta<c(Theta_Min))]=Theta_Min
        out=LL_after_hidden_level_reestimation_NOcov(Data, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi, P,maxP,Theta,Q,maxQ,spec,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
        envir$LL<-out[[1]]
        envir$A<-out[[2]]
        envir$Pi<-out[[3]]
        envir$picoef<-out[[4]]
        
        
        
        
        
        
    } else{
        
        LL_seq<-array(0,c(nrow(Data),1))
        for(seq in 1:nrow(Data)){
            
            
            
            answers_fp<-forward_procedure_cov(Data[seq,],  X_cov[,,seq], seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec)
            
            LL_seq[seq,]<-answers_fp[[3]]
            
        }
        
        LLAlpha=sum(LL_seq)
        
        # 3 # MODIFICATION DES PARAMETRES ----------------------------------------------------------------------------------------------
        
        lengthTheta=length(Theta)
        lengthPhi=length(Phi) 
        lengthPhicov=length(Phi_cov)
        
        # create a vector of all visible parameters:
        VECTPAR=c(Theta, Phi)   
        lengthVECTPAR=length(VECTPAR)
        
        
        
        # Limits of the parameters:
        #####      
        Theta_Min= 0.5*min(apply(Data,1,function(x) var(x[!is.na(x)])))
        if(Theta_Min==0)Theta_Min=0.0001
        Theta_Max= 1.5*max(apply(Data,1,function(x) var(x[!is.na(x)])))
        Phi1_Min= min(Data[!is.na(Data)])-0.5*abs(max(Data[!is.na(Data)])-min(Data[!is.na(Data)]))#*(min(Data)/(2*max(Data)))   # min(Data)/2   
        Phi1_Max= 1.5*max(Data[!is.na(Data)])#*((2*max(Data))/min(Data))  # max(Data)*2
        Phi_Min=-2
        Phi_Max=2
        covmin
        covmax
        VECT_limits_min=c(rep(Theta_Min,k),rep(Phi1_Min,k),rep(Phi_Min, k*(ncol(Phi)-1)),if(is.na(covmin)){rep(-6,length(Phi_cov[Phi_cov!=0]) )}else{covmin} )                        
        VECT_limits_max=c(rep(Theta_Max,k),rep(Phi1_Max,k),rep(Phi_Max, k*(ncol(Phi)-1)),if(is.na(covmax)){rep(6,length(Phi_cov[Phi_cov!=0]) )}else{covmax} )                         
        #----    
        
        # Function to optimise (logL)
        
        LL_10<-function(VECTPAR){
            if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again in LL_10...") }
            Phi=matrix(VECTPAR[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
            if(any(is.na(Phi))==TRUE){stop("Phi has NA again in LL_10...") }
            Phi_cov=array(0,c(k,nrow(X_cov)))
            Phi_cov[phi_cov_places] = VECTPAR[(lengthTheta+lengthPhi+1):length(VECTPAR)]    #,k,ncol(Phi_cov))
            Theta=matrix(VECTPAR[1:(lengthTheta)],k,max(Q)+1)
            Theta[Theta<c(Theta_Min)]=Theta_Min
            VECTPAR[which(Theta<c(Theta_Min))]=Theta_Min
            output=LL_after_hidden_level_reestimation_cov(Data,  X_cov, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi,  Phi_cov  ,P,maxP,Theta,Q,maxQ=maxQ,spec=spec,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
            A<<-output[[2]]
            Pi<<-output[[3]]
            
            return=output[[1]] 
            
        }
        
        LL_11<-function(VECTPAR){
            if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again in LL_11...") }
            Phi=matrix(VECTPAR[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
            if(any(is.na(Phi))==TRUE){stop("Phi has NA again in LL_11...") }
            Phi_cov=array(0,c(k,nrow(X_cov)))
            Phi_cov[phi_cov_places] = VECTPAR[(lengthTheta+lengthPhi+1):length(VECTPAR)]    #,k,ncol(Phi_cov))
            Theta=matrix(VECTPAR[1:(lengthTheta)],k,max(Q)+1)
            Theta[Theta<c(Theta_Min)]=Theta_Min
            VECTPAR[which(Theta<=c(Theta_Min))]=Theta_Min
            output=LL_after_hidden_level_reestimation_cov(Data,  X_cov, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi,  Phi_cov  ,P,maxP,Theta,Q,maxQ=maxQ,spec=spec,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
            A<<-output[[2]]
            Pi<<-output[[3]]
            
            return=-output[[1]]
            
        }
        
        
        
        if (procedure=="SEM_v2"){
            SEM_v2(LL_10,VECTPAR,lengthPhi, lengthTheta,k,0, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 211          
        }
        else if (procedure=="SEMP_v2"){
            SEMP_v2(LL_10,VECTPAR,lengthPhi, lengthTheta,k,0, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 211          
        }
        else if (procedure=="SEMP"){
            nouvelle_appr_consecutive_modification_par_importance_ordre_v2_maxiterations(LL_10,VECTPAR, lengthPhi, lengthTheta, k, 0, 0.3, 0.005, it, envir,lower= VECT_limits_min ,upper= VECT_limits_max )### 0.1!!!! vient de new 221          
        }
        else if (procedure=="NM"){
            rslt <- optim(VECTPAR,LL_11, control=list(maxit=it))# Nelder-Mead
            envir$param=rslt$par
            envir$LL=-rslt$value
            envir$fct_calls=rslt$counts[1]
        }
        else if (procedure=="L-BFGS-B"){
            
            rslt <- optim(VECTPAR,LL_11,method ="L-BFGS-B", control=list(maxit=it),lower= VECT_limits_min ,upper= VECT_limits_max)#
            envir$param=rslt$par
            envir$LL=-rslt$value
            envir$fct_calls=rslt$counts[1]
        }
        else if (procedure=="DE"){
            library(DEoptim)
            rslt <- DEoptim(LL_11,lower= VECT_limits_min ,upper= VECT_limits_max, control=list(itermax=it/(10*length(lower)))) # (10*length(lower)) is the default number of population members
            envir$param=rslt$optim$bestmem
            envir$LL=-rslt$optim$bestval
            envir$fct_calls=-rslt$optim$nfeval
        }
        else if (procedure=="GA"){
            library(GA)
            rslt <- ga(type = "real-valued",fitness = function(VECTPAR) LL_10(VECTPAR),min = VECT_limits_min , max = VECT_limits_max, maxiter = it/50) # 50 function calls in one iteration
            envir$param=rslt@solution
            envir$LL=rslt@fitnessValue
            envir$fct_calls=rslt@iter*rslt@popSize
        }
        else {
            library("pso")
            rslt <- psoptim(rep(NA,length(VECTPAR)),function(x) LL_10(x),lower= VECT_limits_min ,upper= VECT_limits_max,control=list(fnscale=-1,maxf=it)) # 15 function calls in one iteration
            envir$param=rslt[[1]]
            envir$LL=rslt[[2]]
            envir$fct_calls=rslt[[3]]
        }
        
        
        # Compute A and Pi with the optimal VECTPAR
        Phi=matrix(envir$param[(lengthTheta+1):(lengthTheta+lengthPhi)],k,max(P)+1)
        Phi_cov=matrix(envir$param[(lengthTheta+lengthPhi+1):length(envir$param)],k,ncol(Phi_cov))
        Theta=matrix(envir$param[1:(lengthTheta)],k,max(Q)+1)
        Theta[Theta<c(Theta_Min)]=Theta_Min
        envir$param[which(Theta<c(Theta_Min))]=Theta_Min
        out=LL_after_hidden_level_reestimation_cov(Data,  X_cov, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi,  Phi_cov  ,P,maxP,Theta,Q,maxQ,spec,picov) # "_ESSAI": Si forward_procedure_cov fait une erreur, on rend des 0s pour son output!!
        envir$LL<-out[[1]]
        envir$A<-out[[2]]
        envir$Pi<-out[[3]]
        envir$picoef<-out[[4]]
        
    }
    print("Done")
    
}