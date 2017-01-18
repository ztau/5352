### Calculer LL apprès avoir réestimé le niveau caché ###

LL_after_hidden_level_reestimation_cov<-function(Data2,  X_cov2, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi,  Phi_cov  ,P,maxP,Theta,Q,maxQ,spec,...){
    # Data2=Data
    # X_cov2=X_cov
    
    
    #picoef<- array(0,20)
    
   
    ########## Viterbi si A non diagonale
    
    if(exists("picov")){
        
        Vit=array(0,c(nrow(as.matrix(Data2)),ncol(as.matrix(Data2))))
        for(i in 1:nrow(Data2)){
            Vit[i,1:length(Data2[i,][!is.na(Data2[i,])])]<-Viterbi(Data2[i,], if(length(X_cov)>1){X_cov[,,i]}else{0}, length(Data2[i,][!is.na(Data2[i,])]),
                                                                   r,k,OHC,k^OHC,Pi,A,Phi,  
                                                                   Phi_cov  ,
                                                                   P,Theta ,Q)
        }
        
        if(dim(table(Vit[,n]))>1){
            #picov=cbind(Vit[,n]*0.2*rnorm(length(Vit[,n])),Vit[,n]*0.1*rnorm(length(Vit[,n])))
            multi<-multinom(Vit[,n]~picov,trace=FALSE)
            picoef<<-as.array(coef(multi))
            mostprobablestate<-t(t(apply(multi$fitted.values,1,which.is.max))) # find the most probable state for each sequence
            
            # Calculate Pi: always two dimentional in clustering (OHC=1)!
            Pi[1,,1]<-apply(array(c(1:k),c(1,k)),2,function(x) sum(mostprobablestate==x))/length(mostprobablestate)
            
            if(any(Pi<0.1)){Pi[which(Pi<0.1)]=0.1; Pi=round(aperm(apply(Pi,c(1,3),function(x) x/sum(x)), c(2,1,3)),3)}
        }
    }
    
    
    ######### Reestimate A and Pi 
    
    if(!all(A[lower.tri(A)] == 0, A[upper.tri(A)] == 0)){ #  IS A NOT DIAGONAL?
        
        
        #if(any(is.na(Phi))==TRUE){stop("Phi has NA again in beginning of LL...") }
        
        A1=A
        pi_re=Pi
        Atot=A*0
        RAtot=RA*0
        sumseqlengths=sum(seqlengths)
        for(seq in 1:nrow(Data2)){
            Data=Data2[seq,]
            X_cov=X_cov2[,,seq]
            #             if(any(is.na(Pi))==TRUE){print(seq) 
            #                 stop("Pi has NA again in LL...") }
            answers_fp<-forward_procedure_cov(Data,  X_cov, seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec=0)
            answers_bp<-backward_procedure_cov(Data,  X_cov, seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec=0)
            
            SAlpha=answers_fp[[1]]
            SAlog=answers_fp[[2]]
            LLAlpha=answers_fp[[3]]
            SBeta=answers_bp[[1]]
            SBlog=answers_bp[[2]]
            
            Epsilon<-epsilon_cov(Data,  X_cov,seqlengths[seq],r,k,OHC,kOHC,Pi,A,SAlpha,SAlog,SBeta,SBlog,LLAlpha,Phi,  Phi_cov,Theta,P,Q,spec=0)
            Gamma <- gamma1(seqlengths[seq], r, k, OHC, kOHC, SAlpha, SAlog, SBeta, SBlog, LLAlpha, Epsilon)
            
            Epsilon=array(unlist(Epsilon), dim=c(nrow(Epsilon[[1]]), ncol(Epsilon[[1]]), seqlengths[seq]-1))  # transformer de list en array
            Gamma=array(unlist(Gamma), dim=c(kOHC, seqlengths[seq]))
            
            for (g in (maxP+1):max((maxP+OHC),1)){                                                           ### réestimer Pi
                pi_re[,,g-maxP]=t(matrix(Gamma[,g],k,k^(OHC-1)))
                Pi[,,g-maxP]=t(as.matrix(pi_re[,,g-maxP]))/rowSums(t(as.matrix(pi_re[,,g-maxP,drop = FALSE])))
            }
            Pi=replace(Pi,is.nan(Pi),.Machine$double.xmin)
            
            Epssums=sapply(1:dim(Epsilon)[2], function(i) rowSums(Epsilon[,i,]))                # calculer les sommes des elements pour obtenir A
            Gamsums=sapply(1, function(i) rowSums(Gamma[,1:(seqlengths[seq]-1)]))
            
            RA_new=RA                                                                           ### réestimation de A
            for (cl in 1:k){ 
                RA_new[,cl]=t(t(Epssums[,cl])) / max(Gamsums,1e-200)
            }
            for(b in 1:kOHC){
                for (co in 1:k){
                    A1[b,(co-1)*k^(OHC-1)+ceiling(b/(k))]=RA_new[b,co]
                }
            }
            Atot=Atot+A1*(seqlengths[seq]/sumseqlengths)
            RAtot=RAtot+RA_new*(seqlengths[seq]/sumseqlengths)
            
            if(any(is.na(A1))==TRUE){print(seq)
                stop("A has NA after this sequence") }
            
        }
        for (i in 1:OHC){
            Pi[1:(k^(i-1)),,i]=(kronecker(diag(k^(i-1)),array(1,c(1,k^(OHC-i+1))))%*%RAtot)/(k^(OHC-i+1))
        }
        
        
        #### Make Pi and A sum to 1 (on all lines):
        
        Pi<-aperm(apply(Pi,c(1,3),function(x) x/sum(x)), c(2,1,3))     # works for any order of matrices Pi!
        Pi[is.nan(Pi)]=0 # because Nan results from the command before (0s)
        
        #Pi=Pi/rowSums(Pi)       # normaliser, mais ça marche pas !
        A<-Atot/rowSums(Atot)
        
    }
    
    
    
    ######### Plug in the new Pi and A to calculate the new LL   
    
    #     cc <- try(forward_procedure_cov(Data,  X_cov, n,r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec), silent=T) ##### LA DIFFERENCE EST ICI!!! Si forward_procedure fait une erreur, on rend des 0s pour son output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #     if(is(cc,"try-error")) {SAlpha<-array(0, c(kOHC,n))
    #                             SAlog<-array(0, c(1,n))} 
    #     else {
    
    
    LLlastpop=0
    for(seq in 1:nrow(Data2)){
        Data=Data2[seq,]
        X_cov=X_cov2[,,seq]
        
        answers_fp<-forward_procedure_cov(Data,  X_cov, seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec=0)
        #         SAlpha=answers_fp[[1]]
        #         SAlog=answers_fp[[2]]   #}
        #         
        #         Alpha=array(0,c(kOHC,seqlengths[seq]))
        #         for (t in 2:seqlengths[seq]){
        #             for (i in 1:kOHC){
        #                 Alpha[i,t]=SAlpha[i,t]*exp(SAlog[2:t]%*%array(1,c(t-1,1)))
        #             }
        #         }
        LLlastpop1<-answers_fp[[3]] 
        LLlastpop=LLlastpop+LLlastpop1
    }
    
    
    
    
    # ca fait les lignes jusqua 390 du code matlab hhmtd_em_fga (le reste de son code modifie les Phis et Thetas i.e. optimisation procedure)
    
    
    
    if(abs(LLlastpop)<.Machine$double.xmin)LLlastpop=.Machine$double.xmin       ### if LLlastpop < minimal number of the machine => it's = the min of the machine (not 0, because otherwise we have log(0)=-Inf)
    return=list(LLlastpop,A,Pi,if(exists("picov")){picoef}else{0})   # on prnd L et pas LOG(L)!!!!!!!!!!!!!!!!!!!!!!!!!
}