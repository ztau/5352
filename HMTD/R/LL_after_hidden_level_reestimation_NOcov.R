### Calculer LL apprès avoir réestimé le niveau caché ###

LL_after_hidden_level_reestimation_NOcov<-function(Data2, n,seqlengths,r,k,OHC,kOHC,Pi,A,RA,Phi, P,maxP,Theta,Q,picov=NA,...){#maxQ,spec
    # Data2=Data
    
    
    # Only PROBLEM: we cannot have Picovariates if A is not diagonal!!!
    # The new Pi and A are not saved!!!
    
    
    ############################################################################################################
    ############ IF A IS DIAG & PICOV EXIST -> multinomial regression
    ############################################################################################################ 
    
    if(!is.na(picov)){
        
        Vit=array(0,c(nrow(as.matrix(Data2)),ncol(as.matrix(Data2))))
        for(i in 1:nrow(Data2)){
            
            Pi[,,1]=Pi_withCovariates[i,]
            Vit[i,1:length(Data2[i,][!is.na(Data2[i,])])]<-Viterbi(Data2[i,], if(length(X_cov)>1){X_cov[,,i]}else{0}, length(Data2[i,][!is.na(Data2[i,])]),
                                                                   r,k,OHC,k^OHC,Pi,A,Phi,  
                                                                   0 ,
                                                                   P,Theta  ,Q)
        }
        
        if(dim(table(Vit[,n]))>1){    # if more than one different clusters are found
            #picov=Vit[,n]*0.2*rnorm(length(Vit[,n]))
            multi<-multinom(Vit[,n]~picov,trace=FALSE)
            picoef<<-as.array(coef(multi))
            
            z <- summary(multi)$coefficients/summary(multi)$standard.errors # 2-tailed Wald z tests to test significance of coefficients
            picoefpval <<- (1 - pnorm(abs(z), 0, 1)) * 2
            #         mostprobablestate<-t(t(apply(multi$fitted.values,1,which.is.max))) # find the most probable state for each sequence
            #         
            #         # Calculate Pi: always two dimentional in clustering (OHC=1)!
            #         Pi[1,,1]<-apply(array(c(1:k),c(1,k)),2,function(x) sum(mostprobablestate==x))/length(mostprobablestate)
            if(ncol(multi$fitted.values)==k){Pi_withCovariates<-multi$fitted.values}else{
                zeroos=array(0,c(nrow(Data2),k))
                colnames(zeroos)=as.character(1:k)
                zeroos[,c(intersect(colnames(zeroos),colnames(multi$fitted.values)))]=multi$fitted.values
                Pi_withCovariates<-zeroos
            }
            
            pithreshold<-0.08
            if(any(Pi_withCovariates<pithreshold)){Pi_withCovariates[which(Pi_withCovariates<pithreshold)]=pithreshold;Pi_withCovariates<<-round(t(apply(Pi_withCovariates,1,function(x) x/array(sum(x),c(1,length(x))))),3)}
            #if(any(Pi<0.1)){Pi[which(Pi<0.1)]=0.1; Pi=round(aperm(apply(Pi,c(1,3),function(x) x/sum(x)), c(2,1,3)),3)}
        }
        
        ############################################################################################################
        ############ ELSE IF A IS DIAG & NO PICOV -> USE VITERBI TO "AJUST" Pi (make it proportionnal to the current size of the groups):
        ############################################################################################################        
    } else if(all(A[lower.tri(A)] == 0, A[upper.tri(A)] == 0)){
        Vit=array(0,c(nrow(as.matrix(Data2)),ncol(as.matrix(Data2))))
        for(i in 1:nrow(Data2)){
            
            Vit[i,1:length(Data2[i,][!is.na(Data2[i,])])]<-Viterbi(Data2[i,], if(length(X_cov)>1){X_cov[,,i]}else{0}, length(Data2[i,][!is.na(Data2[i,])]),
                                                                   r,k,OHC,k^OHC,Pi,A,Phi, 0 ,P,Theta  ,Q)
        }
        
        tab<-table(factor(Vit[,4], levels = 1:k))#table(Vit[,4])
        Pi[1,,1]<-as.numeric(apply(tab,1,function(x){x=x/nrow(Data2)}))
        if(any(Pi<0.1)){Pi[which(Pi<0.1)]=0.1; Pi=round(aperm(apply(Pi,c(1,3),function(x) x/sum(x)), c(2,1,3)),3)}
        
        
        
        ############################################################################################################
        ######### IS A NOT DIAGONAL? -> Reestimate A and Pi 
        ############################################################################################################
        
    } else{    
        
        
        
        A1=A
        pi_re=Pi
        Atot=A*0
        RAtot=RA*0
        sumseqlengths=sum(seqlengths)
        for(seq in 1:nrow(Data2)){
            Data=Data2[seq,]
            #if(any(is.na(Pi))==TRUE){print(seq) 
            #                         stop("Pi has NA again in LL...") }
            
            answers_fp<-forward_procedure_NOcov(Data, seqlengths[seq],r,k,OHC,kOHC,Pi[.....],A,Phi, P,Theta,Q,spec=0)
            answers_bp<-backward_procedure_NOcov(Data,  seqlengths[seq],r,k,OHC,kOHC,Pi[.....],A,Phi, P,Theta,Q,spec=0)
            
            SAlpha=answers_fp[[1]]
            SAlog=answers_fp[[2]]
            LLAlpha=answers_fp[[3]]
            SBeta=answers_bp[[1]]
            SBlog=answers_bp[[2]]
            
            Epsilon<-epsilon_NOcov(Data, seqlengths[seq],r,k,OHC,kOHC,Pi,A,SAlpha,SAlog,SBeta,SBlog,LLAlpha,Phi, Theta,P,Q,spec=0)
            Gamma <- gamma1(seqlengths[seq], r, k, OHC, kOHC, SAlpha, SAlog, SBeta, SBlog, LLAlpha, Epsilon)
            
            Epsilon=array(unlist(Epsilon), dim=c(nrow(Epsilon[[1]]), ncol(Epsilon[[1]]), seqlengths[seq]-1))  # transformer de list en array
            Gamma=array(unlist(Gamma), dim=c(kOHC, seqlengths[seq]))
            
            for (g in (maxP+1):(maxP+OHC)){                                                           ### réestimer Pi
                pi_re[,,g-maxP]=t(matrix(Gamma[,g],k,k^(OHC-1)))
                Pi[,,g-maxP]=t(as.matrix(pi_re[,,g-maxP]))/rowSums(t(as.matrix(pi_re[,,g-maxP,drop = FALSE])))
            }
            Pi=replace(Pi,is.nan(Pi),0)
            
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
            
            #print(Epssums) print(Gamma) print(Gamsums) print(RA_new)
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
    
    ############################################################################################################    
    ######################### Plug in the new Pi and A to calculate the new LL  
    ############################################################################################################
    LLlastpop=0
    for(seq in 1:nrow(Data2)){
        Data=Data2[seq,]
        if(!is.na(picov)){Pi[,,1]=Pi_withCovariates[i,]}
        answers_fp<-forward_procedure_NOcov(Data, seqlengths[seq],r,k,OHC,kOHC,Pi,A,Phi, P,Theta,Q,spec=0)
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
    return=list(LLlastpop,A,Pi,if(!is.na(picov)){picoef}else{0},if(!is.na(picov)){Pi_withCovariates}else{0})   # on prnd L et pas LOG(L)!!!!!!!!!!!!!!!!!!!!!!!!!
}