# VITERBI algorithm
# finds the most probable sequence of hidden states

Viterbi<-function(Data, X_cov, n,r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q){
    
#     Data=Data2[i,]
#     X_cov=if(length(X_cov)>1){X_cov[,,i]}else{0}
#     n=length(Data2[i,][!is.na(Data2[i,])])
#     kOHC=k^OHC  
#     Phi_cov=0 
    
    if(length(X_cov)>1){mean=hmtd_mu_cov}else{mean=hmtd_mu_NOcov}
    # Data - one oserved sequence
    # output = PXknowingZ
    # Mt = 
    # avgMt = average of Mt
    #
    # OHC=0 treated like OHC=1, with A=Pi in this case
    
    # Initialisation
    Mt<-array(0,c(kOHC,n))
    avgMt<-array(0,c(1,n))
    SeqHS<-array(0,c(1,n))
    
    t=r+1
    
    output=array(0,c(1,k))
    for (j in 1:k){
        mu=mean(Phi[j,1:(P[j]+1)],Data[(t-P[j]):(t-1)],   Phi_cov[j,], X_cov[,t])
        sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t-Q[j]):t)],0,mu)
        output[1,j]=1/(sqrt(2*pi)*sigma)*exp(-(Data[t]-mu)^2/(2*sigma^2))
    }
    for (i in 1:k){
        Mt[i,r+1]=Pi[1,i,1]*output[i]   # =P(Zr)*P(Xr|Zr)
    }
    avgMt[r+1]=array(1,c(1,k))%*%(Mt[1:k,r+1]/k)
    
    if(avgMt[r+1]==0){
        Mt[1:k,r+1]=(1/k)*array(1,c(k,1))
    } else{
        Mt[,r+1]=Mt[,r+1]/avgMt[r+1]
    }
    
    
    # t=r+2...r+OHC
    
    if(OHC>1){
        for (t in (r+2):(r+OHC)){
            selPi=1:k^(t-r-1)
            output=array(0,c(1,k))
            
            for (j in 1:k){
                mu=mean(Phi[j,1:(P[j]+1)],Data[(t-P[j]):(t-1)],   Phi_cov[j,], X_cov[,t])
                sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t-Q[j]):t)],0,mu)
                output[1,j]=1/(sqrt(2*pi)*sigma)*exp(-(Data[t]-mu)^2/(2*sigma^2))
            }
            for (j in 1:k^(t-r)){
                j0=floor((j-1)/(k^(t-r-1)))+1
                calc=Mt[1:k^(t-r-1),t-1]*Pi[selPi,j0,t-r]
                calc2=max(calc)
                Mt[j,t]=output[j0]*calc2
            }
            avgMt[t]=(array(1,c(1,k^(t-r)))%*%Mt[1:k^(t-r),t])/k^(t-r)
            if(avgMt[t]==0){
                Mt[1:k^(t-r),t]=(1/k^(t-r))*array(1,c(k^(t-r),1))
            } else{
                Mt[1:k^(t-r),t]=Mt[1:k^(t-r),t]/avgMt[t]
            }
        }
    }
    
    
    
    # t = r+OHC+1 ... n
    
    for(t in (r+OHC+1):n){
        
        selA=1:k^OHC
        output=array(0,c(1,k))
        
        for (j in 1:k){
            mu=mean(Phi[j,1:(P[j]+1)],Data[(t-P[j]):(t-1)],   Phi_cov[j,], X_cov[,t])
            sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t-Q[j]):t)],0,mu)
            output[1,j]=1/(sqrt(2*pi)*sigma)*exp(-(Data[t]-mu)^2/(2*sigma^2))
        }
        for (j in 1:k^OHC){
            j0=floor((j-1)/(k^(OHC-1)))+1
            calc=Mt[,t-1]*A[selA,j]
            calc2=max(calc)
            Mt[j,t]=output[j0]*calc2
        }
        avgMt[t]=(array(1,c(1,k^OHC))%*%Mt[,t])/k^OHC
        
        if(avgMt[t]==0){
            Mt[,t]=(1/k^OHC)*array(1,c(k^OHC,1))
        } else{
            Mt[,t]=Mt[,t]/avgMt[t]
            Mt[,t][is.infinite(Mt[,t])] <-.Machine$double.xmax   ################################### donne parfois Inf sinon
        }
        
    }
    
    
    
    
    
    
    # Sequence of hidden states
    
    # t=n
    calc4=max(Mt[,n])
    calc5=which(calc4==Mt[,n])
    calc6=calc5[1]
    calc7=floor((calc6-1)/k^(OHC-1))+1
    SeqHS[n]=calc7
    
    for(t in (n-1):(r+OHC)){
        selA=1:k^OHC
        calc4=max(Mt[,t]*A[selA,calc6])
        calc5=which(Mt[,t]*A[selA,calc6]==calc4)
        calc6=calc5[1]
        calc7=floor((calc6-1)/k^(OHC-1))+1
        SeqHS[t]=calc7
    }
    
    
    if(r+OHC>r+1){                   # !!!
        for(t in (r+OHC-1):(r+1)){   #before: for(t in (r+OHC-1):min((r+1),1)){
            selPi=1:k^(t-r)
            calc4=max(Mt[1:k^(t-r),t]*Pi[selPi,calc7,t-r+1])
            calc5=which(Mt[1:k^(t-r),t]*Pi[selPi,calc7,t-r+1]==calc4)
            calc6=calc5[1]
            calc7=floor((calc6-1)/k^(t-r-1))+1 #before: calc7=floor((calc6-1)/k^(OHC-1))+1
            SeqHS[t]=calc7
        }
    }
    
    
    return(SeqHS)
}
