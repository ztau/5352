backward_procedure_cov<-function(Data,  X_cov, n,r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec,...){
    
    Ofs=1
    
    SBeta=array(0,c(kOHC,n))
    SBlog=array(0,c(1,n))
    
#     if (OHC==0){
#         A=kronecker(array(1,dim=c(k,1)),Pi)          #   if OHC=0, répète k fois Pi (juste 2 dim dans ce cas)
#     }
    
    t=n                                             #   on commence par le dernier période
    SBeta[,n]=array(1,c(kOHC,1))
    
    for (t in (n-1):max((OHC+r),1)){                       #   pour chaque période de la fin jusqua OHC+r
        
        output=array(0,c(1,k))
        
        for (j in 1:k){                              #   Pr chaque composante faire les mêmes calculs de mu et sigma comme forward...
            mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t+1-P[j]):t],   Phi_cov[j,], X_cov[,t])
            if(spec!=4){
                sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t+1-Q[j]):t+1)],spec,mu)
            } else{
                sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[((t+1-P[j]-Q[j]):t)])
            }
            output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t+1]-mu)^2/(2*sigma^2)),1e-100)   # densité de la loi normale avec les param de la compos. k.
        }
        for (i in 1:kOHC){
            for (j in 1:k){
                col=k^(OHC-1)*(j-1) + floor((i-1)/k) + 1                   # if k=3 OHC=3 <-1,10,13, 1,10,13, 1,10,13, 2,11,20, 2,11,20, 2,11,20, 3,12,21...
                SBeta[i,t]=SBeta[i,t] + if(OHC==0){A[1,col]}else{A[i-1+Ofs,col]}%*%SBeta[col,t+1]*output[j]
            }
        }
        mB=array(1,c(1,kOHC))%*%SBeta[,t]/kOHC
        
        if (mB==0){
            SBeta[,t]=(1/kOHC)*array(1,c(kOHC,1))
        } else{
            SBeta[,t]=SBeta[,t]/mB
            SBlog[t]=log(mB)
        }
    }
    
    
    
    if (OHC>1){
        for (t in (r+OHC-1):(r+1)){
            output=array(0,c(1,k))
            for (j in 1:k){
                mu=hmtd_mu_cov(Phi[j,],Data[(t+1-P[j]):t],   Phi_cov[j,], X_cov[,t]) # avant: hmtd_mu_nocov(Phi[j,],Data[(t+1-P[j]):t])!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if(spec!=4){
                    sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[(t+1-Q[j]):(t+1)],spec,mu)
                } else{
                    sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[(t+1-P[j]-Q[j]):t])
                }
                output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t+1]-mu)^2/(2*sigma^2)),1e-100)
            }
            for (i in 1:k^(t-r)){
                for (j in 1:k){
                    col=k^(t-r)*(j-1)+i
                    SBeta[i,t]=SBeta[i,t]+Pi[i-1+Ofs,j,t-r+1]%*%SBeta[col,t+1]*output[1,j]
                }
            }
            mB=array(1,kOHC)%*%SBeta[,t]/(k^(t-r))
            
            if (mB==0){
                SBeta[1:k^(t-r),t]=(1/(k^(t-r)))*array(1,c(k^(t-r),1))
            } else {
                SBeta[,t]=SBeta[,t]/mB
                SBlog[t]=log(mB)
            }
        }
    }
    
    return=list(SBeta,SBlog)
}
