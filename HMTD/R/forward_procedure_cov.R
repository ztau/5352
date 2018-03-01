#'Estimates the Forward probabilities
#'
#'@param Data is a single univariate data sequence
#'@param X_cov contains the sequence(s) of visible level covariates (in a matrix form)
#'@return 1. scaled forward probabilities, 2. scaling constants and 3. incomplete log-Likelihood

### Forward procedure ###

forward_procedure_cov<-function(Data,  X_cov, n,r,k,OHC,kOHC,Pi,A,Phi,  Phi_cov  ,P,Theta,Q,spec,...){

    Ofs=1

    SAlpha<-array(0, c(kOHC,n))
    SAlog<-array(0, c(1,n))


    #<----------------------------------------------------------------------------------------------------------------------
    t=r+1             # calcul pour la colone r+1 de SAlpha


    for(j in 1:k){                                      ### calculer mu et sigma de la seq. "Data", utilisant Phi et Theta de chaque composante.
        t=max(r+1,2)
        mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t-P[j]):(t-1)],   Phi_cov[j,], X_cov[,t-1])
        if(spec!=4){
            sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[(t-Q[j]):t],spec,mu)
        } else{
            sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[(t-P[j]-Q[j]):(t-1)])
        }
        output=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t]-mu)^2/(2*sigma^2)),1e-100)  # pr chaq compos calculer la pdf de la loi normale (avec son mu et sigma) de la donn?e actuelle
        SAlpha[j,t]=Pi[Ofs,j,1]*output                # sans les covar Ofs=1, sauf si OHC=1 -> Ofs=0!
    }

    mA=array(1,dim=c(1,k))%*%SAlpha[1:k,t]/k

    if(mA==0){
        SAlpha[1:k,t]=(1/k)*array(1,dim=c(k,1))
    } else{
        SAlpha[1:k,t]=SAlpha[1:k,t]/mA

        SAlog[t]=log(mA)
    }


    if(OHC>1){
        for(t in (r+2):(r+OHC)){  #<------------------------------------------------------------------------------------------------
                                  # calcul pour les colones (r+1):(r+OHC) de SAlpha
                                  output=array(0,dim=c(1,k))
                                  for(j in 1:k){
                                      mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t-P[j]):(t-1)],   Phi_cov[j,], X_cov[,t-1])
                                      if(spec!=4){
                                          sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t-Q[j]):t)],spec,mu)
                                      } else{
                                          sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[((t-P[j]-Q[j]):(t-1))])
                                      }
                                      output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t]-mu)^2/(2*sigma^2)),1e-100)
                                  }
                                  for(i in 1:(k^(t-1-r))){
                                      for(j in 1:k){
                                          row=k^(t-1-r)*(j-1)+i
                                          SAlpha[row,t]=output[1,j]*Pi[(i-1) +Ofs,j,t-r]*SAlpha[i,t-1]   # parfois =0 car valeurs < min de la machine
                                      }
                                  }
                                  mA=array(1,dim=c(1,k^(t-r)))%*%SAlpha[1:k^(t-r),t]/(k^(t-r))  # moyenne des SAlphas de la p?riode
                                  if(mA==0){
                                      SAlpha[1:(k^(t-r)),t]=(1/(k^(t-r)))*array(1,dim=c(k^(t-r),1))
                                  } else{
                                      SAlpha[1:(k^(t-r)),t]=SAlpha[1:(k^(t-r)),t]/mA
                                      SAlog[t]=log(mA)
                                  }
        }
    }


    for(t in max((r+OHC),2):(n-1)){ #<-------------------------------------------------------------------------------------------------
                             # calcul pour les colones (r+OHC+1):(n-1) de SAlpha
                             # le calcul commence en r+OHC car il se fait pour chaque periode suivante

                             selA=1:kOHC
                             if(OHC==0){selA=1}

                             output=array(0,dim=c(1,k))            # calculer mu et sigma de la seq. "Data", utilisant Phi et Theta de chaque composante.
                             for(j in 1:k){
                                 mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t+1-P[j]):t],   Phi_cov[j,], X_cov[,t])
                                 if(spec!=4){
                                     sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[((t+1-Q[j]):t+1)],spec,mu)
                                 } else{
                                     sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[((t+1-P[j]-Q[j]):t)])
                                 }
                                 output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t+1]-mu)^2/(2*sigma^2)),1e-100) # parfois =0 car valeurs < min de la machine
                             }
                             for(j in 1:kOHC){
                                 j0=floor((j-1)/(k^(OHC-1)))+1   #number of state
                                 if(OHC==0){j0=j}
                                 SAlpha[j,t+1]=t(A[selA,j])%*%SAlpha[,t]*output[1,j0]  # parfois =0 car valeurs < min de la machine
                             }
                             mA=array(1,dim=c(1,kOHC))%*%SAlpha[,t+1]/kOHC

                             if(mA==0){
                                 SAlpha[,t+1]=(1/(kOHC))*array(1,dim=c(kOHC,1))
                             } else{
                                 SAlpha[,t+1]=SAlpha[,t+1]/mA
                                 SAlog[t+1]=log(mA)          # SAlog - constante de normalisation car SAlpha devient parfois trop petit et pose de problemes numeriques
                             }
    }

    LLAlpha=log(array(1,dim=c(1,kOHC)) %*% SAlpha[,n]) + SAlog %*% array(1,dim=c(n,1)) # LL = log( Somme_SAlpha_last * exp(tout_les_SAlog) )

    return=list(SAlpha,SAlog,LLAlpha)
}
