#'Estimates the          probabilities
#'
#'@param Data is a single univariate data sequence
#'@param X_cov contains the sequence(s) of visible level covariates (in a matrix form)
#'@return Epsilon probabilities


epsilon_cov<-function(Data,  X_cov, n,r,k,OHC,kOHC,Pi,A,SAlpha,SAlog,SBeta,SBlog,LLAlpha,Phi,  Phi_cov ,Theta,P,Q,spec){


    # joint proba of OHC successive hidden states for each t ???


    Epsilon<-array(0,c(kOHC,k,n-1))
    ZEpsilon<-array(1,c(kOHC,k,n-1))

    if (OHC==1){
        A=kronecker(array(1,c(k,1)),Pi)
    }
    A=matrix(A,nrow=kOHC,ncol=kOHC)


    ### trouver les 0 dans Epsilon:
    ###      la ou   SAlpha, SBeta, Pi, output   sont nuls

    for (t in (r+1):(n-1)){
        output= array(0,c(1,k))
        for(j in 1:k){
            mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t+1-P[j]):t],   Phi_cov[j,], X_cov[,t])
            if(spec!=4){
                sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[(t+1-Q[j]):(t+1)],spec,mu)
            } else{
                sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[((t+1-P[j]-Q[j]):t)])
            }
            output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t+1]-mu)^2/(2*sigma^2)) ,.Machine$double.xmin)    ### avant sans MAX!!!!!!!!!!!!!!!!!!!!
        }

        if (t-r < OHC){
            for (i in 1:k^(t-r)){     # de k^0 ? k^(OHC-1)
                for (j in 1:k){
                    col=k^(t-r)*(j-1)+i
                    if (SAlpha[i,t]==0 || SBeta[col,t+1]==0 || Pi[i,j,t-r+1]==0 || output[j]==0){
                        ZEpsilon[i,j,t]=0
                    }
                }
            }
        } else{
            for (i in 1:kOHC){
                for (j in 1:k){
                    col=k^(OHC-1)*(j-1)+floor((i-1)/k)+1
                    if (SAlpha[i,t]==0 || SBeta[col,t+1]==0 || A[i,col]==0 || output[j]==0){
                        ZEpsilon[i,j,t]=0
                    }
                }
            }

        }
    }

    ### Computation of Epsilon

    LSAfact=0  # sum of all SAlog over time until t
    for (t in (r+1):(n-1)){
        output= array(0,c(1,k))
        for(j in 1:k){
            mu=hmtd_mu_cov(Phi[j,1:(P[j]+1)],Data[(t+1-P[j]):t],   Phi_cov[j,], X_cov[,t])
            if(spec!=4){
                sigma=hmtd_sigma_nocov(Theta[j,1:(Q[j]+1)],Data[(t+1-Q[j]):(t+1)],spec,mu) # avant: Data[((t+1-Q[j]):t+1)]
            } else{
                sigma=hmtd_sigma_arch_nocov(Theta[j,1:(Q[j]+1)],Phi[j,],P[j],Data[((t+1-P[j]-Q[j]):t)])
            }
            output[1,j]=max(1/(sqrt(2*pi)*sigma)*exp(-(Data[t+1]-mu)^2/(2*sigma^2)) ,.Machine$double.xmin)    ### avant sans MAX!!!!!!!!!!!!!!!!!!!!
        }
        LSAfact=LSAfact+SAlog[t]
        if (t-r < OHC){
            for (i in 1:k^(t-r)){     # de k^0 ? k^(OHC-1)
                for (j in 1:k){
                    if (ZEpsilon[i,j,t]!=0){
                        Epsilon [i,j,t]=LSAfact + log(SAlpha[i,t]) + log(Pi[i,j,t-r+1]) + log(output[1,j])
                    }
                }
            }
        }else{
            for (i in 1:kOHC){
                for (j  in 1:k){
                    col=k^(OHC-1)*(j-1)+floor((i-1)/k)+1
                    if (ZEpsilon[i,j,t]!=0){
                        Epsilon[i,j,t]=LSAfact + log(SAlpha[i,t]) + log(A[i,col]) + log(output[1,j])   #### partie 1 dela formule     Alpha*A*output...
                    }
                }
            }

        }

    }

    LSBfact=0

    for (t in (n-1):(r+1)){

        LSBfact=LSBfact+SBlog[t+1]  # sum of all SBlog over time from the end until t+1   =   constantes de normalisation qui se sont accumul?es
        if (t-r >= OHC){
            for (i in 1:kOHC){
                for (j in 1:k){
                    col=k^(OHC-1)*(j-1)+floor((i-1)/k)+1
                    if (ZEpsilon[i,j,t]!=0){
                        Epsilon[i,j,t]= Epsilon[i,j,t] + LSBfact + log(SBeta[col,t+1]) - LLAlpha
                        Epsilon[i,j,t]= exp(Epsilon[i,j,t])
                    }
                }
            }
        }else{
            for (i in 1:k^(t-r)){     # de k^0 ? k^(OHC-1)
                for (j in 1:k){
                    col=k^(t-r)*(j-1)+i
                    if (ZEpsilon[i,j,t]!=0){
                        Epsilon[i,j,t]= Epsilon[i,j,t] + LSBfact + log(SBeta[col,t+1]) - LLAlpha    ####  partie 2 de la formule   ... *Beta / LL
                        Epsilon[i,j,t]= exp(Epsilon[i,j,t])
                    }
                }
            }
        }
    }

    if( diff(range(Epsilon)) < 1e-99){Epsilon=Epsilon+1e-199}   # ".Machine" to see the smallest possible number for each computer. here: 2.225074e-308


    return=list(Epsilon)

}
