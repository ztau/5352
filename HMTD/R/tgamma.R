tgamma <- function(Gamma, n, r, k, OHC){

    # Computes the unconditionnal probability of each state at time t. P(X=i) no matter after which state i appears.

    # TGamma (k,n) - matrix showing the uncondit probability of each state at each time.
    
    # Gamma (k^(OHC), n) - Proba of OHC successive states at each time.
    
TGamma=array(0,c(k,n))
indic=kronecker(diag(1,k,k),array(1,c(1,k^(OHC-1))))

for (t in (r+1):n){
TGamma[,t]=indic%*%Gamma[[1]][,t]
}

return=list(TGamma)
}