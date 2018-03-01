hmtd_mu_cov<-function(Phi,Data, Phi_cov, X_cov){
    
    # Phi: (1 x p+1): phi(g,0) to phi(g,p) coeficients
    # Data: (p x 1): observations from t-p to t-1                       !!!
    # p: number of lags to compute the expectation for the mean
    # mu: expectation of component g at time t
    
    
    #P=length(Data)
    #mu=Phi[1]
    
    #for(i in 1:p){3
    #  mu=mu+Phi(i+1)*Data(p-i+1)  # -> Data x Phi
    # }
    
    mu=Phi[1] + Phi_cov[which(Phi_cov!=0)]%*%X_cov[which(Phi_cov!=0)] + if(length(Phi)>1){Phi[2:(length(Data)+1)]%*%t(t(rev(Data)))}else{0}  #length(Phi)
    
    
    return(mu)
}