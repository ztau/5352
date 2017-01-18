hmtd_sigma_nocov<-function(Theta,Data,spec,mu){
  
  # Theta: (1 x q+2) theta(g,0) to theta(g,q)
  # Data: (q+1 x 1) observations from t-q to t
    # Specification of sigma:
          #1. use of X(t)
          #2. use of the expectation
          #3. use of the mean of lag values
          #4. ARCH  (see hmtd sigma arch)
  
  q=length(Data)-1
  
  if(spec==1){
    Data[q+1]=0          #?????
  } else if(spec==2){
    Data[q+1]=mu
  } else if(spec==3){
    Data[q+1]=mean(Data[1:q])
  }
 # sigma=Theta[1]
  
 # for(j in 1:q){
 #  sigma=sigma+Theta[i+1]*(Data[q+1]-Data[q-j+1])^2 # Sum (thetas x (sq diff btw "1,2OR3" and each of othe past obs))
 # }
  
  if(length(Theta)>1){
  sigma=Theta[1]+Theta[2:length(Theta)]%*%(Data[q+1]-t(t(rev(Data[1:q]))))^2
  } else{
    sigma=Theta[1]    ########### LE BOUCLE IF EST AJOUTE PAR MOI. EST-IL JUSTE??????????????????????????????????????????????????????????????????
  }
    
  sigma=sqrt(abs(sigma))
  
  return(sigma)
}