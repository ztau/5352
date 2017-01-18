hmtd_sigma_arch_nocov<-function(Theta,Phi,p,Data){
  
  # Si pas de covariables -> les mêmes notations qu'avant
  # Theta (1 x q+1): theta(g,0) to theta(g,q)
  # Phi (1 x pp): phi(g,0) to phi(g,p)
  # Data(p+q x 1) vector of zise p+q: observations from t-p-q to t-1
  
  # pp length of Phi
  
  q=length(Theta)-1
  pp=length(Phi)
  sigma=Theta[1]

for(j in 1:q){
  perror=Data[p+q-j+1]-Phi[1] # 
  
  for(i in 1:p){
    perror=perror-Phi[i+1]*Data[p+q-j+1-i]
  }
  sigma=sigma+Theta[j+1]*(perror^2)
}
  sigma=sqrt(sigma)
  
  return(sigma)
}

# to modelise  the variance using ARCH