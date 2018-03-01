#'Maximization procedure SEMP1
#'
#'@param fct function to be maximized (logLikelihood)
#'@param VECTPAR vector of visiible parameters Phi and Theta
#'@return logL, parameters, iterations and number of function calls

nouvelle_appr_consecutive_modification_par_importance_ordre_v2_maxiterations<- function(fct, VECTPAR, lengthPhi, lengthTheta, k, target, maxdelta, mindelta, it, envir...){

    # SEMP

    # min_param_change!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # maxdelta=0.3 - max percentage change (ex: 0,25)
    # mindelta=0.005 - min percentage change (ex: 0,01)
    # fct=LL_10
    # target=0
    # it=50

    lengthVECTPAR=length(VECTPAR)
    LLAlpha_vis_last=fct(VECTPAR)
    LL_last_max=LLAlpha_vis_last
    LL_delta=array(0, c(1,lengthVECTPAR))
    VECTPAR_elite=cbind(t(VECTPAR),0)
    ii=3
    jj=0

    changes<-array(0.2, c(1,lengthVECTPAR))
    VECThist=array(0,c(5000,lengthVECTPAR))
    LLAlphahist_up<- array(0,c(1,lengthVECTPAR))
    LLAlphahist_down<- array(0,c(1,lengthVECTPAR))
    Lvisible=array(0,c(1,it))
    Lvisible[ii]=-1
    Lvisible[ii-1]=-0.5
    nbr_fcalls=0
    importance_order=1:lengthVECTPAR

    min_param_change<-cbind(array(VECTPAR[1]/300,c(1,lengthTheta)),array(VECTPAR[lengthTheta+1]/200,c(1,k)),array(0.01,c(1,lengthVECTPAR-lengthTheta-k)))

    while(nbr_fcalls<it){
        LLAlpha_vis_last<-LL_last_max


        VECTPAR_up=VECTPAR
        VECTPAR_down=VECTPAR

        for (m in importance_order){


            if(abs(VECTPAR[m])<0.05 & m>lengthTheta){
                VECTPAR[m]=-(VECTPAR[m])                    # if parameter too small AND IS NOT a VARIANCE PARAMETER, try NEGATIVE value
            }


            VECTPAR_up[m]=VECTPAR[m]+max(VECTPAR[m]*changes[m],min_param_change[m])                                                       # positive change
            VECTPAR_down[m]=VECTPAR[m]-max(VECTPAR[m]*changes[m],min_param_change[m])


            if(any(is.na(VECTPAR_up))==TRUE){stop("VECTPAR_up has NA again in SEMP 1...") }
            if(any(is.na(VECTPAR_down))==TRUE){stop("VECTPAR_down has NA again in SEMP 1...") }

            LL_up=fct(VECTPAR_up)     # vecteur de L pour chaque param modifi? vers le haut
            jj=jj+1
            if(LL_up>LL_last_max){
                LL_delta[m]=LL_up-LL_last_max
                LL_last_max=LL_up
                VECTPAR[m]=VECTPAR_up[m]
                changes[m]=min(changes[m]*2,maxdelta)
            } else{
                LL_down=fct(VECTPAR_down)       # vecteur de L pour chaque param modifi? vers le bas
                jj=jj+1
                if(LL_down>LL_last_max){
                    LL_delta[m]=LL_down-LL_last_max
                    LL_last_max=LL_down
                    VECTPAR[m]=VECTPAR_down[m]
                    changes[m]=min(changes[m]*2,maxdelta)
                } else{
                    changes[m]<-max(changes[m]/3,mindelta)
                }
            }
        }

        importance_order<-order(-(LL_delta))




        #         negative_deltaL_place<-which((LLAlphahist_up-LLAlpha_vis_last)<0 & (LLAlphahist_down-LLAlpha_vis_last)<0)     # when the L doesn't grow up by changing 1 param
        #         changes[negative_deltaL_place]<-max(changes[negative_deltaL_place]/6,mindelta)                      # diminish its change step
        #
        #         positive_deltaL_place_up<-which((LLAlphahist_up-LLAlpha_vis_last)>0)                                             #
        #         positive_deltaL_place_down<-which((LLAlphahist_down-LLAlpha_vis_last)>0)                                         #
        #
        #         # we modify only the parameters that raise L (once upwards and downwards)                       #
        #         VECTPAR[positive_deltaL_place_up]=VECTPAR[positive_deltaL_place_up]*(1+changes[positive_deltaL_place_up])
        #         VECTPAR[positive_deltaL_place_down]=VECTPAR[positive_deltaL_place_down]*(1-changes[positive_deltaL_place_down])
        #         # we increase the  changes in this case
        #         changes[positive_deltaL_place_up]=min(changes[positive_deltaL_place_up]*3,maxdelta)                              #
        #         changes[positive_deltaL_place_down]=min(changes[positive_deltaL_place_down]*3,maxdelta)                          #

        Lvisible[ii]=LL_last_max
        print(LL_last_max)
        print(VECTPAR)

        if(Lvisible[ii]==Lvisible[ii-1]){VECTPAR_elite=rbind(VECTPAR_elite,cbind(t(VECTPAR),ii)) # save the last VECTPAR before the "jitter" AND the nbr of its corresponding function call!!!
        VECTPAR=VECTPAR+rnorm(length(VECTPAR))*0.3*VECTPAR
        VECTPAR[1:lengthTheta]=abs(VECTPAR[1:lengthTheta])
        #if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again in SEMP 1 (at the end)...") }

        LL_last_max=fct(VECTPAR)
        jj=jj+1 }

        ii=ii+1
        nbr_fcalls=jj
        if ((ii/10)%%1==0){    # si l'iteration est divisable par 10 sans d?cimales (chaque 10e iteration on note VECTPAR)
            VECThist[ii/10,]=VECTPAR
        }
        # print(changes)
        if(any(is.na(VECTPAR))==TRUE){stop("VECTPAR has NA again...") }
    }
    #return=list(max(Lvisible[Lvisible<(-0.5)]), VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ], ii-3, nbr_fcalls) # car le 2e est tjs =-0.5
    envir$LL<-max(Lvisible[Lvisible<(-0.5)])
    envir$param<-if(length(as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] )  )==0){as.array(VECTPAR)} else{as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] ) }
    envir$iter<-ii-3
    envir$fct_calls<- nbr_fcalls
}
