#'Maximization procedure SEM
#'
#'@param fct function to be maximized (logLikelihood)
#'@param VECTPAR vector of visiible parameters Phi and Theta
#'@return logL, parameters, iterations and number of function calls

nouvelle_appr_consecutive_modification_v2_maxiterations<- function(fct, VECTPAR,lengthPhi,lengthTheta,k, target, maxdelta, mindelta, it, envir...){

    # SEM

    # min_param_change!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # maxdelta=0.3 - max percentage change (ex: 0,25)
    # mindelta=0.005 - min percentage change (ex: 0,01)
    # fct=LL_10
    # target=0

    LLAlpha_vis_last=fct(VECTPAR)
    LL_last=LLAlpha_vis_last
    VECTPAR_elite=cbind(t(VECTPAR),0)
    ii=3     # iterations
    jj=0     # nbr of calls
    lengthVECTPAR=length(VECTPAR)
    changes<-array(0.2, c(1,lengthVECTPAR))
    VECThist=array(0,c(5000,lengthVECTPAR))
    LLAlphahist_up<- array(0,c(1,lengthVECTPAR))
    LLAlphahist_down<- array(0,c(1,lengthVECTPAR))
    Lvisible=array(0,c(1,it))
    Lvisible[ii]=-1
    Lvisible[ii-1]=-0.5
    nbr_fcalls=0

    min_param_change<-cbind(array(VECTPAR[1]/300,c(1,lengthTheta)),array(VECTPAR[lengthTheta+1]/200,c(1,k)), array(0.01,c(1,lengthVECTPAR-lengthTheta-k)))

    while(nbr_fcalls<it){
        LLAlpha_vis_last<-LL_last


        VECTPAR_up=VECTPAR
        VECTPAR_down=VECTPAR

        for (m in 1:(lengthVECTPAR)){


            if(abs(VECTPAR[m])<0.05 & m>lengthTheta){
                VECTPAR[m]=-VECTPAR[m]                    # if parameter too small AND IS NOT a VARIANCE PARAMETER, try NEGATIVE value
            }


            VECTPAR_up[m]=VECTPAR[m]+max(VECTPAR[m]*changes[m],min_param_change[m])                                                       # positive change
            VECTPAR_down[m]=VECTPAR[m]-max(VECTPAR[m]*changes[m],min_param_change[m])

            LL_up=fct(VECTPAR_up)     # vecteur de L pour chaque param modifi? vers le haut
            jj=jj+1
            if(LL_up>LL_last){
                LL_last=LL_up
                VECTPAR[m]=VECTPAR_up[m]
                changes[m]=min(changes[m]*2,maxdelta)
            } else{
                LL_down=fct(VECTPAR_down)       # vecteur de L pour chaque param modifi? vers le bas
                jj=jj+1
                if(LL_down>LL_last){
                    LL_last=LL_down
                    VECTPAR[m]=VECTPAR_down[m]
                    changes[m]=min(changes[m]*2,maxdelta)
                } else{
                    changes[m]<-max(changes[m]/3,mindelta)
                }
            }




        }

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

        Lvisible[ii]=LL_last
        print(LL_last)
        print(VECTPAR)


        if(Lvisible[ii]==Lvisible[ii-1]){VECTPAR_elite=rbind(VECTPAR_elite,cbind(t(VECTPAR),ii)) # save the last VECTPAR before the "shake" AND the nbr of its corresponding function call!!!
        VECTPAR=VECTPAR+rnorm(length(VECTPAR))*0.3*VECTPAR
        VECTPAR[1:lengthTheta]=abs(VECTPAR[1:lengthTheta])
        LL_last=fct(VECTPAR)
        jj=jj+1 }

        ii=ii+1
        nbr_fcalls=jj
        if ((ii/10)%%1==0){    # si l'iteration est divisable par 10 sans d?cimales (chaque 10e iteration on note VECTPAR)
            VECThist[ii/10,]=VECTPAR
        }
        # print(changes)



    }
    # return=list(max(Lvisible[Lvisible<(-0.5)]), VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  which(Lvisible==max(Lvisible[Lvisible<(-0.5)])) ), 1:length(VECTPAR) ], jj, nbr_fcalls)  # car le 2e est tjs =-0.5
    envir$LL<-max(Lvisible[Lvisible<(-0.5)])
    envir$param<-if(length(as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] )  )==0){as.array(VECTPAR)} else{as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] ) }
    envir$iter<-ii-3
    envir$fct_calls<- nbr_fcalls

}
