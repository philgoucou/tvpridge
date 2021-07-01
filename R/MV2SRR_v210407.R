#library(pracma)
#library(fGarch)


tvp.ridge <-function(X,Y,lambda.candidates=exp(linspace(-6,20,n=15)),oosX=c(),
                      lambda2=.1,kfold=5,CV.plot=TRUE,CV.2SRR=TRUE,
                      sig.u.param=.75,sig.eps.param=.75,ols.prior=0){
  
  ###############################################
  ############### translation ####################
  lambdavec=lambda.candidates
  #CV.2SRR=MV.CV.again
  CV.again=CV.2SRR
  #type=I(!plain.ridge)+1
  olsprior=ols.prior
  silent = as.numeric(!CV.plot)
  sweigths=1
  sv.param=sig.eps.param
  homo.param=sig.u.param
  
  #################################################
  ########### legend ##############################
  
  #type: plain ridge vs 2SRR
  # X,Y: matrices of data. Y can either be a vector or a matrix. But it should be formatted as a matrix.
  #lambdavec: candidates for ridge'lambda driving the amount of time-variation
  #,lambda2=.1,
  #CV.again: for 2SRR, do we CV again after updating weights.
  # sweigths=1,
  #oosX: a vector of X for out-of-sample prediction
  #kfold : number of folds for CV
  #silent, if 1, no graphs shows up, if 0, some do
  #olsprior, if =1, then shrinks beta_0's to OLS rather than zero
  #CV.2SRR : when doing multivariate, should we do the 2nd CV step for each equation separately (can take time)
  #homo.param: how much should we shrink to constant variance across u's (between 0 and 1, 0 means constant variance, 1 means plain 2-step)
  #sv.param: how much should we shrink to constant variance through time (between 0 and 1, 0 means constant variance, 1 means plain 2-step)
  
  ##################################################
  ##################################################
  
  
  #TYPE=type

  ######### Data pre-processing ###################

  mse<-c() #
  grrats<-list()
  grra<- list()
  fcast<- c()

  #Standardize variances
  Y=as.matrix(Y)
  Yoriginal=as.matrix(Y)

  M=dim(Y)[2]
  bigT=dim(Y)[1]
  if(is.null(M)==TRUE){M=1}
  scalingfactor = array(1,dim=c(M,dim(X)[2]))
  sdy = 1:M

  for(j in 1:dim(X)[2]){
    for(m in 1:M){scalingfactor[m,j]=sd(Y[,m])/sd(X[,j])}
    X[,j]=(X[,j])/sd(X[,j])
  }
  for(m in 1:M){sdy[m]=sd(Y[,m])
  Y[,m] = Y[,m]/sdy[m]}

  ZZ <- Zfun(X) #Basis expansion
  YY <- as.matrix(Y) #univariate model
  yy = vec(YY)
  dimX= dim(as.matrix(X))[2]+1 # + constant

    ######### TYPE 1 (FOR THE BASE MODEL OR ESTIMATION GIVEN HPS) ############
    #CV
    BETAS_GRR = array(0,dim=c(M,dimX,dim(Y)[1]))
    BETAS_GRRATS =BETAS_GRR
    LAMBDAS = rep(0,M)
    EWmat = array(0,dim=dim(Y))
    YHAT_VAR=Y
    YHAT_VARF=Y



    #if(CV.2SRR==TRUE){
      if(length(lambdavec)>1){
        l.subset = 1:length(lambdavec) #c(round(length(lambdavec)/4),round(2*length(lambdavec)/4),round(3*length(lambdavec)/4))
        cvlist <- CV.KF.MV(ZZ=ZZ,Y=Y,dimX = dimX,lambdavec = lambdavec[l.subset],lambda2=lambda2,
                           k=kfold,plot=abs(silent-1),sweigths=sweigths) #,eweigths = abs(yy)/mean(abs(yy)))
        lambdas_list =cvlist$lambdas_het
      }
    else{
      lambdas_list = rep(lambdavec,M)
    }
    #}
    #else{
    #  lambdas_list = rep(lambdavec[round(length(lambdavec)/2)],M)
    #  mse=100
    #}

    for(m in 1:M){
      if(M>1){print(paste('Computing/tuning equation ',m,sep=''))}
      #lambda1 <- cvlist$minimizer #lambdavec[round(length(lambdavec)/4)]
      lambda1 <-  lambdas_list[m]
      mse <- 1000 #cvlist$minima

      #Final estimation
      grr<-dualGRR(ZZ,Y[,m],dimX,lambda1,olsprior = olsprior,lambda2=lambda2,sweigths=sweigths) #,eweigths = sqrt(abs(Y[,m]))/mean(abs(Y[,m])))

      BETAS_GRR[m,,] <- grr$betas_grr[1,,]
      LAMBDAS[m]=lambda1
      YHAT_VAR[,m]=grr$yhat

      e = Y[,m] - grr$yhat
      e[abs(e)>50*mad(e)]=0 #glicth detector
      #plot.ts(e)
      garchstuff = hush(garchFit(data=e,~ garch(1,1)))
      EW = (garchstuff@sigma.t)^sv.param
      EW = (EW/mean(EW)) #technically, both .75 should be CV-ed values, the latter should be faster to implement
      EWmat[,m] = EW

      #NEW ###############
      betas_grr = BETAS_GRR[m,,]
      umat = betas_grr[,2:nrow(ZZ)]-betas_grr[,1:(nrow(ZZ)-1)]
      sigmasq= diag(umat%*%t(umat))^homo.param
      sigmasq=(sigmasq/mean(sigmasq))

      if(CV.again){cvlist <- cv.univariate(ZZ=ZZ,Y=Y[,m],dimX = dimX,lambda2 = lambda2,eweigths = EW,
                                                     lambdavec = lambdavec,k=kfold,plot=abs(silent-1),sweigths=sigmasq)
      usethislambda1=cvlist$minimizer
      }
      else{usethislambda1=lambda1}

      if(homo.param>0 | sv.param>0){
      grrats<-dualGRR(ZZ,Y[,m],dimX,lambda1=usethislambda1,eweigths = EW,olsprior = olsprior,
                      lambda2=lambda2,sweigths = sigmasq)
      betas_grrats <- grrats$betas_grr
      BETAS_GRRATS[m,,]=grrats$betas_grr[1,,]
      YHAT_VARF[,m]=grrats$yhat

      }else{
        betas_grrats <- grr$betas_grr
        BETAS_GRRATS[m,,]=grr$betas_grr[1,,]
        YHAT_VARF[,m]=YHAT_VAR[,m]
      }


    }
    EWvec = vec(EWmat)/mean(EWmat)

  lambda1_step2 =lambda1

  fcast <- c()
  BETAS_VARF=BETAS_GRRATS
  BETAS_VARF_STD=BETAS_VARF

  #re-scale
  for(m in 1:M){for(j in 1:(dimX-1)){
    BETAS_VARF[m,j+1,]=BETAS_VARF[m,j+1,]*scalingfactor[m,j]
    BETAS_GRR[m,j+1,]=BETAS_GRR[m,j+1,]*scalingfactor[m,j]
  }
    BETAS_VARF[m,1,]=BETAS_VARF[m,1,]*sdy[m]
    BETAS_GRR[m,1,]=BETAS_GRR[m,1,]*sdy[m]
    YHAT_VARF[,m]= YHAT_VARF[,m]*sdy[m]
    YHAT_VAR[,m]= YHAT_VAR[,m]*sdy[m]

  }

  if(length(oosX)>1){
    fcast <-  BETAS_VARF[,,dim(BETAS_VARF)[3]] %*% as.matrix(append(1,oosX))  #basically the last beta_t
  }

  closeAllConnections()

  return(list(betas.rr=BETAS_GRR,betas.2srr=BETAS_VARF,
              lambdas=LAMBDAS,forecast=fcast, 
              yhat.rr=YHAT_VAR,yhat.2srr=YHAT_VARF,sig.eps=EWvec))
}


#to silence the GARCH function
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}



#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################



dualGRR <-function(Zprime,y,dimX,lambda1,lambda2=0.001,CI=0,sweigths=1,olsprior=1,eweigths=1,GCV=0,calcul_beta=1,nf=0){
  ### Inputs are: ####
  #Zprime is the expanded X's
  #y is self-explanatory
  #lambda[1] is the the ratio of sigma_u and sigma_epsilon
  #lambda[2] is the penalty for non-u's parameters
  # d is the derivative being penalized (1 or 2)
  #####################


  ############ PREPARE PRIOR MATRIX ############
  lambdas = c(lambda1,lambda2)
  if(nf==0){nf=dimX}

  #lambda[2]  cannot be too small. let's impose a safety guard.
  if(lambdas[2]<0.0000001){
    print('Lambda_2 imposed to be at least 0.0001')
    lambdas[2]<-0.0000001
  }

  #create inv(Lambda_p) implicitely
  ncolZ = ncol(Zprime)
  T=(dim(Zprime)[2]-dimX)/nf +1

  if(length(sweigths)==1){sweigths=rep(1,1,nf)}
  Kmat_half<-t(Zprime)
  for(m in 1:(nf)){
    begin=(m-1)*(T-1)+1
    end=(m)*(T-1)
    Kmat_half[begin:end,]= (1/lambdas[1])*sweigths[m]*Kmat_half[begin:end,]
    if(-nf>-dimX){Kmat_half[begin,]=1000*Kmat_half[begin,]}
  }
  Kmat_half[(ncolZ-dimX+1):ncolZ,]=(1/lambdas[2])*Kmat_half[(ncolZ-dimX+1):ncolZ,]

  ############ DUAL GRR ############

  Lambda_T<-diag(nrow(Zprime))
  if(length(eweigths)>1){Lambda_T<- diag(as.numeric(eweigths))}

  param=nf*(T-1) + dimX
  obs = nrow(Zprime)

  if(param>obs){ #Go Dual
    Kmat_du <- Zprime%*%Kmat_half

    if(olsprior==0){
      alpha <- solve(Kmat_du + Lambda_T,y)
      uhat <- Kmat_half%*%alpha
    }
    if(olsprior==1){
      X= Zprime[,(ncolZ-dimX+1):ncolZ]
      beta_ols<-solve(crossprod(X),crossprod(X,y))
      alpha <- solve(Kmat_du + Lambda_T,(y-X%*%beta_ols))
      uhat <- Kmat_half%*%alpha
      uhat[(ncolZ-dimX+1):ncolZ]= uhat[(ncolZ-dimX+1):ncolZ]+beta_ols
    }
    yhat <- Kmat_du%*%alpha
  }

  else{ #Go Primal
    for(tt in 1:obs){Kmat_half[,tt]=Kmat_half[,tt]*eweigths[tt]^-1}
    Kmat_pri <- Kmat_half%*%Zprime

    if(olsprior==0){
      uhat <- solve(Kmat_pri + diag(param),Kmat_half%*%y)
    }
    if(olsprior==1){
      X= Zprime[,(ncolZ-dimX+1):ncolZ]
      beta_ols<-solve(crossprod(X),crossprod(X,y))
      uhat <- solve(Kmat_pri + diag(param),Kmat_pri%*%(y-X%*%beta_ols))
      uhat[(ncolZ-dimX+1):ncolZ]= uhat[(ncolZ-dimX+1):ncolZ]+beta_ols
    }
    yhat <- Zprime%*%uhat
  }

  ############ RECOVER BETA'S ############
  betas_grr <- array(0,dim=c(dim(as.matrix(y))[2],nf,T))

  if(calcul_beta==1){

    for(eq in 1:dim(as.matrix(y))[2]){
      for(k in 1:nf){betas_grr[eq,k,1] = uhat[(dim(uhat)[1]-dimX+k),eq]}
      for(t in 2:(T)){
        begin=(k-1)*(T-1)+1
        positions=c()
        for(k in 1:nf){positions = append(positions,(k-1)*(T-1)+(t-1))}
        #print(end)
        #if(d==1){begin=end}
        betas_grr[eq,,t] = betas_grr[eq,,t-1]+uhat[positions,eq]
      }
    }
  }
  ############ CI'S ############
  if(CI>0){

    Kmh = Kmat_half
    if(param>obs){ #Go Dual
      for(tt in 1:obs){Kmh[,tt]=Kmat_half[,tt]*eweigths[tt]^-1}
      #Kmat_pri <- Kmh%*%Zprime
    }

    if(dimX>nf){
      Kmh=Kmh[1:(dim(Kmat_half)[1]-dimX),]
      #Kmat_pri = Kmat_pri[1:(dim(Kmat_half)[1]-dimX),1:(dim(Kmat_half)[1]-dimX)]
    }
    else{

      #Putting back Kmat_half into order
      #Put back the first col of each
      for(k in 1:dimX){
        Kmh[1+T*(k-1),] = Kmh[(dim(Kmat_half)[1]-dimX+k),]
        begin = 2+T*(k-1)
        end = T+T*(k-1)
        Kmh[begin:end,]= Kmh[(begin-k):(end-k),]
      }}
    #Getting #Vb

    invCov = chol2inv(chol(tcrossprod(Kmh)+diag(dim(Kmh)[1])))

    CT= array(1,dim=c(T-1,T-1))
    CT= lower.triangle(CT)
    C=kronecker(diag(nf),CT)
    #print(dim(invCov))
    #print(dim(C))



    #Vb=diag(C%*%tcrossprod(invCov,C))*sd(y-yhat)/sqrt(dim(Kmh)[2])
    #sandwich = invCov%*%Kmat_pri%*%invCov
    sandwich=invCov*(sd(y-yhat)^2)*(lambdas[1])^(-1)
    Vb=sqrt(diag(C%*%tcrossprod(sandwich,C))) #/(dim(Kmh)[2])

    betas_grr_ci <- array(0,dim=c(dim(as.matrix(y))[2],nf,T,2))
    for(c in 1:2){
      #temp<- uhat + ((-1)^c)*1.96*Vu
      for(eq in 1:dim(as.matrix(y))[2]){
        for(k in 1:nf){ #variables
          #betas_grr_ci[eq,k,1,c] = temp[(dim(uhat)[1]-dimX+k),eq]
          for(t in 1:(T)){
            #begin=(k-1)*(T-1)+1
            #end =(k-1)*(T-1)+(t-1)
            betas_grr_ci[eq,k,t,c] = betas_grr[eq,k,t]+((-1)^c)*CI*Vb[(k-1)*(T-1)+1]
          }
        }
      }
    }

    return(list(uhat=uhat,betas_grr=betas_grr,yhat=yhat,betas_grr_ci=betas_grr_ci,lambdas=lambdas,
                sigmasq=sweigths))
  }
  else{return(list(uhat=uhat,betas_grr=betas_grr,yhat=yhat,sigmasq=sweigths,lambdas=lambdas))
  }
}




Zfun <- function(data){

  #Generate matrices
  X= cbind(matrix(1,(nrow(data)),1),data)
  Z=array(0,dim=c((nrow(data)),(nrow(data)-1),(ncol(data)+1)))
  X<-as.matrix(X)

  #The sequence of tau's

  #New regressors
  for(tt in 2:(dim(data)[1])){
    Z[tt,1:(tt-1),]= repmat(X[tt,],tt-1,1)
  }

  Zprime <- c()
  for(kk in 1:(ncol(data)+1)){
    Zprime = cbind(Zprime,Z[,,kk])
  }

  Zprime <- cbind(Zprime,X)


  return(Zprime)
}

cv.univariate <- function(Y,ZZ,k,lambdavec,lambda2=0.001,dimX,plot=0,sweigths=1,eweigths=1,nf=dimX){
  #what takes time is to compute the K matrix each time. Using the FWL trick
  #and having the k loop being the outter loop rather than lambda (as GA) function,
  #we can speed things up dramatically. However, CV interpretation will change.
  #Instead of finding the best model, it will find the best improvment wrt to constant
  #paramater VAR.

  set.seed(1071) # Seed for replicability


  # k - number of folds
  Y<- as.matrix(Y)
  T<-dim(ZZ)[1]
  data <- cbind(Y,ZZ)
  dfTrain<- data

  index <- sample(1:k, nrow(data), replace = TRUE)

  #randomseq <- sample(nrow(data))
  #dfTrain <- data[randomseq,]
  #folds <- cut( seq(1,nrow(data)), breaks = k, labels=FALSE)

  seqY = 1:dim(Y)[2] #for MV model
  PMSE <- array(0,dim=c(k,length(lambdavec)))

  # Do the cross validation for K fold
  for(i in 1:k){
    #Segment your data by fold using the which() function
    testIndexes <- which(index == i, arr.ind=TRUE)
    testData <- dfTrain[testIndexes, ]
    trainData <- dfTrain[-testIndexes, ]

    Lambda_T<-diag(nrow(trainData))
    if(length(eweigths)>1){Lambda_T<- diag(as.numeric((eweigths[-testIndexes])))}

    #Correction for dropouts (ineficient for now)
    bindex = index
    bindex[index==i]=0
    bindex[bindex!=0]=1
    DOsigfactor <- (1+cumul_zeros(bindex))

    #Kmat (Train)
    Z = trainData[,-seqY]
    ncolZ = ncol(Z)
    X = Z[,(ncolZ-dimX+1):ncolZ]
    MX = diag(nrow(X)) - X%*%solve(crossprod(X)+lambda2*diag(ncol(X)))%*%t(X)
    MXZ=MX%*%Z

    for(m in 1:nf){
      begin=(m-1)*(T-1)+1
      end=(m)*(T-1)
      if(length(sweigths)>1){MXZ[,begin:end]= (sweigths[m])*MXZ[,begin:end]}
      for(tt in 1:(T-1)){
        MXZ[,(m-1)*(T-1)+tt]= MXZ[,(m-1)*(T-1)+tt]*DOsigfactor[tt+1]
      }
    }

    Kmat = MXZ%*%t(MXZ)

    #kmat (for TEST)
    z = testData[,-seqY]
    ncolz = ncol(z)
    x = z[,(ncolz-dimX+1):ncolz]
    Mx = diag(nrow(x)) - x%*%solve(crossprod(x)+lambda2*diag(ncol(x)))%*%t(x)
    mxz=Mx%*%z
    if(length(sweigths)==1){kmat = mxz%*%t(MXZ)}
    else{
      for(m in 1:length(sweigths)){
        begin=(m-1)*(T-1)+1
        end=(m)*(T-1)
        mxz[,begin:end]= sweigths[m]*mxz[,begin:end]
      }
      kmat = mxz%*%t(MXZ)
    }

    #OOS pred
    for(j in 1:length(lambdavec)){
      #print(j)
      pred <- kmat%*%solve(Kmat+lambdavec[j]*Lambda_T,MX%*%trainData[,seqY])
      PMSE[i,j] <- mean((pred-Mx%*%testData[,1:dim(Y)[2]])^2)
    }
  }

  #Find min
  score <- colMeans(PMSE)
  lambdastarpos <- which(score == min(score))
  finalmse <- score[lambdastarpos]
  lambdastar <- lambdavec[lambdastarpos]

  if(length(lambdavec)>1){
    #one SE rule
    SE = array(NA,dim=c(1,length(lambdavec)))
    for(j in 1:length(lambdavec)){
      SE[1,j]=sd(PMSE[,j])/sqrt(k)
    }
    se <- SE[lambdastarpos] #mean(SE)
    scoreub = score + as.numeric(SE)
    scorelb = score - as.numeric(SE)

    lambda1sepos <- lambdastarpos
    repeat {
      lambda1sepos <- lambda1sepos + 1
      if(lambda1sepos>=length(score)) {break}
      if(score[lambda1sepos]>(finalmse+se)) {break}
    }
    lambda1se <- lambdavec[lambda1sepos]

    lambda2sepos <- lambdastarpos
    repeat {
      lambda2sepos <- lambda2sepos + 1
      if(lambda2sepos>=length(score)) {break}
      if(score[lambda2sepos]>(finalmse+2*se)) {break}
    }
    lambda2se <- lambdavec[lambda2sepos]

    #Plot
    if(plot==1){
      limit2 <- lambda2sepos
      repeat {
        limit2 <- limit2 + 1
        if(limit2>=length(score)) {break}
        if(score[limit2]>(finalmse+20*se)) {break}
      }
      limit1 <- lambdastarpos
      repeat {
        limit1 <- limit1 - 1
        if(limit1<=1) {break}
        if(score[limit1]>(finalmse+20*se)) {break}
      }
      ts.plot(score[limit1:limit2],lwd=4,xlab='lambda') #,scoreub,scorelb)
      abline(h=finalmse+se, col="blue",lwd=4)
      abline(h=finalmse+2*se, col="purple",lwd=4)
      title('Cross-Validation Results')
      legend("topright",  legend=c('MSE','Minimum + 1 se','Minimum + 2 se'),
             col=c('black','blue','purple') ,lty=1,lwd=4)
    }
    return(list(minimizer=lambdastar,minima=finalmse,
                minimizer1se=lambda1se,minimizer2se=lambda2se))
  }
  else{return(list(minimizer=lambdastar,minima=finalmse))}
}


cumul_zeros <- function(x)  {
  x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}


CV.KF.MV <- function(Y,ZZ,k,lambdavec,lambda2=0.001,dimX,plot=0,sweigths=1,eweigths=1,nf=dimX,beta0_given=FALSE){
  #what takes time is to compute the K matrix each time. Using the FWL trick
  #and having the k loop being the outter loop rather than lambda (as GA) function,
  #we can speed things up dramatically. However, CV interpretation will change.
  #Instead of finding the best model, it will find the best improvment wrt to constant
  #paramater VAR.

  set.seed(1071) # Seed for replicability

  # k - number of folds
  Y<- as.matrix(Y)
  T<-(dim(ZZ)[2]-dimX)/nf +1
  if(beta0_given==TRUE){ZZ=ZZ[,1:(dim(ZZ)[2]-dimX)]}
  data <- cbind(Y,ZZ)
  dfTrain<- data

  index <- sample(1:k, nrow(data), replace = TRUE)

  #randomseq <- sample(nrow(data))
  #dfTrain <- data[randomseq,]
  #folds <- cut( seq(1,nrow(data)), breaks = k, labels=FALSE)

  seqY = 1:dim(Y)[2] #for MV model
  PMSE <- array(0,dim=c(k,length(lambdavec)))
  PMSE.het <- array(0,dim=c(k,length(lambdavec),length(seqY)))

  # Do the cross validation for K fold
  for(i in 1:k){
    #Segment your data by fold using the which() function
    testIndexes <- which(index == i, arr.ind=TRUE)
    testData <- dfTrain[testIndexes, ]
    trainData <- dfTrain[-testIndexes, ]

    Lambda_T<-diag(nrow(trainData))
    if(length(eweigths)>1){Lambda_T<- diag(as.numeric((eweigths[-testIndexes])))}

    #Correction for dropouts (ineficient for now)
    bindex = index
    bindex[index==i]=0
    bindex[bindex!=0]=1
    DOsigfactor <- (1+cumul_zeros(bindex))
    #Kmat (Train)
    Z = trainData[,-seqY]
    ncolZ = ncol(Z)
    if(beta0_given==FALSE){
      X = Z[,(ncolZ-dimX+1):ncolZ]
      MX = diag(nrow(X)) - X%*%solve(crossprod(X)+lambda2*diag(ncol(X)))%*%t(X)
      MXZ=MX%*%Z
    }
    else{MXZ=Z} #This saves lotsa time

    for(m in 1:nf){
      begin=(m-1)*(T-1)+1
      end=(m)*(T-1)
      if(length(sweigths)>1){MXZ[,begin:end]= (sweigths[m])*MXZ[,begin:end]}
      if(-nf>-dimX){MXZ[,begin]= 1000*MXZ[,begin]}
      for(tt in 1:(T-1)){
        MXZ[,(m-1)*(T-1)+tt]= MXZ[,(m-1)*(T-1)+tt]*DOsigfactor[tt+1]
      }
    }

    param=nf*(T-1) #+ dimX
    obs = nrow(ZZ)

    if(param>obs){ #GO DUAL

      Kmat = tcrossprod(MXZ)

      #kmat (for TEST)
      z = testData[,-seqY]
      ncolz = ncol(z)
      if(beta0_given==FALSE){
        x = z[,(ncolz-dimX+1):ncolz]
        Mx = diag(nrow(x)) - x%*%solve(crossprod(x)+lambda2*diag(ncol(x)))%*%t(x)
        mxz=Mx%*%z
      }
      else{mxz=z}

      #if(length(sweigths)==1){kmat = tcrossprod(mxz,MXZ)}
      #else{
      #  for(m in 1:length(sweigths)){
      #    begin=(m-1)*(T-1)+1
      #    end=(m)*(T-1)
      #    mxz[,begin:end]= sweigths[m]*mxz[,begin:end]
      #  }
      #  kmat = tcrossprod(mxz,MXZ)
      #}
      for(m in 1:nf){
        begin=(m-1)*(T-1)+1
        if(-nf>-dimX){mxz[,begin]= 1000*mxz[,begin]}
      }

      if(length(sweigths)==1){kmat = tcrossprod(mxz,MXZ)}
      else{
        for(m in 1:nf){
          begin=(m-1)*(T-1)+1
          end=(m)*(T-1)
          mxz[,begin:end]= sweigths[m]*mxz[,begin:end]
        }
        kmat = tcrossprod(mxz,MXZ)
      }

      #OOS pred
      if(beta0_given==FALSE){
        MXY=MX%*%trainData[,seqY]
        real=Mx%*%testData[,1:dim(Y)[2]]
      }
      else{
        MXY=trainData[,seqY]
        real=testData[,1:dim(Y)[2]]
      }

      for(j in 1:length(lambdavec)){
        pred <- kmat%*%solve(Kmat+lambdavec[j]*Lambda_T,MXY)
        PMSE[i,j] <- mean((pred-real)^2)
        for(mm in 1:dim(Y)[2]){PMSE.het[i,j,mm]= mean((pred[,mm]-as.matrix(real)[,mm])^2)}
      }

    } # end dual
    else{ #GO PRIMAL

      MXZ2=MXZ
      if(length(eweigths)>1){for(tt in 1:nrow(Lambda_T)){MXZ2[tt,]=MXZ[tt,]*Lambda_T[tt,tt]^-1}}
      Kmat = crossprod(MXZ2,MXZ)


      #kmat (for TEST)
      z = testData[,-seqY]
      ncolz = ncol(z)
      if(beta0_given==FALSE){
        x = z[,(ncolz-dimX+1):ncolz]
        Mx = diag(nrow(x)) - x%*%solve(crossprod(x)+lambda2*diag(ncol(x)))%*%t(x)
        mxz=Mx%*%z
      }
      else{mxz=z}

      for(m in 1:nf){
        begin=(m-1)*(T-1)+1
        if(-nf>-dimX){mxz[,begin]= 1000*mxz[,begin]}
      }

      if(length(sweigths)>1){
        for(m in 1:nf){
          begin=(m-1)*(T-1)+1
          end=(m)*(T-1)
          mxz[,begin:end]= sweigths[m]*mxz[,begin:end]
        }
      }

      #OOS pred
      if(beta0_given==FALSE){
        MXY=crossprod(MXZ2,MX%*%trainData[,seqY])
        real=Mx%*%testData[,1:dim(Y)[2]]
      }
      else{
        MXY=crossprod(MXZ2,trainData[,seqY])
        real=testData[,1:dim(Y)[2]]
      }

      for(j in 1:length(lambdavec)){
        pred <- mxz%*%solve(Kmat+lambdavec[j]*diag(nrow(Kmat)),MXY)
        #pred=abs(fft(pred))
        #real=abs(fft(real))
        PMSE[i,j] <- mean((pred-real)^2)
        for(mm in 1:dim(Y)[2]){PMSE.het[i,j,mm]= mean((pred[,mm]-as.matrix(real)[,mm])^2)}
      }
    } # end primal
  }

  lambdas_het = seqY
  score.het = apply(PMSE.het,c(2,3),mean)
  for(mm in 1:dim(Y)[2]){lambdas_het[mm]=lambdavec[which(score.het[,mm]==min(score.het[,mm]))]}

  #for now
  #Find min
  score <- colMeans(PMSE)
  lambdastarpos <- which(score == min(score))
  finalmse <- score[lambdastarpos]
  lambdastar <- lambdavec[lambdastarpos]

  if(length(lambdavec)>1){
    #one SE rule
    SE = array(NA,dim=c(1,length(lambdavec)))
    for(j in 1:length(lambdavec)){
      SE[1,j]=sd(PMSE[,j])/sqrt(k)
    }
    se <- SE[lambdastarpos] #mean(SE)
    scoreub = score + as.numeric(SE)
    scorelb = score - as.numeric(SE)

    lambda1sepos <- lambdastarpos
    repeat {
      lambda1sepos <- lambda1sepos + 1
      if(lambda1sepos>=length(score)) {break}
      if(score[lambda1sepos]>(finalmse+se)) {break}
    }
    lambda1se <- lambdavec[lambda1sepos]

    lambda2sepos <- lambdastarpos
    repeat {
      lambda2sepos <- lambda2sepos + 1
      if(lambda2sepos>=length(score)) {break}
      if(score[lambda2sepos]>(finalmse+2*se)) {break}
    }
    lambda2se <- lambdavec[lambda2sepos]

    #Plot
    if(plot==1){
      limit2 <- lambda2sepos
      repeat {
        limit2 <- limit2 + 1
        if(limit2>=length(score)) {break}
        if(score[limit2]>(finalmse+20*se)) {break}
      }
      limit1 <- lambdastarpos
      repeat {
        limit1 <- limit1 - 1
        if(limit1<=1) {break}
        if(score[limit1]>(finalmse+20*se)) {break}
      }
      ts.plot(score[limit1:limit2],lwd=4,xlab='lambda') #,scoreub,scorelb)
      abline(h=finalmse+se, col="blue",lwd=4)
      abline(h=finalmse+2*se, col="purple",lwd=4)
      title('Cross-Validation Results')
      legend("topright",  legend=c('MSE','Minimum + 1 se','Minimum + 2 se'),
             col=c('black','blue','purple') ,lty=1,lwd=4)
      #abline(h=finalmse+2*se, col="purple")
    }
    return(list(minimizer=lambdastar,minima=finalmse,
                minimizer1se=lambda1se,minimizer2se=lambda2se,lambdas_het=lambdas_het))
  }
  else{return(list(minimizer=lambdastar,minima=finalmse,lambdas_het=lambdas_het))}
}


