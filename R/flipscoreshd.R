
flipscoreshd <- function(model0, X1,  w=200, scoretype="basic", V="identity"){
  
  
  if(missing(model0)){
    stop('specify model0')
  }
  
  if(missing(X1)){
    stop('specify X1')
  }
  
  if(dim(X1)[1]!=dim(model0$x)[1]){
    stop('the nr of rows of X1 should be n, the nr of rows of the design matrix of model0.')
  }
  
  if((!is.matrix(V))&(!is.character(V))){
    stop("V should be of type matrix or character")
  }

  n = dim(model0$x)[1]
  
  ###namesbeta <- setdiff(colnames(model1$x),colnames(model0$x))  #names of covariates of interest
  namesnuis <- colnames(model0$x)                                 #names of all other (nuisance) covs
  
  d <- dim(X1)[2]
  
  X1 <- as.matrix(X1)
  X.full <- as.matrix(cbind(X1,model0$x))  #full design matrix with tested and nuisance covariates
  
  if(is.matrix(V)){
    if((dim(V)[1]!=d)|(dim(V)[1]!=d)){
      stop("V should be a d by d matrix, where d is the dimension of the parameter of interest.")
    }
  }
  
  for(i in 1:d){
    if( length(levels(X1[,i]))>2 ){
      stop('A (categorical) covariate of interest has more than 2 levels, so its score is not defined.')
    }
  }
  
  for(i in 1:d){
    if(is.factor(X1[,i])){
      X1 <- 2*(as.numeric(X1)-1.5)
    }
  }

  if(class(model0)[1]=="lm"){
    a <- 1
  } else{
    
    #The following lines for obtaining 'a' are taken from the mdscore package on CRAN
    #by Antonio Hermes M. da Silva-Junior, Damiao N. da Silva and Silvia L. P. Ferrari
    
    
    mu.est <- model0$fitted.values
    eta.est <- model0$family$linkfun(mu.est)
    
    
    
    
    .V <- if(model0$family[[1]] == "gaussian") quote(1) else
      as.list(model0$family$variance)[[2]]
    
    
    if(model0$family[[2]] %in% c("log", "cloglog", "logit")){
      mu <- switch(model0$family[[2]],
                   log     = quote(exp(eta)),
                   cloglog = quote(1 - exp(-exp(eta))),
                   logit   = quote(exp(eta)/(1 + exp(eta))))
    }else mu <- as.list(model0$family$linkinv)[[2]]
    
    Dmu <- D(mu,"eta")
    
    a <- eval(.V, list(mu= mu.est)) / eval(Dmu, list(eta= eta.est))
  }
  
  
  
  if(is.character(V)){
    if(V=="identity"|V=="id") V=diag(d)
  }
  
  if(is.character(V)){
    if(V=="invinfo"){  #V is inverse of effective Fisher info
      if(class(model0)[1]=="lm"){
        wei <- numeric(n)+1
      } else{
        wei <- as.numeric(model0$weights)
      }
      infMat <- t(X.full*wei)%*%X.full  
      I11 <- infMat[1:d,1:d]
      I12 <- infMat[namesnuis,1:d]
      I22 <- infMat[namesnuis,namesnuis]
      effI <- I11 - t(I12) %*% solve(I22) %*% I12  #formula for the effective Fisher info
      inv.effI=try(solve(effI),silent=TRUE) 
      inv.effI=if (is(inv.effI,"try-error")) diag(dim(effI)[1]) else inv.effI
      V=inv.effI #V is the inverse of the effective information
    }
  }
  
  
  #BASIC SCORE
  if(scoretype=="basic"){
    scores <- X1*(residuals(model0,"response"))/a   #n by d matrix
    flips <- matrix( (rbinom(n*w, 1, 0.5))*2-1, nrow=w,ncol=n ) #w by n matrix
    flips[1,] <- numeric(n)+1
    flipScores<-  flips%*%scores   #w by d matrix
    
    T<-numeric(w)
    for(j in 1:w){
      T[j] <- t(flipScores[j,])%*%V%*%flipScores[j,]
    }
    
    pv <- sum(T[j]<=T)/w  #p-value
    
  }
  
  
  #EFFECTIVE SCORE
  if(scoretype=="eff" | scoretype=="effective"){
    
    #Compute I12 and I22 (needed for effective score)
    if(class(model0)[1]=="lm"){
      wei <- numeric(n)+1
    } else{
      wei <- as.numeric(model0$weights)
    }
    infMat <- t(X.full*wei)%*%X.full  
    I12 <- infMat[namesnuis,1:d]
    I22 <- infMat[namesnuis,namesnuis]
    
    #Get the summands underlying the effective scores. 
    #By using the following formula (Marohn, 2002) we don't need to invert the whole info matrix. Hence beta can be high-dimensional.
    nu <- t(X1*(residuals(model0,"response"))/a)   -t(I12)%*%solve(I22)%*%t(model0$x*(residuals(model0,"response"))/a)  #d by n matrix
    
    if(dim(nu)[2]==1){   
      nu<-t(nu)
    }
    
    #Flip and compute test statistics:
    flippedtt <- numeric(w)
    S <- numeric(d)
    flip.nu <- matrix(nrow=d,ncol=n)
    for(k in 1:d){S[k]<- sum(nu[k,])}
    flippedtt[1] <-   as.numeric(t(S)%*%V%*%S)       #the original, unflipped statistic
    for(j in 2:w){
      for(i in 1:n){
        flip.nu[,i] <- nu[,i]*(2*rbinom(1,1,0.5)-1)  #flip the nu's
      }
      for(k in 1:d){S[k]<- sum(flip.nu[k,])}         #score after flipping
      flippedtt[j] <-   as.numeric(t(S)%*%V%*%S)     #statistic after flipping
    }
    pv <- sum(flippedtt[1]<=flippedtt)/w   #p-value
  }
  
  pv  
}


