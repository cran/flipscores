
flipscores <- function(model0, model1, X1, alternative = "two.sided",  w=1E5, scoretype="basic"){

  if(!missing(X1)){.X1 <- X1}


  if(missing(X1) & missing(model1)){
    stop('specify either X1 or model1')
  }


  if(!missing(model1)){
    if ( length(model1$coefficients) != 1+length(model0$coefficients)  ){
      stop('model1 should contain exactly one variable more than model0')
    }
  }


  #Get covariate of interest X1:
  if(missing(X1)){
    nvars <- dim(model.matrix(model1))[2]

    vecs0 <- vector("list", nvars-1)
    vecs1 <- vector("list", nvars)

    for(i in 1:(nvars-1) ){
      vecs0[[i]]<- model.matrix(model0)[,i]
      vecs1[[i]]<- model.matrix(model1)[,i]
    }

    vecs1[[nvars]]<- model.matrix(model1)[,nvars]


    .X1 = setdiff(vecs1,vecs0)[[1]]
  }


  if( length(levels(.X1))>2 ){
    stop('The covariate of interest has more than 2 levels, so the score is not defined.')
  }


  if(is.factor(.X1)){
    .X1 <- 2*(as.numeric(.X1)-1.5)
  }


  n <- length(.X1)


  if(class(model0)[1]=="lm"){
    a <- 1
  } else{

    #The following lines for obtaining 'a' are taken from the mdscore package on CRAN
    #by Antonio Hermes M. da Silva-Junior, Damiao N. da Silva and Silvia L. P. Ferrari


    mu.est <- model0$fitted.values
    eta.est <- model0$family$linkfun(mu.est)




    V <- if(model0$family[[1]] == "gaussian") quote(1) else
      as.list(model0$family$variance)[[2]]


    if(model0$family[[2]] %in% c("log", "cloglog", "logit")){
      mu <- switch(model0$family[[2]],
                   log     = quote(exp(eta)),
                   cloglog = quote(1 - exp(-exp(eta))),
                   logit   = quote(exp(eta)/(1 + exp(eta))))
    }else mu <- as.list(model0$family$linkinv)[[2]]

    Dmu <- D(mu,"eta")

    a <- eval(V, list(mu= mu.est)) / eval(Dmu, list(eta= eta.est))
  }




  #BASIC SCORE
  if(scoretype=="basic"){
    scores <- .X1*(model0$residuals)/a
    flips <- matrix( (rbinom(n*w, 1, 0.5))*2-1, nrow=w,ncol=n )
    flips[1,] <- numeric(n)+1
    flipScores <-  flips %*% scores
    pv <- sum(flipScores >= flipScores[1]) /w


    if(alternative=="greater" | alternative== "larger"){
      pv <- sum(flipScores >= flipScores[1]) /w
    }
    if(alternative=="less" | alternative== "smaller"){
      pv <- sum(flipScores <= flipScores[1]) /w
    }
    if(alternative=="two.sided"){
      pv <- 2 * min(  sum(flipScores >= flipScores[1]),
                      sum(flipScores <= flipScores[1])) / w
    }

  }


  #EFFECTIVE SCORE
  if(scoretype=="eff" | scoretype=="effective"){
    X <- cbind(.X1,model.matrix(model0))
    if(class(model0)[1]=="lm"){
      wei <- numeric(n)+1
    } else{
    wei <- as.numeric(model0$weights)
    }
    invInfMat=solve(t(X*wei)%*%X,silent=TRUE)

    scores <- X*(model0$residuals)/a


    effScores=invInfMat%*% t(scores)
    effScores=as.vector(effScores[1,]) # see Marohn 2002


    flips <- matrix( (rbinom(n*w, 1, 0.5))*2-1, nrow=w,ncol=n )
    flips[1,] <- numeric(n)+1
    flipEffScores <-  flips %*% effScores




    if(alternative=="greater" | alternative== "larger"){
      pv <- sum(flipEffScores >= flipEffScores[1]) /w
    }
    if(alternative=="less" | alternative== "smaller"){
      pv <- sum(flipEffScores <= flipEffScores[1]) /w
    }
    if(alternative=="two.sided"){
      pv <- 2 * min(  sum(flipEffScores >= flipEffScores[1]),
                     sum(flipEffScores <= flipEffScores[1])) / w
    }

  }

  pv
}





