#' Integrative finite population mean (IntegrativeFPM)
#'
#'Implements integrative analyses for the finite population mean parameters
#'combining a non-probability sample with a probability sample which provides
#'high-dimensional representative covariate information of the target population.
#'
#'The non-probability sample contains observations on (X,Y); however, the sampling mechanism is unknown.
#'On the other hand, although the probability sample with sampling weights represents
#'the finite population, it contains observations only on X.
#'
#'The function "IntegrativeFPM" implements a two-step approach for integrating the probability sample
#'and the nonprobability sample for finite population inference.
#'
#'The first step selects important variables in the sampling score model and
#'the outcome model by penalization.
#'
#'The second step conducts a doubly robust inference for the finite population
#'mean of interest. In this step, the nuisance model parameters are refitted based on the selected variables
#' by minimizing the asymptotic squared bias of the doubly robust estimator.
#' This estimating strategy mitigates the possible first-step selection error and
#' renders the doubly robust estimator root-n consistent if either the sampling probability
#' or the outcome model is correctly specified.
#'
#' @param y is a vector of outcome; NA for the probability sample (n x 1)
#' @param x is a matrix of covariates without intercept  (n x p)
#' @param deltaB is a vector of the binary indicator of belonging to the nonprobability sample;
#' i.e., 1 if the unit belongs to the nonprobability sample, and 0 otherwise  (n x 1)
#' @param sw is a vector of weights: 1 for the unit in the nonprobability sample and the design weight for the unit in the probability sample  (n x 1)
#' @param family specifies the outcome model
#'
#' \code{"gaussian"}: a linear regression model for the continuous outcome
#'
#' \code{"binomial"}: a logistic regression model for the binary outcome
#'
#' @param lambda_a is a scalar tuning parameter in the penalized estimating equation for alpha (the sampling score parameter)
#'
#' The sampling score is a logistic regression model for the probability of selection into the nonprobability sample given X
#'
#' @param lambda_b is a scalar tuning parameter in the penalized estimation for beta (the outcome model parameter)
#'
#' The outcome model is a linear regression model for continuous outcome or a logistic regression model for binary outcome
#'
#'
#'
#' @return
#' \itemize{
#'
#' \item \code{est}: estimate of the finite population mean of Y
#'
#' \item \code{se}:  standard error estimate for \code{est}
#'
#' \item \code{alpha.selected}: the column(s) of X selected for the sampling score model
#'
#' \item \code{beta.selected}: the column(s) of X selected for the outcome model
#' }
#'
#' @details
#'
#'
#' @import MASS rootSolve ncvreg stats
#'
#'@references
#'
#'Yang, S., Kim, J.K, and Song, R. (2019). Doubly Robust Inference when Combining Probability and Non-probability Samples with High-dimensional Data.
#'
#' @examples
#'
#'set.seed(1234)
#'
#'## population size
#'
#'N <- 10000
#'
#'## x is a p-dimensional covariate
#'
#'p <- 50
#'x <- matrix( rnorm(N*p,0,1),N,p)
#'
#'## y is a continuous outcome
#'
#'beta0 <- c(1,1,1,1,1,rep(0,p-4))
#'y <- cbind(1,x)%*%beta0 + rnorm(N,0,1)
#'true <- mean(y)
#'
#'## y2 is a binary outcome
#'
#'ly2 <- (cbind(1,x)%*%beta0)
#'ply <- exp(ly2)/(1+exp(ly2))
#'y2 <- rbinom(N,1,ply)
#'true2 <- mean(y2)
#'
#'## A.set is a prob sample: SRS
#'## sampling probability into A is known when estimation
#'
#'nAexp <- 1000
#'probA <- rep(nAexp/N,N)
#'A.index <- rbinom(N,size = 1,prob = probA)
#'A.loc <- which(A.index == 1)
#'nA <- sum(A.index == 1)
#'sw.A <- 1/probA[A.loc]
#'x.A <- x[A.loc,]
#'y.A <- rep(NA,nA) # y is not observed in Sample A
#'y2.A <- rep(NA,nA)
#'
#'## B.set is a nonprob sample
#'## sampling probability into B is unknown when estimation
#'
#'nBexp <- 2000
#'alpha0 <- c(-2,1,1,1,1,rep(0,p-4))
#'probB <- (1+exp(-cbind(1,x)%*%alpha0))^(-1)
#'B.index <- rbinom(N,size = 1,prob = probB)
#'B.loc <- which(B.index == 1)
#'nB <- sum(B.index)
#'x.B <- x[B.loc,]
#'y.B <- y[B.loc]
#'y2.B <- y2[B.loc]
#'
#'## combined dataset
#'
#'y.AB <- c(y.A,y.B)
#'y2.AB <- c(y2.A,y2.B)
#'x.AB <- rbind(x.A,x.B)
#'deltaB <- c(rep(0,nA),rep(1,nB))
#'sw <- c(sw.A,rep(1,nB))
#'
#'## specify tuning parameters
#'
#'lambda_a <- 1.25
#'lambda_b <- 0.25
#'lambda2_b <- 0.02
#'
#'true
#'IntegrativeFPM(y=y.AB, x=x.AB, deltaB, sw, family="gaussian",lambda_a,lambda_b)
#'
#'true2
#'IntegrativeFPM(y=y2.AB, x=x.AB, deltaB, sw, family="binomial",lambda_a, lambda2_b)
#'
#' @export

IntegrativeFPM <- function(y, x, deltaB, sw, family, lambda_a, lambda_b){
  x.AB<-x
  y.AB<-y
  y.AB[which(is.na(y.AB))]<-0
  nA<-sum(1-deltaB)
  nB<-sum(  deltaB)
  loc.A<-which(deltaB==0)
  loc.B<-which(deltaB==1)
  n<-dim(x.AB)[1]
  p<-dim(x.AB)[2]
  x.B<-x.AB[loc.B,]
  y.B<-y.AB[loc.B]
  N<-sum(sw[loc.A])

  ##################################################
  ## Step 1: variable selection
  ## penalized estimating equations for alpha and beta
  ##################################################

  alpha.ini<-rep(0,p+1)
  beta.ini <-rep(0,p+1)
  par0<-c(alpha.ini)
  En<-matrix(0,p+1,p+1)
  cnt1<-0
  epsilon=10^(-6)
  for(jj in 1:100){
    cnt1<-cnt1+1
    Uee0<-Uee_a(par0,deltaB,x.AB,y.AB,sw)
    Aee0<-Aee_a(par0,deltaB,x.AB,y.AB,sw)
    diag(En)<-abs(q1(par0,lambda_a))/(epsilon+abs(par0))
    par<-  par0+ ginv(Aee0+En)%*%(Uee0-En%*%par0)
    if( sum( abs(par-par0)  )<10^(-6) )break;
    if( sum( abs(par-par0)  )>1000 )break;
    par0<-par
  }
  par[which( abs(par) <10^(-3))]<-0
  alpha.est<-par
  alpha.selected<-which(alpha.est!=0)

  if(family=="gaussian"){

    beta <- ncvreg::ncvreg(x.B, y.B, lambda=lambda_b, penalty='SCAD', family="gaussian")
    beta.est <- beta$beta
    beta.selected<-as.numeric(which(beta.est!=0))

    ##################################################
    ## Step 2: joint ee
    ##################################################

    Cindex<- unique( c(beta.selected[-1]-1,alpha.selected[-1]-1) )
    pjee<-length(Cindex)
    x.ABC<-x.AB[,Cindex]

    par0<-rep(0,2*(pjee+1))
    jeepar<-rootSolve::multiroot(Uee,par0,deltaB=deltaB,x.AB=x.ABC,y.AB=y.AB,sw=sw)$root
    jeepar1<-jeepar[ 1:(pjee+1)]
    jeepar2<-jeepar[(pjee+2):(2*pjee+2)]

    lpi.est <- exp(as.vector(as.matrix( cbind(1,x.ABC) ) %*% as.matrix(jeepar1)))
    piB.est  <- lpi.est/(1+lpi.est)
    m.est <- as.vector( as.matrix( cbind(1,x.ABC) )%*% as.matrix(jeepar2))

    ##################################################
    ## doubly robust estimation for FPM
    ##################################################

    est.pdr    <- sum( (y.AB-m.est)*deltaB/piB.est )/N+(sum(m.est*(1-deltaB)*sw ) )/N
    sw.A<-sw[loc.A]
    m.A<-m.est[loc.A]
    sigmasqhat<-mean( (y.AB[loc.B]-m.est[loc.B] )^2)
    ve.pdr <- sum((sw.A^2 - sw.A) * (m.A)^2)/(N^2) + (sum(deltaB*(1-2*piB.est )/piB.est ^2) + (N )) * sigmasqhat/(N^2)
    se.pdr<-sqrt(ve.pdr)

  }

  if(family=="binomial"){

    beta <- ncvreg::ncvreg(x.B, y.B, lambda=lambda_b, penalty='SCAD', family="binomial")
    beta.est <- beta$beta
    beta.selected<-as.numeric(which(beta.est!=0))

    ##################################################
    ## Step 2: joint ee
    ##################################################

    Cindex<- unique(beta.selected[-1]-1,alpha.selected[-1]-1)
    pjee<-length(Cindex)
    x.ABC<-x.AB[,Cindex]

    par0<-rep(0,2*(pjee+1))
    jeepar<-rootSolve::multiroot(Uee_binary,par0,deltaB=deltaB,x.AB=x.ABC,y.AB=y.AB,sw=sw)$root
    jeepar1<-jeepar[ 1:(pjee+1)]
    jeepar2<-jeepar[(pjee+2):(2*pjee+2)]

    lpi.est <- exp(as.vector(as.matrix( cbind(1,x.ABC) ) %*% as.matrix(jeepar1)))
    piB.est  <- lpi.est/(1+lpi.est)
    lm.est <- as.vector( as.matrix( cbind(1,x.ABC) )%*% as.matrix(jeepar2))
    m.est  <- expoit(lm.est)

    ##################################################
    ## doubly robust estimation for FPM
    ##################################################

    est.pdr    <- sum( (y.AB-m.est)*deltaB/piB.est )/N+(sum(m.est*(1-deltaB)*sw ) )/N

    sw.A<-sw[loc.A]
    m.A<-m.est[loc.A]
    sigmasqhat.A<-m.est[loc.A]*(1-m.est[loc.A])
    infl1<-  (y.AB-m.est)^2*deltaB/(piB.est^2)
    infl2<-  (y.AB-m.est)^2*deltaB/piB.est
    ve.pdr <-sum((sw.A^2-sw.A)*(m.A)^2)/(N^2)+ (sum((infl1)-2*infl2)+sum( sw.A*sigmasqhat.A  )) /(N^2)
    se.pdr<-sqrt(ve.pdr)
  }

  alpha.selected<-alpha.selected[-1]-1
  beta.selected<-beta.selected[-1]-1
  return(list(est=est.pdr,se=se.pdr,alpha.selected=alpha.selected,beta.selected=beta.selected))
}


####################################################################################
## function list
####################################################################################


q1<-function(par,lambda){
  ## SCAD penalty derivative
  a=3.7
  penaltyd<-(abs(par)<lambda)*lambda +(abs(par)>=lambda)*( (a*lambda) > abs(par) )*((a*lambda) - abs(par))/(a-1)
  penaltyd[1]<-0 # no penalty on the intercept
  penaltyd
}


Uee<-function(par,deltaB,x.AB,y.AB,sw){
  ## joint score equation for alpha and beta
  p<-dim(x.AB)[2]
  nAB<-length(deltaB)
  alpha<-par[1:(p+1)]
  beta <-par[(p+2):(2*p+2)]
  x.AB0<-cbind(1,x.AB)
  lpiB<-x.AB0%*%alpha
  piB <-exp(lpiB)/(1+exp(lpiB))
  deltaA<-1-deltaB
  y.AB[which(is.na(y.AB))]<-0
  resB<-(y.AB-x.AB0%*%beta)
  piB<-as.vector(piB)
  sw<-as.vector(sw)
  resB<-as.vector(resB)
  c(apply(x.AB0*deltaB/piB-x.AB0*deltaA*sw,2,sum),
    apply(x.AB0*deltaB*(1/piB-1)*resB,2,sum))/nAB
}


Uee_a<-function(par,deltaB,x.AB,y.AB,sw){
  ## score equation for alpha
  alpha<-par
  nAB<-length(deltaB)
  x.AB0<-cbind(1,x.AB)
  lpiB<-x.AB0%*%alpha
  piB <-exp(lpiB)/(1+exp(lpiB))
  deltaA<-1-deltaB
  piB<-as.vector(piB)
  sw<-as.vector(sw)
  c(apply(x.AB0*deltaB/piB-x.AB0*deltaA*sw,2,sum))/nAB
}


Aee_a<-function(par,deltaB,x.AB,y.AB,sw){
  ## derivative of score equation for alpha
  alpha<-par
  p<-dim(x.AB)[2]
  nAB<-length(deltaB)
  x.AB0<-cbind(1,x.AB)
  lpiB<-x.AB0%*%alpha
  piB <-exp(lpiB)/(1+exp(lpiB))
  deltaA<-1-deltaB
  piB<-as.vector(piB)
  sw<-as.vector(sw)
  negA11<-matrix(0,(p+1),(p+1))
  for(ii in 1:(p+1) ){
    for(jj in ii:(p+1) ){
      negA11[ii,jj]<-negA11[jj,ii]<- sum( deltaB*(1-piB)/piB*x.AB0[,ii]*x.AB0[,jj])
    }
  }
  negA11/nAB
}


loss_alpha<-function(par,deltaB,x.AB,y.AB,sw){
  ## loss function for alpha using the square distance of X between A and B
  alpha<-par
  nAB<-length(deltaB)
  x.AB0<-cbind(1,x.AB)
  lpiB<-x.AB0%*%alpha
  piB <-exp(lpiB)/(1+exp(lpiB))
  deltaA<-1-deltaB
  piB<-as.vector(piB)
  sw<-as.vector(sw)
  NhatA<-sum(deltaA*sw)
  NhatB<-sum(deltaB/piB)
  sum ( apply( ( x.AB0*deltaB/piB/NhatB - x.AB0*deltaA*sw/NhatA ),2,sum)^2 )
}


Uee_binary<-function(par,deltaB,x.AB,y.AB,sw){
  ## joint score equation for alpha and beta
  p<-dim(x.AB)[2]
  if(is.null(p))p=1
  nAB<-length(deltaB)
  alpha<-par[1:(p+1)]
  beta <-par[(p+2):(2*p+2)]

  x.AB0<-cbind(1,x.AB)
  lpiB<-x.AB0%*%alpha
  piB <-exp(lpiB)/(1+exp(lpiB))
  deltaA<-1-deltaB
  mB<-expoit(x.AB0%*%beta)
  resB<-(y.AB- mB)
  piB<-as.vector(piB)
  sw<-as.vector(sw)
  resB<-as.vector(resB)
  mB<-as.vector(mB)
  c(apply(x.AB0*deltaB/piB*mB*(1-mB) -x.AB0*deltaA*sw*mB*(1-mB),2,sum),
    apply(x.AB0*deltaB*(1/piB-1)*resB,2,sum))/nAB
}

expoit<-function(x){
  exp(x)/(1+exp(x))
}
