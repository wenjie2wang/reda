#############################
# Author: Haoda Fu
# Update Time: Dec 31, 2012
# Function: fit piece wise gamma frailty model for recurrent events
# Data Format: 
# > dataM.test  <- rbind(c(1,20,1,5),c(1,70,1,5),c(1,100,0,5),c(2,60,1,8),c(2,100,0,8))
# > colnames(dataM.test)  <- c("ID","Time","Event","X");
# > BL.pieces  <- c(50,100);
# > dataM.test
# ID Time Event X
# [1,]  1   20     1 5
# [2,]  1   70     1 5
# [3,]  1  100     0 5
# [4,]  2   60     1 8
# [5,]  2  100     0 8
#############################

#' fit heart model.
#'
#' \code{functionname} returns fitted model results.
#'
#' This is a test Roxygen comments
#'
#' @param ... Numeric, complex, or logical vectors.
#' @param na.rm A logical scalar. 
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org/wiki/Integer_overflow} occurs, the output
#'   will be NA with a warning. Otherwise it will be a length-one numeric or
#'   complex vector.
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F, F, F, T, T)
HEART.PW.gammaFrailty  <- function(dataM,BaselinePieces,ini){
  fit  <- nlm(logL.HEART,ini,dataM=dataM,BaselinePieces=BaselinePieces,hessian=T);
  length.par  <- length(ini);
  length.alpha  <- length(BaselinePieces);
  length.beta  <- length.par - length.alpha -1;  
  est.beta  <- matrix(NA,nrow= length.beta, ncol=3);
  colnames(est.beta)  <- c("beta","se(beta)","two sided p-value");
  
  se.vec  <- sqrt(diag(solve(fit$hessian)));
  
  est.beta[,1]  <- fit$estimate[1:length.beta];
  est.beta[,2]  <- se.vec[1:length.beta];
  est.beta[,3]  <- 2*pnorm(-abs(est.beta[,1]/est.beta[,2]));
  
  est.theta  <- matrix(NA,nrow=1,ncol=2);
  colnames(est.theta)  <- c("theta","se(theta)");  
  est.theta[1,]  <- c(fit$estimate[length.beta+1],se.vec[length.beta+1]);
  
  est.alpha  <- matrix(NA,nrow=length.alpha,ncol=2);
  colnames(est.alpha)  <- c("alpha","se(alpha)");
    
  est.alpha[,1]  <- fit$estimate[(length.beta+2):length.par];
  est.alpha[,2]  <- se.vec[(length.beta+2):length.par];  
  results <- list("beta"=est.beta,"theta"=est.theta,"alpha"=est.alpha,"Convergence"=fit$code, "Fisher Info Matrix"=fit$hessian);
  return(results);  
}

plot.baseline  <- function(est.results,BaselinePieces,CI = T, CI.level = 0.95){  
 n.xx  <- 1000;
 n.pieces  <- length(BaselinePieces);
 BL.segments  <- c(BaselinePieces[1],diff(BaselinePieces));
 xx  <- seq(0,max(BaselinePieces),length=n.xx);
 indx  <- sapply(xx,whereT,BaselinePieces=BaselinePieces);
 LinCom.M  <- t(apply(array(indx),1,function(ind.indx) (c(BL.segments[1:ind.indx],rep(0,n.pieces -ind.indx)))));
 CMF.B4.indx  <- c(0,BaselinePieces)[indx]; 
 LinCom.M[(indx-1)*n.xx+1:n.xx]  <- xx-CMF.B4.indx; 
 n.par  <- dim(est.results$'Fisher Info Matrix')[1];
 Cov.M  <- solve(est.results$'Fisher Info Matrix')[c((n.par-n.pieces+1):n.par),c((n.par-n.pieces+1):n.par)];
 CI.band <- qnorm((1+CI.level)/2)*sqrt(diag(LinCom.M%*%Cov.M%*%t(LinCom.M)));
 baseline.mean  <- LinCom.M%*%est.results$alpha[,1];
 ymax  <- max(baseline.mean+CI.band);
 plot(xx, baseline.mean,type="l",lwd=2,ylim=c(0,ymax),xlim=c(0,max(BaselinePieces)), xlab="Time",ylab="MCF",main="Mean Cumulative Function");
 lines(xx,baseline.mean+CI.band, lty=3);
 lines(xx,baseline.mean-CI.band, lty=3); 
}

whereT  <- function(tt,BaselinePieces){
  min(which(tt<=BaselinePieces));
}

rho_0  <- function(par.BaselinePW, BaselinePieces,Tvec)
{
  indx  <- apply(as.array(Tvec),1,whereT,BaselinePieces);
  return(par.BaselinePW[indx]);
}

mu0  <- function(par.BaselinePW, BaselinePieces,Tvec)
{
  indx  <- apply(as.array(Tvec),1,whereT,BaselinePieces); 
  BL.segments  <- c(BaselinePieces[1],diff(BaselinePieces));
  CumMean.Pieces  <- diffinv(BL.segments*par.BaselinePW)[-1]; #The CMF at each time point  
  return(CumMean.Pieces[indx]-(BaselinePieces[indx]-Tvec)*par.BaselinePW[indx]);
}

dmu0_alpha <- function(tt,BaselinePieces){
  BL.segments  <- c(BaselinePieces[1],diff(BaselinePieces));
  indx  <- min(which(tt<=BaselinePieces));
  value  <- BL.segments;
  n.pieces  <- length(BaselinePieces);
  if(indx == n.pieces){
    value[n.pieces]  <- tt-BaselinePieces[n.pieces-1];
  } else if (indx > 1){
    value[(indx+1):n.pieces]  <- 0;
    value[indx]  <- tt-BaselinePieces[indx-1]
  } else {
    value[(indx+1):n.pieces]  <- 0;
    value[indx]  <- tt;
  }    
  return(value);
}

logL.HEART  <- function(par,dataM,BaselinePieces){
  npieces  <- length(BaselinePieces);  
  if(BaselinePieces[npieces] < max(dataM$Time)){
    BaselinePieces[npieces] <- max(dataM$Time)+0.00000001;
    warning("Extend the Baseline Pieces to adjust the data");
  }
  nbeta  <- dim(dataM)[2]-3;
  par.beta  <- par[1:nbeta];
  par.theta  <- par[nbeta+1];
  par.alpha  <- par[(nbeta+2):length(par)];
  m  <- length(unique(dataM$ID));
  expXBeta  <- exp(as.matrix(dataM[,4:(3+nbeta)])%*%(as.matrix(par.beta)));
  rho_0_ij  <- rho_0(par.BaselinePW=par.alpha, BaselinePieces=BaselinePieces,dataM$Time[dataM$Event==1]);
  rho_i  <- expXBeta[dataM$Event==1]*rho_0_ij;
  rho_i[rho_i <1e-100]  <- 1e-100;
  sum.log.rho_i  <- sum(log(rho_i));
  n_ij  <- table(dataM$ID)[order(unique(dataM$ID))]-1; # these codes to make sure that the order will not change if the patient ID is not ordered
  theta_j_1  <- par.theta+sequence(n_ij)-1; #if there is a subject with 0 event, the sequence will not be generated for this subject
  theta_j_1[theta_j_1<1e-100] <- 1e-100;
  sum.log_theta_j_1  <- sum(log(theta_j_1)); 
  mu0i  <- mu0(par.BaselinePW=par.alpha, BaselinePieces=BaselinePieces,dataM$Time[dataM$Event==0])
  mui  <- mu0i*expXBeta[dataM$Event==0];
  mui_theta <- par.theta+mui;
  mui_theta[mui_theta < 1e-100 ] <- 1e-100;
  sum.log_theta_mui  <- sum((n_ij+par.theta)*log(mui_theta));
  if (par.theta <1e-100){    
    par.theta <- 1e-100
  }
  logLH  <- m*par.theta*log(par.theta)+sum.log.rho_i+sum.log_theta_j_1 - sum.log_theta_mui;
  penal  <- ifelse(par.theta < 0 | min(par.alpha) < 0, 1e+50,0);  
  negLH  <- -logLH+penal;  
  ### Calculate the gradient  
  dl_dbeta  <- apply(diag((n_ij-mui)/(par.theta+mui)*par.theta)%*%as.matrix(dataM[dataM$Event==0,4:(3+nbeta)]),2,sum);
  dl_dtheta  <- m+m*log(par.theta)+sum(1/(par.theta+sequence(n_ij)-1))-sum((n_ij+par.theta)/(par.theta+mui))-sum(log(mui_theta));
  indx  <- apply(as.array(dataM$Time[dataM$Event==1]),1,whereT,BaselinePieces);
  if(length(unique(indx))<length(par.alpha)){
    stop("Some segements have zero events!")   
  }   
  indx_taui <- apply(as.array(dataM$Time[dataM$Event==0]),1,whereT,BaselinePieces);
  dl_dalpha.part2  <- diag((n_ij+par.theta)/(par.theta+mui)*expXBeta[dataM$Event==0])%*%t(apply(array(dataM$Time[dataM$Event==0]),1,dmu0_alpha,BaselinePieces));
  dl_dalpha  <- 1/par.alpha*table(indx)- apply(dl_dalpha.part2,2,sum);
  attr(negLH,"gradient") <- -c(dl_dbeta, dl_dtheta,dl_dalpha);  
  return(negLH)
}

