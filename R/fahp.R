optim_function_y <- function(f,y, LA, LB, zeromat,poimat, i ) {
  lnlik=matrix(nrow = p, ncol = 2)
  for (j in 1:p) {
    # 计算表达式的第一部分
    lnlik[j,1] <- zeromat[i,j]*(c(LA[j,][1],1,0,LA[j,][-1])%*%c(1,f,y[i])) - log(1 + exp(c(LA[j,][1],1,0,LA[j,][-1])%*%c(1,f,y[i])))
    # 计算表达式的第二部分
    lnlik[j,2] <- ifelse(is.na((zeromat[i,j]* ((poimat[i,j] - 1) * (c(LB[j,][1],0,1,LB[j,][-1])%*%c(1,f,y[i])) - exp(c(LB[j,][1],0,1,LB[j,][-1])%*%c(1,f,y[i])) ))), 0, (zeromat[i,j]* ((poimat[i,j] - 1) * (c(LB[j,][1],0,1,LB[j,][-1])%*%c(1,f,y[i])) - exp(c(LB[j,][1],0,1,LB[j,][-1])%*%c(1,f,y[i])) )))
  }
  ss=sum(colMeans(lnlik))
  # if (any(is.infinite(ss) | is.nan(ss) | is.na(ss))) {
  #   return(1e16)
  # }
  if (is.infinite(ss) && ss < 0) {
    return(1e16)
  } else{
    return(- ss)
  }
}

objfunc_y <- function(FA, mu1,mu2, LA, LB,y, zeromat, poimat){
  n <- nrow(zeromat); p <- ncol(zeromat)
  q=ncol(FA)
  none=rep(1, n)
  pone=rep(1, p)
  eps1 <- 1e-20
  pa_1 <- 1 / (1 + exp(-FA%*%t(LA[,c(-1,-(q+2))])-outer(none,LA[,1])-outer(mu1,pone)-outer(y,LA[,(q+2)])))
  lik_1=sum(-zeromat*log(pa_1+eps1)-(1-zeromat)*log(1-pa_1 + eps1))

  pa_2=exp(FA%*%t(LB[,c(-1,-(q+2))])+outer(mu2,pone)+outer(none,LB[,1])+outer(y,LB[,(q+2)]))
  lik_2=sum(-log(dpois(poimat-1, pa_2)+eps1)[!is.na(-log(dpois(poimat-1, pa_2)+eps1))])
  obj <- lik_1/(n*p)+lik_2/(n*p)
  return(obj)
}

proj <- function(a, C) {
  # 检查并处理非正常值
  abnormal_indices <- which(!is.finite(a))
  if (length(abnormal_indices) > 0) {
    cat("警告: 在位置", abnormal_indices, "发现非正常值，已用0替换\n")
    a[abnormal_indices] <- 0
  }

  norm_a <- sqrt(sum(a^2))
  if (norm_a > C) {
    a <- a * (C / norm_a)
  }
  return(a)
}


#' FAHP
#' @description Factor augmented hurdle Poisson inverse model. This method can conduct supervised dimensionality reduction on high dimensional count data containing a large number of zeros.
#'
#' @param zipdata Observation matrix, commonly characterized by zero inflation.
#' @param y A response vector with dimension equal to the number of rows in zipdata.
#' @param q The number of latent factors. Default as 2.
#' @param maxit Maximum number of iterations within optim function, defaults to 300.
#' @param constraint Constraint constant. Default as 5.
#' @param para The number of cores used for parallel processing. Defaults to 0, indicating that parallel processing is not performed.
#'
#' @return
#'  \item{FA}{Factor score matrix.}
#'  \item{mu1}{The intercept vector for hurdle part}
#'  \item{mu2}{The intercept vector for positive count part}
#'  \item{t}{Number of iterations}
#'  \item{obj}{The loglikelihhod results (omitting the constant), which can be used for model selection}
#' @export
#' @author Zhijing Wang
#'
#' @examples
#' n=50
#' p=20
#' x=matrix(rpois(n*p,lambda = 1),ncol = p)
#' y=rnorm(n)
#' fahp(x,y,q=1,maxit = 30,para = 10)

fahp=function(zipdata,y, q=2, maxit=300, constraint=5, para=0){
  if(q==1){
    ss=GFM::gfm(list(zipdata),types = 'poisson', q=q, algorithm = "AM", verbose = FALSE)
  } else{
    ss=GFM::gfm(list(zipdata),types = 'poisson', q=q, algorithm = "VEM", verbose = FALSE)
  }
  y=y-mean(y)
  FA=ss$hH
  n=nrow(zipdata)
  p=ncol(zipdata)
  t=0
  c=1
  dc=Inf
  zeromat <- ifelse(zipdata != 0, 1, 0)
  poimat=ifelse(zipdata == 0, NA, zipdata)
  mu1=log(rowSums(zeromat)/(p-rowSums(zeromat))+1e-8)
  if(any(!is.finite(mu1))){
    warning("mu1 contains non-finite values even after correction. Check input data or logic.")
    mu1[!is.finite(mu1)] <- 0
  }
  mu2=log(rowSums(poimat, na.rm = TRUE)-rowSums(zeromat)+1e-8)
  FAMU=cbind(mu1,mu2,FA)
  LA=matrix(0,nrow = p,ncol = (q+2))
  LB=matrix(0,nrow = p,ncol = (q+2))
  while ((t < maxit)&&(dc>10e-3)) {
    t=t+1
    LA_old=LA
    LB_old=LB
    for (j in 1:p) {
      LA[j,]= tryCatch({
        proj(glm(zeromat[,j]~cbind(FA,y)+offset(mu1),family=binomial)$coefficients, constraint)
      }, error = function(e) {
        message(sprintf("Error in optimization for LA index %d: %s", j, e$message))
        LA_old[j,]+0.01
      })
    }
    for (j in 1:p) {
      LB[j,]=tryCatch({
        proj(glm((poimat[,j]-1)~cbind(FA,y)+offset(mu2), family=poisson)$coefficients, constraint)
      }, error = function(e) {
        message(sprintf("Error in optimization for LB index %d: %s", j, e$message))
        LB_old[j,]+0.01
      })
    }
    FAMU_old=FAMU
    # 判断 para 是否为 TRUE,是否进行并行运算
    if (para>0) {
      cl <- makeCluster(para)
      registerDoParallel(cl)
      FAMU <-foreach(i = 1:n, .combine = rbind, .packages = c("stats"), .export = c("proj", "optim_function_y", "objfunc_y","p"))  %dopar% {
        tryCatch({
          proj(optim(FAMU_old[i,], optim_function_y, y = y, LA = LA, LB = LB, zeromat = zeromat, poimat = poimat, i = i)$par, constraint)
        }, error = function(e) {
          message(sprintf("Error in optimization for FAMU index %d: %s", i, e$message))
          runif((q+2), min = -2, max = 2)  # 返回一个在-2到2之间的随机数
        })
      }
      stopCluster(cl)

    } else {
      for (i in 1:n) {
        FAMU[i,] <- tryCatch({
          proj(optim(FAMU_old[i,], optim_function_y, y = y, LA = LA, LB = LB, zeromat = zeromat, poimat = poimat, i = i)$par, constraint)
        }, error = function(e) {
          message(sprintf("Error in optimization for FAMU index %d: %s", i, e$message))
          runif((q+2), min = -2, max = 2)
        })
      }
    }
    FA=as.matrix(FAMU[,-(1:2)])
    mu1=FAMU[,1]
    mu2=FAMU[,2]
    c_new=objfunc_y(FA, mu1,mu2, LA, LB,y, zeromat, poimat)
    dc=abs((c-c_new)/c)
    c=c_new
  }
  return(list(FA=FA,mu1=mu1,mu2=mu2, LA=LA, LB=LB, t, obj=c))
}

