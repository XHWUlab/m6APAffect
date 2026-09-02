#' @importFrom stats dnbinom predict pnorm
NULL


## codes were took from the QNB packages due to the failure of package installation

.getBaseMean <- function(meth1,meth2,unmeth1,unmeth2,s_t1,s_t2,s_c1,s_c2){
  
  # estimate probability of methylation under a condition
  p1 <- .estimateP(meth1, unmeth1, s_t1, s_c1)
  p2 <- .estimateP(meth2, unmeth2, s_t2, s_c2)
  p0 <- .estimateP(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2))
  
  # estimate the abundance of feature
  q0 <- .estimateQ(cbind(meth1,meth2), cbind(unmeth1,unmeth2), c(s_t1,s_t2), c(s_c1,s_c2),p0)
  q1 <- .estimateQ(meth1,unmeth1,s_t1,s_c1,p0)
  q2 <- .estimateQ(meth2,unmeth2,s_t2,s_c2,p0)
  
  # estimate size e
  e1 <- .estimateE(meth1,unmeth1,s_t1,s_c1,q0)
  e2 <- .estimateE(meth2,unmeth2,s_t2,s_c2,q0)
  
  res <- list(p0=p0,p1=p1,p2=p2,q0=q0,q1=q1,q2=q2,e1=e1,e2=e2)
  return(res)
}

.baseFitPer <- function(meth1,meth2,unmeth1,unmeth2,p1,p2,q,e1,e2,s_t1,s_t2,s_c1,s_c2){
  w_t1 <-.calculateW(meth1,s_t1,e1)
  w_t2 <-.calculateW(meth2,s_t2,e2)
  w_c1 <-.calculateW(unmeth1,s_c1,e1)
  w_c2 <-.calculateW(unmeth2,s_c2,e2)
  
  # locfit 
  fit_t1 <- .locfitW(p1,q,w_t1)
  fit_t2 <- .locfitW(p2,q,w_t2)
  fit_c1 <- .locfitW(p1,q,w_c1)
  fit_c2 <- .locfitW(p2,q,w_c2)
  
  res <- list(fit_t1,fit_t2,fit_c1,fit_c2)
  return(res)
}

.baseFitCondition <- function(meth,unmeth,p,q,e,s_t,s_c){
  w_t <-.calculateW(meth,s_t,e)
  w_c <-.calculateW(unmeth,s_c,e)
  
  # locfit 
  fit_t <- .locfitW(p,q,w_t)
  fit_c <- .locfitW(p,q,w_c)
  
  res <- list(fit_t,fit_c)
  return(res)
}

.getPoolMean <- function(meth,unmeth,s_t,s_c)
{
  # estimate probability of methylation under a condition
  p <- .estimateP(meth,unmeth, s_t,s_c)
  
  # estimate the abundance of feature
  q <- .estimateQ(meth,unmeth,s_t,s_c,p,useAll=TRUE)
  
  # estimate size e
  e <- .estimateE(meth,unmeth,s_t,s_c,q)
  
  res <- list(p=p,q=q,e=e)
  return(res)
}

.poolFit <- function(meth,unmeth,p,q,e,s_t,s_c){
  # calculate methylation reads count variance on common scale
  w_t <-.calculateW(meth,s_t,e)
  w_c <-.calculateW(unmeth,s_c,e)
  
  # locfit 
  fit_t <- .locfitW(p,q,w_t)
  fit_c <- .locfitW(p,q,w_c)
  
  res <- list(fit_t,fit_c)
  return(res)
  
}

.sizeFactor2 <- function(n) {
  temp <- log(colSums(n))
  temp <- temp - mean(temp)
  s_size <- exp(temp)
  return(s_size)
}

.sizeFactor <- function(n, useTotal=FALSE) {
  if (useTotal) {
    temp <- log(colSums(n))
    temp <- temp - mean(temp)
    s_size <- exp(temp)
  } else {
    n <- pmax(n,1e-5)
    log_n <- log(n)
    pseudo <- rowMeans(log_n)
    ratio <- log_n-pseudo
    s_size <- exp(apply(ratio,2,median)) }
  return(s_size)
}

.calculateW <- function(meth,size_t,e_t){
  temp_t <- t(t(meth)/size_t)
  q <- temp_t/e_t
  #  q[is.na(q)] <- 0
  resi <- q-rowMeans(q)
  w <- rowSums(resi^2)/(length(size_t)-1)
  w <- pmax(w,1e-8)
  return(w)
}

.calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)
  
  norow <- length(q)
  nocol <- length(size)
  temp2 <- matrix(1,nrow=norow, ncol=nocol )
  temp3 <- t(t(temp2)/size)/e
  
  z <- rowSums(temp3)*temp
  
  #  z[is.infinite(z)] <- 0
  #  z[is.na(z)] <- 0
  return(z)
}

.estimateE <- function(meth,unmeth,size_t,size_c,q){
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  e <- rowSums(temp_n)/length(size_t)/q
}

.estimateP <- function(meth, unmeth, size_t, size_c) {
  temp_t <- t(t(meth)/size_t)
  temp_c <- t(t(unmeth)/size_c)
  temp_n <- temp_t+temp_c
  p <- rowSums(temp_t)/rowSums(temp_n)
  p[is.na(p)] <- 0.5
  p[is.infinite(p)] <- 0.5
  return(p)
  # which(is.na(p))
  # which(is.infinite(p))  
}

.estimateQ <- function(meth,unmeth,size_t,size_c,p,useAll=TRUE){
  if (useAll) {
    temp_t <- t(t(meth)/size_t)
    temp_c <- t(t(unmeth)/size_c)
    temp_n <- temp_t+temp_c
    q <- rowSums(temp_n)/length(size_t)
  } else {
    temp_c <- t(t(unmeth)/size_c)
    qc <- rowMeans(temp_c)
    q <- qc/(1-p)
  }
  
  q[is.na(q)] <- 0
  return(q)
  # which(is.na(q))
  # which(is.infinite(q))  
}


.locfitW <- function(p,q,w) {
  l <- log(q+1)
  c=c(rep(1,length(p)))
  data=data.frame(cbind(p,l,c,w))
  ID <- which(rowSums(is.na(data))>0)
  #  data <- data[-ID,]
  fit=locfit::locfit(w~locfit::lp(p,l,c),data=data,family="gamma")
  return(fit)
}


.fittedW <- function(p,q,fit){ 
  l <- log(q+1)
  c <- c(rep(1,length(p)))
  data=data.frame(cbind(p,l,c))
  w_fit <- predict(fit,data)
  return(w_fit)
}

.quadNBtest <- function(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c){
  nrows <- length(t)
  pval4 <- rep(1,nrows)
  
  t2 <- t-t1
  
  for (irow in 1:nrows) {
    
    trip <- t[irow]
    
    if (trip<1) {p <- NA} else {
      
      trip_t1 <- 0:t[irow]
      trip_t2 <- t[irow] - trip_t1
      trip_c1 <- n1[irow] - trip_t1
      trip_c2 <- n2[irow] - trip_t2
      
      
      p1 <- dnbinom(x=trip_t1, size=size1_t[irow], mu=mu1_t[irow], log = TRUE)
      p2 <- dnbinom(x=trip_t2, size=size2_t[irow], mu=mu2_t[irow], log = TRUE)
      p3 <- dnbinom(x=trip_c1, size=size1_c[irow], mu=mu1_c[irow], log = TRUE)
      p4 <- dnbinom(x=trip_c2, size=size2_c[irow], mu=mu2_c[irow], log = TRUE)
      options(warn=-1)
      #p4
      p <- p1+p2+p3+p4
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      p <- min(1,2*min(sum(p[1:(t1[irow]+1)]),sum(p[(t1[irow]+1):(t[irow]+1)])))- p[(t1[irow]+1)]/2  
      
      
      pval4[irow] <- p
    }
    
  }
  
  mu <- (mu1_t+mu2_t+mu1_c+mu2_c)/4
  
  res <- data.frame(pval4)
  
  return(res)
}

QNBtest <-
  function(control_ip,treated_ip,control_input,treated_input,
           size.factor=NA,
           mode="auto") {
    options(warn =-1)
    control_ip <- data.frame(control_ip)
    treated_ip <- data.frame(treated_ip)
    control_input <- data.frame(control_input)
    treated_input <- data.frame(treated_input)
    if( any( ncol(control_ip)!=ncol(control_input), ncol(treated_ip)!=ncol(treated_input)) ){
      stop( "IP sample and control sample must be the same replicates" )
    } 
    print("Estimating dispersion for each RNA methylation site, this will take a while ...")
    if(mode=="per-condition"){
      if(anyNA(size.factor)){
        s <- .sizeFactor2(cbind(control_ip,treated_ip,control_input,treated_input))
        s_t1 <- s[1:length(control_ip[1,])]
        s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
        s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
        s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
      }else{
        s_t1 <- size.factor$control_ip
        s_t2 <- size.factor$treated_ip
        s_c1 <- size.factor$control_input
        s_c2 <- size.factor$treated_ip
      }
      mean <- .getBaseMean(control_ip,treated_ip,control_input,treated_input,s_t1,s_t2,s_c1,s_c2)
      p0 <- mean[[1]]
      p1 <- mean[[2]]
      p2 <- mean[[3]]
      q0 <- mean[[4]]
      q1 <- mean[[5]]
      q2 <- mean[[6]]
      e1 <- mean[[7]]
      e2 <- mean[[8]]
      res <- .baseFitPer(control_ip,treated_ip,control_input,treated_input,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
      fit_t1<-res[[1]]
      fit_t2<-res[[2]]
      fit_c1<-res[[3]]
      fit_c2<-res[[4]]
      
    }else if(mode=="pooled"){
      meth<-rbind(control_ip,treated_ip)
      unmeth<-rbind(control_input,treated_input)
      if(anyNA(size.factor)){
        s <- .sizeFactor2(cbind(meth,unmeth))
        s_t <- s[1:length(meth[1,])]
        s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
        s_t1 <- s_t
        s_t2 <- s_t
        s_c1 <- s_c
        s_c2 <- s_c
      }else{
        s_t1 <- size.factor$control_ip
        s_t2 <- size.factor$treated_ip
        s_c1 <- size.factor$control_input
        s_c2 <- size.factor$treated_ip
        s_t <- c(s_t1,s_t2)
        s_c <- c(s_c1,s_c2)
      }
      mean <- .getPoolMean(meth,unmeth,s_t,s_c)
      p1 <- mean[[1]]
      p2 <- p1
      p0 <- p1
      q0 <- mean[[2]]
      q1 <- q0
      q2 <- q0
      e1 <- mean[[3]]
      e2 <- e1
      res <- .baseFitCondition(meth,unmeth,p0,q0,e1,s_t,s_c)
      fit_t1<-res[[1]]
      fit_t2<-res[[1]]
      fit_c1<-res[[2]]
      fit_c2<-res[[2]]
      
    }else if(mode=="blind"){
      meth<-cbind(control_ip,treated_ip)
      unmeth<-cbind(control_input,treated_input)
      if(anyNA(size.factor)){
        s <- .sizeFactor2(cbind(meth,unmeth))
        s_t <- s[1:length(meth[1,])]
        s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
        s_t1 <- s[1:length(control_ip[1,])]
        s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
        s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
        s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
      }else{
        s_t1 <- size.factor$control_ip
        s_t2 <- size.factor$treated_ip
        s_c1 <- size.factor$control_input
        s_c2 <- size.factor$treated_ip
        s_t <- c(s_t1,s_t2)
        s_c <- c(s_c1,s_c2)
      }
      mean <- .getPoolMean(meth,unmeth,s_t,s_c)
      p1 <- mean[[1]]
      p2 <- p1
      p0 <- p1
      q0 <- mean[[2]]
      q1 <- q0
      q2 <- q0
      e1 <- mean[[3]]
      e2 <- e1
      res <- .poolFit(meth,unmeth,p0,q0,e1,s_t,s_c)
      fit_t1<-res[[1]]
      fit_t2<-res[[1]]
      fit_c1<-res[[2]]
      fit_c2<-res[[2]]
    }else if(mode=="auto"){
      rep1 <- ncol(control_ip)
      rep2 <- ncol(treated_ip)
      if((rep1==rep2&&rep1>1)||(rep1!=rep2&&min(rep1,rep2)>1)){
        if(anyNA(size.factor)){
          s <- .sizeFactor2(cbind(control_ip,treated_ip,control_input,treated_input))
          s_t1 <- s[1:length(control_ip[1,])]
          s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
          s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
          s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
        }else{
          s_t1 <- size.factor$control_ip
          s_t2 <- size.factor$treated_ip
          s_c1 <- size.factor$control_input
          s_c2 <- size.factor$treated_ip
        }
        mean <- .getBaseMean(control_ip,treated_ip,control_input,treated_input,s_t1,s_t2,s_c1,s_c2)
        p0 <- mean[[1]]
        p1 <- mean[[2]]
        p2 <- mean[[3]]
        q0 <- mean[[4]]
        q1 <- mean[[5]]
        q2 <- mean[[6]]
        e1 <- mean[[7]]
        e2 <- mean[[8]]
        res <- res <- .baseFitPer(control_ip,treated_ip,control_input,treated_input,p1,p2,q0,e1,e2,s_t1,s_t2,s_c1,s_c2)
        fit_t1<-res[[1]]
        fit_t2<-res[[2]]
        fit_c1<-res[[3]]
        fit_c2<-res[[4]]
        
      }else if((rep1==rep2&&rep1==1)||(rep1!=rep2&&min(rep1,rep2)<2)){
        meth=cbind(control_ip,treated_ip)
        unmeth=cbind(control_input,treated_input)
        
        if(anyNA(size.factor)){
          s <- .sizeFactor2(cbind(meth,unmeth))
          s_t <- s[1:length(meth[1,])]
          s_c <- s[(length(meth[1,])+1):(length(cbind(meth,unmeth)[1,]))]
          s_t1 <- s[1:length(control_ip[1,])]
          s_t2 <- s[(length(control_ip[1,])+1):(length(cbind(control_ip,treated_ip)[1,]))]
          s_c1 <- s[(length(cbind(control_ip,treated_ip)[1,])+1):(length(cbind(control_ip,treated_ip,control_input)[1,]))]
          s_c2 <- s[(length(cbind(control_ip,treated_ip,control_input)[1,])+1):(length(cbind(control_ip,treated_ip,control_input,treated_input)[1,]))]
        }else{
          s_t1 <- size.factor$control_ip
          s_t2 <- size.factor$treated_ip
          s_c1 <- size.factor$control_input
          s_c2 <- size.factor$treated_ip
          s_t <- c(s_t1,s_t2)
          s_c <- c(s_c1,s_c2)
        }
        mean <- .getPoolMean(meth,unmeth,s_t,s_c)
        p1 <- mean[[1]]
        p2 <- p1
        p0 <- p1
        q0 <- mean[[2]]
        q1 <- q0
        q2 <- q0
        e1 <- mean[[3]]
        e2 <- e1
        
        res <- .poolFit(meth,unmeth,p0,q0,e1,s_t,s_c)
        fit_t1<-res[[1]]
        fit_t2<-res[[1]]
        fit_c1<-res[[2]]
        fit_c2<-res[[2]]
        
      }
    }
    
    
    # if (is.na(output.dir)) {
    #   output.dir <- getwd()
    # }
    
    # path <- paste(output.dir,"dispersion.pdf",sep = '/')
    # if(plot.dispersion){
    #   .plotDispersion(fit_t1,fit_c1,fit_t2,fit_c2,path)
    # }
    
    
    # calculate z
    z_t1 <- .calculateZ(q0,p1,s_t1,e1)
    z_t2 <- .calculateZ(q0,p2,s_t2,e2)
    z_c1 <- .calculateZ(q0,(1-p1),s_c1,e1)
    z_c2 <- .calculateZ(q0,(1-p2),s_c2,e2)
    
    # get estimate w
    w_fit_t1 <- .fittedW(p0,q0,fit_t1)
    w_fit_t2 <- .fittedW(p0,q0,fit_t2)
    w_fit_c1 <- .fittedW(p0,q0,fit_c1)
    w_fit_c2 <- .fittedW(p0,q0,fit_c1)
    
    # get estimate of upi
    ups_t1 <- pmax(w_fit_t1 - z_t1, 1e-8)
    ups_t2 <- pmax(w_fit_t2 - z_t2, 1e-8)
    ups_c1 <- pmax(w_fit_c1 - z_c1, 1e-8)
    ups_c2 <- pmax(w_fit_c2 - z_c2, 1e-8)
    
    # get all means
    mu_t1 <- (e1*q0*p0)%*%t(as.numeric(s_t1))
    mu_t2 <- (e2*q0*p0)%*%t(as.numeric(s_t2))
    mu_c1 <- (e1*q0*(1-p0))%*%t(as.numeric(s_c1))
    mu_c2 <- (e2*q0*(1-p0))%*%t(as.numeric(s_c2))
    
    # get all variance
    raw_t1 <- (e1%*%t(s_t1))^2*ups_t1
    raw_t2 <- (e2%*%t(s_t2))^2*ups_t2
    raw_c1 <- (e1%*%t(s_c1))^2*ups_c1
    raw_c2 <- (e2%*%t(s_c2))^2*ups_c2
    
    # put mu together
    mu1_t <- rowSums(mu_t1)
    mu2_t <- rowSums(mu_t2)
    mu1_c <- rowSums(mu_c1)
    mu2_c <- rowSums(mu_c2)
    
    # put size together
    size1_t <- (mu1_t^2)/rowSums(raw_t1)
    size2_t <- (mu2_t^2)/rowSums(raw_t2)
    size1_c <- (mu1_c^2)/rowSums(raw_c1)
    size2_c <- (mu2_c^2)/rowSums(raw_c2)
    
    # observation together
    t1 <- rowSums(control_ip)
    t2 <- rowSums(treated_ip)
    c1 <- rowSums(control_input)
    c2 <- rowSums(treated_input)
    t <- t1 + t2
    n1 <- t1 + c1
    n2 <- t2 + c2
    
    raw <- (rowSums(raw_t1)+rowSums(raw_t2)+rowSums(raw_c1)+rowSums(raw_c2))/4
    # go to test
    res <- .quadNBtest(t1,t,n1,n2,mu1_t,mu2_t,mu1_c,mu2_c,size1_t,size2_t,size1_c,size2_c)
    
    # add fc
    p1 <- .estimateP(control_ip,control_input,s_t1,s_c1)
    p2 <- .estimateP(treated_ip,treated_input,s_t2,s_c2)
    log2.RR <- log2(p2/p1)
    
    p.treated <- p2
    p.control <- p1
    
    log2.OR <- log2((rowSums(t(t(treated_ip)/s_t2))/rowSums(t(t(treated_input)/s_c2)))/(rowSums(t(t(control_ip)/s_t1))/rowSums(t(t(control_input)/s_c1))))
    m1 <- rowSums(t(t(control_ip)/s_t1))
    m2 <- rowSums(t(t(treated_ip)/s_t2))
    u1 <- rowSums(t(t(control_input)/s_c1))
    u2 <- rowSums(t(t(treated_input)/s_c2))
    mfc <- log2(m1)-log2(m2)
    ufc <- log2(u1)-log2(u2)
    
    padj <- p.adjust( res[,1], method="BH" )
    res <- data.frame(p.treated,p.control,log2.RR,log2.OR,res[,1],q0,padj)
    colnames(res) <- c("p.treated","p.control","log2.RR","log2.OR","pvalue","q","padj")
    
    return(res)}


