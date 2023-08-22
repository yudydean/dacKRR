library(mvtnorm)
library(MASS)
library(randtoolbox)
library(gss)
library(listdtr)
library(CVST)
library(groupdata2)
library(DRR)
library(randomForest)
library(dplyr)

## The ssa function takes the pre-defined smoothing parameter as input. 
## Author: Xiaoxiao Sun
## 10/31/2016
##input parameters. All the parameters except lambda and theta are the same as the those in original
##ssanova function. 
##lambda: log10(nlambda), from ssanova function
##theta: log10(theta), from ssanova funciton

#############################################
## ssa
#############################################
ssa <- function(formula,type=NULL,data=list(),lambda,theta,weights,subset,
                offset,na.action=na.omit,method="v",alpha=1.4,varht=1,
                id.basis=NULL,nbasis=NULL,seed=NULL)
{
  ## Obtain model frame and model terms
  mf <- match.call()
  mf$type <- mf$method <- mf$varht <- NULL
  mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
  mf$lambda <- mf$theta <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf,parent.frame())
  wt <- model.weights(mf)
  ## Generate sub-basis
  nobs <- dim(mf)[1]  
  if (is.null(id.basis)) {
    if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
    if (nbasis>=nobs)  nbasis <- nobs
    if (!is.null(seed))  set.seed(seed)
    id.basis <- sample(nobs,nbasis,prob=wt)
  }
  else {
    if (max(id.basis)>nobs|min(id.basis)<1)
      stop("gss error in ssanova: id.basis out of range")
    nbasis <- length(id.basis)
  }
  ## Generate terms
  term <- mkterm(mf,type)
  ## Generate s, r, and y
  s <- r <- NULL
  nq <- 0
  for (label in term$labels) {
    if (label=="1") {
      s <- cbind(s,rep(1,len=nobs))
      next
    }
    x <- mf[,term[[label]]$vlist]
    x.basis <- mf[id.basis,term[[label]]$vlist]
    nphi <- term[[label]]$nphi
    nrk <- term[[label]]$nrk
    if (nphi) {
      phi <- term[[label]]$phi
      for (i in 1:nphi)
        s <- cbind(s,phi$fun(x,nu=i,env=phi$env))
    }
    if (nrk) {
      rk <- term[[label]]$rk
      for (i in 1:nrk) {
        nq <- nq+1
        r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq))
      }
    }
  }  
  if (is.null(r))
    stop("gss error in ssanova: use lm for models with only unpenalized terms")
  if (qr(s)$rank<dim(s)[2])
    stop("gss error in ssanova: unpenalized terms are linearly dependent")
  ## Prepare the data
  y <- model.response(mf,"numeric")
  offset <- model.offset(mf)
  if (!is.null(offset)) {
    term$labels <- c(term$labels,"offset")
    term$offset <- list(nphi=0,nrk=0)
    y <- y - offset
  }
  if (!is.null(wt)) wt <- sqrt(wt)
  ## Fit the model
  if (nq==1) {
    r <- r[,,1]
    z <- sreg(s,r,r[id.basis,],y,lambda,theta,wt,method,alpha,varht)
  }
  else z <- mreg(s,r,id.basis,y,lambda,theta,wt,method,alpha,varht)
  ## Brief description of model terms
  desc <- NULL
  for (label in term$labels)
    desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
  desc <- rbind(desc,apply(desc,2,sum))
  rownames(desc) <- c(term$labels,"total")
  colnames(desc) <- c("Unpenalized","Penalized")
  ## Return the results
  obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc,alpha=alpha,
                id.basis=id.basis),z)
  class(obj) <- c("ssanova")
  obj
}

## Fit Single Smoothing Parameter (Gaussian) REGression
sreg <- function(s,r,q,y,lambda,theta,wt,method,alpha,varht)
{
  alpha <- abs(alpha)
  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  nxiz <- nxi 
  nn <- nxiz + nnull
  if (!is.null(wt)) {
    y <- wt*y
    s <- wt*s
    r <- wt*r
  }
  
  ## weighted q using lambda and theta   
  q.wk <- 10^(lambda+theta)*q
  
  z <- .Fortran("reg",
                as.double(cbind(s,10^theta*r)), as.integer(nobs), as.integer(nnull),
                as.double(q.wk), as.integer(nxiz), as.double(y),
                as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                as.double(alpha), varht=as.double(varht),
                score=double(1), dc=double(nn),
                as.double(.Machine$double.eps),
                chol=double(nn*nn), double(nn),
                jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
                PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
  if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
  assign("fit",z[c(1:5,7)],inherits=TRUE)
  
  se.q.wk <- 10^theta*q
  se.aux <- seaux(s,10^theta*r,se.q.wk,lambda,fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  c(list(method=method,theta=theta,c=c,d=d,nlambda=lambda),
    fit[-3],list(se.aux=se.aux))
}

## Fit Multiple Smoothing Parameter (Gaussian) REGression
mreg <- function(s,r,id.basis,y,lambda,theta,wt,method,alpha,varht)
{
  
  alpha <- abs(alpha)
  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  nxiz <- nxi
  nn <- nxiz + nnull
  nq <- dim(r)[3]
  ## weighted q using lambda and theta   
  r.wk0 <- 0
  for (i in 1:nq) {
    r.wk0 <- r.wk0 + 10^theta[i]*r[,,i]
  }      
  qq.wk <- r.wk0[id.basis,]
  q.wk <- 10^lambda*qq.wk
  
  if (!is.null(wt)) {
    y.wk <- wt*y
    s.wk <- wt*s
    r.wk0 <- wt*r.wk0
  } else {
    y.wk <- y
    s.wk <- s
  }
  
  z <- .Fortran("reg",
                as.double(cbind(s.wk,r.wk0)), as.integer(nobs), as.integer(nnull),
                as.double(q.wk), as.integer(nxiz), as.double(y.wk),
                as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                as.double(alpha), varht=as.double(varht),
                score=double(1), dc=double(nn),
                as.double(.Machine$double.eps),
                chol=double(nn*nn), double(nn),
                jpvt=as.integer(c(rep(1,nnull),rep(0,nxiz))),
                wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
                PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
  if (z$info) stop("gss error in ssanova: evaluation of GML score fails")
  assign("fit",z[c(1:5,7)],inherits=TRUE)
  
  r.wk <- 0
  for (i in 1:nq) {
    r.wk <- r.wk + 10^theta[i]*r[,,i]
  }
  se.qq.wk <- r.wk[id.basis,]
  se.q.wk <- se.qq.wk
  
  if (!is.null(wt)) {
    s <- wt*s
    r.wk <- wt*r.wk
  }
  
  se.aux <- seaux(s,r.wk,se.q.wk,lambda,fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  c(list(method=method,theta=theta[1:nq],c=c,d=d,nlambda=lambda),fit[-3],list(se.aux=se.aux))
}

## Auxiliary Quantities for Standard Error Calculation
seaux <- function(s,r,q,nlambda,fit)
{
  nnull <- dim(s)[2]
  nn <- nnull +  dim(q)[1]
  zzz <- eigen(q,symmetric=TRUE)
  rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
  val <- zzz$val[1:rkq]
  vec <- zzz$vec[,1:rkq,drop=FALSE]
  if (nnull) {
    wk1 <- qr(s)
    wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),]
  }
  else wk1 <- r%*%vec
  wk2 <- t(t(wk1)/sqrt(val))
  wk2 <- t(wk2)%*%wk2
  wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)
  wk2 <- (wk2+t(wk2))/2
  wk2 <- t(wk2/sqrt(val))/sqrt(val)
  wk2 <- diag(1/val,dim(wk2)[1])-wk2
  z <- .Fortran("regaux",
                as.double(fit$chol), as.integer(nn),
                as.integer(fit$jpvt), as.integer(fit$rkv),
                drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                PACKAGE="gss")[c("drcr","sms")]
  drcr <- matrix(z$drcr,nn,rkq)
  dr <- drcr[1:nnull,,drop=FALSE]
  sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
  wk1 <- matrix(0,nnull+rkq,nnull+rkq)
  wk1[1:nnull,1:nnull] <- sms
  wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
  wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
  z <- chol(wk1,pivot=TRUE)
  wk1 <- z
  rkw <- attr(z,"rank")
  while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
  wk1[row(wk1)>col(wk1)] <- 0
  if (rkw<nnull+rkq)
    wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
  hfac <- wk1
  hfac[,attr(z,"pivot")] <- wk1
  list(vec=vec,hfac=hfac)
}

#############################################
## esps-m
#############################################

esps = function(form, data, sam.size=NULL, r=NULL, iter.c=10, iter.p=5){
  ##Author: Xiaoxiao Sun
  mf = model.frame(form, data = data)
  sam.theta <- NULL
  sam.lamb <- NULL
  obs = dim(data)[1]
  for(t in 1:iter.c){
    set.seed(t)
    sam.indx <- sample(obs, sam.size)
    sam.dat <- data[sam.indx, ]
    sam.fit <- ssanova(form, data=sam.dat, seed=t, alpha=1.0)
    sam.theta <- rbind(sam.theta, sam.fit$theta)
    sam.lamb[t] <- sam.fit$nlambda
  }
  me.lamb <- median(sam.lamb)
  me.theta <- apply(sam.theta, 2, median)
  me.lamb.org = 10^me.lamb/sam.size
  if(is.null(r)){
    cat("r=4 for the univariate case and r=3 for the multivariate case!", "\n")
  }
  p.c=1
  list(lambda=me.lamb.org, theta=me.theta, p=p.c)
}



#############################################
#### Adaptive Distributed Smoothing Spline
# divide
SlDivide_repeat<-function(my.data, nnode, basecut = 1, copy = TRUE,
                          breaks = 10, slicebound = NULL, downsample = TRUE, downsamplesize = 20){
  ##### basecut = c: when obs in slide >= max(slidesize)/c, partition; otherwise, partition + sample 
  ##### my.data: dataframe, x = x, y = y
  SampleSlDV = NULL
  my.y = my.data$y
  N = length(my.y)
  if(is.null(slicebound)){
    my.y.slice <- hist(my.y,breaks=breaks)$breaks
  } else my.y.slice <- c(min(my.y), slicebound, max(my.y))
  
  nslice <- length(my.y.slice)-1
  
  ## number of samples in each slides
  slideSize <- c()
  for(ii in 1:(nslice)){  
    my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
    slideSize <- c(slideSize, length(my.which))
  }
  
  basen = max(slideSize)/basecut
  
  data.list = SamplePerNode = vector("list", nnode)
  #my.sample <- integer()
  for(ii in 1:(nslice)){ 
    my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
    if(length(my.which) > 0){
      if(slideSize[ii] >= max(basen, nnode)){
        df.my.which <- data.frame(lab = my.which)
        ### partition
        n_per_node = floor(slideSize[ii]/nnode)
        pt = partition(df.my.which, p = rep(n_per_node,nnode), list_out = F)
        for(i in 1:nnode){
          pti = pt$lab[pt$.partitions == i]  ## Ai
          if(downsample){
            ptii = sample(pti, min(length(pti),downsamplesize))
            pti = ptii
          }
          SamplePerNode[[i]] = c(SamplePerNode[[i]], pti)
        } 
      } else{
        if(copy){
          ncopy = floor(basen/slideSize[ii])
          if(slideSize[ii] < nnode){
            ncopy = floor(nnode/slideSize[ii])+1
          }
          my.which.rep <- rep(my.which, each = ncopy)
          my.which.rep <- sample(my.which.rep)
        } else{
          my.which.rep <- sample(my.which, floor(basen), replace = T)
        }
        df.my.which.rep <- data.frame(lab = my.which.rep)
        ### partition
        n_per_node = floor(length(my.which.rep)/nnode)
        pt = partition(df.my.which.rep, p = rep(n_per_node,nnode), list_out = F)
        for(i in 1:nnode){
          pti = pt$lab[pt$.partitions == i]  ## Ai
          if(downsample){
            ptii = sample(pti, min(length(pti),downsamplesize))
            pti = ptii
          }
          pti <- pti[!duplicated(pti)]
          SamplePerNode[[i]] = c(SamplePerNode[[i]], pti)
        } 
      }
    }
  }
  
  for(i in 1:nnode){
    data.list[[i]] = my.data[SamplePerNode[[i]], ]
    #SampleSlDV = rbind(SampleSlDV, SamplePerNode[[i]])
  }
  # return(list(data.list = data.list, samples = SampleSlDV, SamplePerSlice = slideSize))
  return(list(data.list = data.list, SamplePerSlice = slideSize))
}



#### Adaptive Distributed Smoothing Spline
# divide
SlDivide<-function(my.data, nnode, basecut = 1, 
                   breaks = 10, slicebound = NULL, downsample = TRUE, downsamplesize = 20){
  ##### basecut = c: when obs in slide >= max(slidesize)/c, partition; otherwise, partition + sample 
  ##### my.data: dataframe, x = x, y = y
  SampleSlDV = NULL
  my.y = my.data$y
  N = length(my.y)
  if(is.null(slicebound)){
    my.y.slice <- hist(my.y,breaks=breaks)$breaks
  } else my.y.slice <- c(min(my.y), slicebound, max(my.y))
  
  nslice <- length(my.y.slice)-1
  
  ## number of samples in each slides
  slideSize <- c()
  for(ii in 1:(nslice)){  
    my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
    slideSize <- c(slideSize, length(my.which))
  }
  
  order.lab <- order(-slideSize)
  nbase = floor(max(slideSize)/nnode/basecut)
  
  M1 = sum(slideSize - max(slideSize)/basecut >=0)
  
  M21 = sum(slideSize - nbase >= 0)
  M22 = sum(slideSize >= nnode)
  M2 = min(M21, M22)
  
  data.list = SamplePerNode = vector("list", nnode)
  #my.sample <- integer()
  for(ii in order.lab[1:M1]){  
    my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
    #my.sample <- union(my.sample, sample(my.which, min(nobs.slice, length(my.which)) ))
    df.my.which <- data.frame(lab = my.which)
    ### partition
    n_per_node = floor(slideSize[ii]/nnode)
    pt = partition(df.my.which, p = rep(n_per_node,nnode), list_out = F)
    for(i in 1:nnode){
      pti = pt$lab[pt$.partitions == i]  ## Ai
      if(downsample){
        ptii = sample(pti, min(length(pti),downsamplesize))
        pti = ptii
      }
      SamplePerNode[[i]] = c(SamplePerNode[[i]], pti)
    } 
  }
  
  if(M1<M2){
    for(ii in order.lab[(M1+1):M2]){  
      my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
      #my.sample <- union(my.sample, sample(my.which, min(nobs.slice, length(my.which)) ))
      df.my.which <- data.frame(lab = my.which)
      ### partition
      n_per_node = floor(slideSize[ii]/nnode)
      pt = partition(df.my.which, p = rep(n_per_node,nnode), list_out = F)
      for(i in 1:nnode){
        pti = pt$lab[pt$.partitions == i]  ## Ai
        if(downsample){
          ptii = sample(pti, min(length(pti),downsamplesize))
          pti = ptii
        }
        labi = sample(setdiff(my.which, pti), nbase - n_per_node, replace = F) ## Bi
        SamplePerNode[[i]] = c(SamplePerNode[[i]], pti, labi)
      } 
    }
  }
    
  if(M2<nslice){
      for(ii in order.lab[(M2+1):nslice]){
        my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
        for(i in 1:nnode) SamplePerNode[[i]] = union(SamplePerNode[[i]], my.which)
      }
  }
  
  for(i in 1:nnode){
      data.list[[i]] = my.data[SamplePerNode[[i]], ]
      #SampleSlDV = rbind(SampleSlDV, SamplePerNode[[i]])
  }
  # return(list(data.list = data.list, samples = SampleSlDV, SamplePerSlice = slideSize))
  return(list(data.list = data.list, SamplePerSlice = slideSize))
}

ADivide<-function(my.data, nnode, breaks = 10, slicebound = NULL, downsample = TRUE, downsamplesize = 20){
  ##### my.data: dataframe, x = x, y = y
  SampleADV = NULL
  my.y = my.data$y
  if(is.null(slicebound)){
    my.y.slice <- hist(my.y,breaks=breaks)$breaks
  } else my.y.slice <- c(min(my.y), slicebound, max(my.y))
  
  nslice <- length(my.y.slice)-1
  nob = floor(length(my.y)/nnode)
  if(downsample) nob = downsamplesize*nslice
  nobs.slice <- floor(nob/nslice)
  
  data.list = vector("list", nnode)
  for(i in 1:nnode){
    my.sample <- integer()
    SamplePerSlice = c()
    for(ii in 1:(nslice)){  
      my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
      my.sample <- union(my.sample, sample(my.which, min(nobs.slice, length(my.which)) ))
      SamplePerSlice = c(SamplePerSlice, length(my.which))
    }
    
    labi = c(my.sample, sample(setdiff(1:length(my.y),my.sample), nob-length(my.sample), replace = F))
    data.list[[i]] = my.data[labi, ]
    SampleADV = rbind(SampleADV, labi)
  }
  return(list(data.list = data.list, samples = SampleADV, SamplePerSlice = SamplePerSlice))
}

UDivide<-function(my.data, nnode, downsample = TRUE, downsamplesize = 20){
  ##### my.data: dataframe, x = x, y = y
  SampleUDV = NULL
  my.y = my.data$y
  nob = floor(length(my.y)/nnode)
  if(downsample) nob = downsamplesize*2
  data.list = vector("list", nnode)
  for(i in 1:nnode){
    labi = sample.int(length(my.y),size=nob,replace=F)
    data.list[[i]] = my.data[labi, ]
    SampleUDV = rbind(SampleUDV, labi)
  }
  return(list(data.list = data.list, samples = SampleUDV))
}

SDivide<-function(my.data, nnode, setseed = 2020){
  ##### my.data: dataframe, x = x, y = y
  SampleSDV = NULL
  my.y = my.data$y
  X = my.data$x
  nob = floor(length(my.y)/nnode)
  data.list = vector("list", nnode)
  for(i in 1:nnode){
    design1<- sobol(n=nob ,dim = 1, init = T, scrambling = 1, seed=setseed+10*i)
    labi <- as.vector(nabor::knn(X, design1, k = 1)$nn.idx)
    data.list[[i]] = my.data[labi, ]
    SampleSDV = rbind(SampleSDV, labi)
  }
  return(list(data.list = data.list, samples = SampleSDV))
}

#################################################################
# Conquer
FuncConquer = function(data.list, x.test, dim, method = "SSANOVA", fast = TRUE,
                       mtype = list("cubic",c(0,1)), sigma = 500, lambda = NULL, basis = "default"){
  ## basis: default / full
  nnode = length(data.list)
  pred = NULL
  for(i in 1:nnode){
    datai = data.list[[i]]
    if(method == "SSANOVA"){
      typelist = vector("list", dim)
      for(l in 1:dim) typelist[[l]] = mtype
      if(dim == 1){
        names(typelist) = "x"
      } else{
        names(typelist) = paste0("x.", 1:dim)
      }
      if(basis == "default"){
        fit1 <- ssanova(y~., type=typelist, data=datai)
      } else{
        fit1 <- ssanova(y~., type=typelist, data=datai, id.basis=1:nrow(datai))
      }
      dTest = data.frame(x.test)
      if(dim == 1){
        names(dTest) = "x"
      } else{
        names(dTest) = paste0("x.", 1:dim)
      }
      pdi<-predict(fit1, newdata=dTest, se.fit=F)
    } else if(method == "KRR"){
      nn = nrow(datai)
      if(is.null(lambda)){
        lambda = nn^(-2/3)
      }
      if(fast){
        d = constructData(x = as.matrix(datai[, 1:dim]), y = datai$y) ## Structure data in CVST format
        krr_learner = constructFastKRRLearner()   ## Build the base learner
        params = list(kernel = 'rbfdot', sigma = sigma, lambda = lambda, nblocks = floor(nn^(1/3)))
        ## Function params; documentation defines lambda as '.1/getN(d)'
      } else{
        d = constructData(x = datai[, 1:dim], y = datai$y) ## Structure data in CVST format
        krr_learner = constructKRRLearner()   ## Build the base learner
        params = list(kernel = 'rbfdot', sigma = sigma, lambda = lambda)
        ## Function params; documentation defines lambda as '.1/getN(d)'
      }
      krr_trained = krr_learner$learn(d, params)
      ## Now to predict, format your test data, 'dTest', the same way as you did in 'd'
      if(fast){
        dTest = constructData(x = as.matrix(x.test[, 1:dim]), y = 0)
      } else{
        dTest = constructData(x = x.test[, 1:dim], y = 0)
      }
      pdi = krr_learner$predict(krr_trained, dTest)
      pdi = as.vector(pdi)
    } else if(method == 'RF'){
        typelist = vector("list", dim)
        for(l in 1:dim) typelist[[l]] = mtype
        names(typelist) = paste0("x", 1:dim)
        fit1 <- randomForest(y~., data=datai)
        dTest = data.frame(x.test)
        names(dTest) = paste0("x", 1:dim)
        pdi<-predict(fit1, newdata=dTest)
    }
    pred = rbind(pred, pdi)
  }
  pd = apply(pred, 2, mean)
  return(list(Pred = pd))
}


######################################################################################
###### true function

Tfunc1d = function(x, v1 = .05, v2 = 10){
  vib = .1/(abs(x - 1/2) + v1)
  y = vib * sin(pi*vib/v2)
  return(y)
}

Tfunc2d = function(x, v1 = .05, v2 = 10){
  vib = .1/(sqrt((x[,1]-1/2)^2 + (x[,2]-1/2)^2) + v1)
  y = vib * sin(pi*vib/v2)
  return(y)
}

l2dist = function(x, y) sqrt(sum((x-y)^2))
Tfunc_dim = function(x, v1 = .05, v2 = 10, dim = 4, cp = 0.5){
  c = rep(cp, dim)
  if(dim == 1){
    d = abs(x - cp)
  } else{
    d = apply(x, 1, l2dist, y = c)
  }
  vib = .1/(d + v1)
  y = vib * sin(pi*vib/v2)
  return(y)
}
# d <= 0.074, y > = .2

DataGenFunc = function(ngrid, dim, cp, peaks = 1, seed = 2023){
  set.seed(seed)
  if(dim==1){
    ###### data generating
    x = seq(0, 1, length = ngrid)
    x.gen = runif(ngrid)
  } else{
    x = sobol(n = 1e4, dim = dim, init = T, scrambling = 1)
    x.gen = matrix(runif(ngrid*dim), ncol = dim)
  }
  #y.t = Tfunc1d(x, v1 = .1, v2 = 1)
  if(peaks == 1){
      y.t = Tfunc_dim(x, dim = dim, cp = cp)
      y.gen.t = Tfunc_dim(x.gen, dim = dim, cp = cp)
  } else if(peaks == 2){
      y.t1 = Tfunc_dim(x, dim = dim, cp = cp)
      y.t2 = Tfunc_dim(x, dim = dim, cp = cp+.3)
      y.t = y.t1 + y.t2
      
      y.gen.t1 = Tfunc_dim(x.gen, dim = dim, cp = cp)
      y.gen.t2 = Tfunc_dim(x.gen, dim = dim, cp = cp+.3)
      y.gen.t = y.gen.t1 + y.gen.t2
  }
  if(dim == 1){
    y = y.gen.t + rnorm(ngrid, sd = .06)
  } else{
    y = y.gen.t + rnorm(ngrid, sd = .02)
  }
  #points(x, y, pch = 16, col = "gray")
    
  my.data = data.frame(x = x.gen, y = y)
  return(list(my.data = my.data, y.t = y.t, x.grid = x))
}

ckmat = function(x, range1, range2){
  c1 = x[1] >= range1[1]; c2 = x[1] <= range1[2]
  c3 = x[2] >= range2[1]; c4 = x[2] <= range2[2]
  return(c1 + c2 + c3 + c4)
}

DataGenFunc_Nonunif = function(ngrid, rate, bound, dim){
  if(dim==1){
    ###### data generating
    x = seq(0, 1, length = ngrid)
    #y.t = Tfunc(x, v1 = .1, v2 = 1)
    y.t = Tfunc1d(x)
    #y.t = sin(x*2*pi)
    #plot(x, y.t, type = "l")
    
    mark = which(y.t >= bound)
    xbound = c(x[min(mark)], x[max(mark)])
    
    n1 = floor(ngrid*(1-rate)/2)
    n0 = ngrid - 2*n1
    
    x.s11 = runif(n1, min = 0, max = xbound[1])
    x.s12 = runif(n1, min = xbound[2], max = 1)
    x.s2 = runif(n0, min = xbound[1], max = xbound[2])
    
    x.gen = c(x.s11, x.s2, x.s12)
    y.gen.t = Tfunc1d(x.gen)
    
    set.seed(2020)
    y = y.gen.t + rnorm(ngrid, sd = .06)
    points(x.gen, y, pch = 16, col = "gray")
    
    my.data = data.frame(x = x.gen, y = y)
  } else if(dim==2){
    x1 = x2 = seq(0, 1, length = 300)
    M = mesh(x1, x2)
    x = cbind(as.vector(M$x), as.vector(M$y))
    
    y.t = Tfunc2d(x)
    
    mark = which(y.t >= bound)
    x1bound = range(x[,1][mark])
    x2bound = range(x[,2][mark])
    
    n1 = floor(ngrid*(1-rate))
    n0 = ngrid - n1
    
    x11 = runif(n1)
    x12 = runif(n1)
    x1 = cbind(x11, x12)
    
    ck = apply(x1, 1, ckmat, range1 = x1bound, range2 = x2bound)
    rmlab = which(ck == 4)
    
    x.s1 = x1[-rmlab, ]
    
    x21 = runif(n0, min = x1bound[1], max = x1bound[2])
    x22 = runif(n0, min = x2bound[1], max = x2bound[2])
    x.s2 = cbind(x21, x22)
    
    x.gen = rbind(x.s1, x.s2)
    y.gen.t = Tfunc2d(x.gen)
    
    set.seed(2020)
    y = y.gen.t + rnorm(length(y.gen.t), sd = .02)
    
    my.data = data.frame(x1 = x.gen[,1], 
                         x2 = x.gen[,2], y = as.vector(y))
  } else{
    x = sobol(n = 1e4, dim = dim, init = T, scrambling = 1)
    #x = matrix(runif(1e6*dim), ncol = dim)
    
    y.t = Tfunc_dim(x, dim = dim)
    
    mark = which(y.t >= bound)
    xbound = range(x[,1][mark])
    
    n1 = floor(ngrid*(1-rate))
    n0 = ngrid - n1
    
    x1 = matrix(runif(n1*dim), ncol = dim)
    
    ck = apply(x1, 1, l2dist, y = rep(1/2, dim))
    rmlab = which(ck <= 0.074)
    
    if(length(rmlab)>0){
      x.s1 = x1[-rmlab, ]
    } else{
      x.s1 = x1
    }
    
    x.s2 = matrix(runif(n0*dim, min = xbound[1], max = xbound[2]), ncol = dim)
    
    x.gen = rbind(x.s1, x.s2)
    y.gen.t = Tfunc_dim(x.gen, dim = dim)
    
    sgm = sd(y.gen.t)
    
    set.seed(2020)
    y = y.gen.t + rnorm(length(y.gen.t), sd = .02)
    
    #testlab = sample(1:1e5, 1e4)
    
    my.data = data.frame(x = x.gen, y = as.vector(y))
  }
  return(list(my.data = my.data, y.t = y.t, x.grid = x))
}


############################################
## Our proposed method
############################################

RDLRT = function(my.data=my.data, x.test = x.test, seed=2023, nnode=100, basecut=1, 
                 breaks=10, downsample=F, sigma=100, method='KRR', downsamplesize=20){
  set.seed(seed)
  dim=ncol(my.data)-1
  sdv = SlDivide_repeat(my.data, nnode, basecut = basecut, breaks = breaks, downsample = F, downsamplesize)
  data.list.sdv = sdv$data.list
  lambda.sdv = nrow(sdv$data.list[[1]])^(-2/3)
  res.sdv = FuncConquer(data.list.sdv, x.test, dim = dim, method = "KRR",
                        sigma = sigma, lambda = lambda.sdv)
  pd = res.sdv$Pred
  return(pd)
}

