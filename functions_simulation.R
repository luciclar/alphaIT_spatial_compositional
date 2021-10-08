################# FUNCTIONS FOR GEOSTATISTICS FOR CODA ##############

######## ALPHA-IT TRANSFORMATION ##########
alpha.IT <- function(x, alpha){ 
  if (alpha == 0) {
    z <- as.data.frame(ilr(x,  V=t(helm(ncol(x)))))
  }
  else {
    D <- ncol(x)
    m.alpha = 1/D*sum(x**alpha)
    u <- (t(x**alpha) - m.alpha)/alpha
    H = helm(D)
    z <- t(H %*% u)
  }
  
  return (z)
}


######## INVERSE ALPHA-IT (NUMERICAL) ##########

alpha.IT.inv = function(z, alpha){
  X = apply(z,MARGIN = 1,FUN = alphaIT1.inv,alpha=alpha)
  return(t(X))
}

alphaIT1.inv = function(z,alpha){
  f = function(y,alpha,b){
    yy = exp(y)/sum(exp(y))
    y.alpha = yy^alpha/alpha
    m.alpha = mean(y.alpha)
    Q = sum((y.alpha-m.alpha-b)^2)
    return(Q)
  }
  D = length(z)+1
  b = as.vector(t(helm(D))%*%z)
  x.i = log(rep(1/D,D))
  x.o = nlminb(x.i, f, alpha=alpha,b=b,lower = 0, upper = Inf)$par
  x.out = exp(x.o)/sum(exp(x.o))
  return(x.out)
}


######## SCALING FOR ALPHA-IT CODOMAIN ##########

inside_shape <- function (alfa, data) {
  D <- ncol(data)+1
  H <- helm(D)
  inside = 0
  s = 1
  while (inside == 0) {
    ok = 1
    for (i in 1:nrow(data)) {
      result <- try(alpha.IT.inv(z[i,],alfa), TRUE)
      if (typeof(result) == "list") {  ## use only the simulations that work (otherwise typeof(result)="character")
        if (sum(result[[1]] == 0) ==1 && result[[2]] > 10^-7) {
          ok = 0
          break
        }
      }
      else {
        ok = 0
        break
      }
    }
    if (ok == 0) {
      s = s *0.99
      # print(s)
      data.scaled = data.frame(scale(data, center = TRUE, scale=c(1/s, 1/s)))
      data[,1] = data.scaled[,1] + me[1]
      data[,2] = data.scaled[,2] + me[2]
      
    }
    else {
      inside = 1
    }
  }
  return(list(s, data))
}


####### LOG-LIKELIHOOD ALPHA-IT ########

LogLik.alphaIT =function(X,alpha){
  # Computes the log-likelihood for comositional data 
  # according to alpha-IT meta Gaussian model
  # Limited to D=3 or D=4
  D = dim(X)[2] ; H = helm(D)
  
  # First, the likelihood is computed on interior points
  f.int  = function(x){sum(x>0)==length(x)}
  id.int = apply(X,FUN = f.int,MARGIN = 1)
  X.int  = X[id.int,]; X.bord = X[!id.int,]
  out = LL.alphaIT(alpha,X.int)
  
  # Loop on parts 
  for (i in 1:D){
    # First, data that have one part=0 only
    Xi = X.bord[X.bord[,i]==0,-i]
    if (!is.null(dim(Xi))){
      if (dim(Xi)[1]>D){
        id.i = apply(Xi,FUN = f.int,MARGIN = 1)
        if (sum(id.i)>(D-1)){
          out = out + LL.alphaIT(alpha,Xi[id.i,])
          # Now, data with txo parts = 0. Note, since D=4, if three parts are 0, 
          # the last part is 1. In this case, Log-Likelihood is 0 for all values of alpha  
          if (D>3){
            for (j in 1:(D-1)){
              Xij = Xi[Xi[,j]==0,-j]
              if(!is.null(dim(Xij))){
                if (dim(Xij)[1] > (D-1)){
                  id.ij = apply(Xij,FUN = f.int,MARGIN = 1)
                  if (sum(id.ij)>(D-2)){
                    out = out + 0.5*LL.alphaIT(alpha,Xij[id.ij,])
                  }
                }
              }
            }
          }
          # End loop on two parts being 0  
        }
      }
    }
  }
  return(out)
}

LL.alphaIT = function(alpha,X){
  LL = 0
  if (any(X==0)){
    stop("LogLik.alphaIT: no values equal to 0 allowed")
  }
  else{
    if(!is.null(dim(X))  & (dim(X)[1]> dim(X)[2]) ){
      ns = dim(X)[1]
      Z = alpha.IT(X,alpha)
      m = apply(Z,MARGIN = 2,FUN = mean)
      S = var(Z)
      prec = solve(S)
      LL = 0.5*ns*determinant(S,logarithm=TRUE)$modulus[1]
      LL = LL + 0.5*sum(apply(Z,MARGIN = 1,FUN = LQ1,m=m,prec=prec))
      LJ = LogJa.alphaIT(X,alpha)
      LL = as.numeric(LL) - LJ
    }
    return(LL)
  }
}


LQ1 = function(z,m,prec){
  out = t(z-m)%*%prec%*%(z-m)
  return(as.numeric(out))
}

LogJa.alphaIT = function(X,alpha){
  out = apply(X, MARGIN = 1, FUN = LogJa1,alpha=alpha)
  return(sum(out))
}

LogJa1 = function(x,alpha){
  if (any(x==0)){
    stop("LogJa1: no values equal to 0 allowed")
  }
  else{
    D = length(x)
    H = helm(D)
    y = x
    y[x>0] = x[x>0]^(alpha-1)
    J = matrix(0,ncol=D-1,nrow=D-1)
    for (i in 1:(D-1)){
      for (j in 1:(D-1)){
        J[i,j] = H[i,j]*y[j]-H[i,D]*y[D]
      }
    }
    out=determinant(J)$modulus
    return(as.numeric(out))
  }
}
