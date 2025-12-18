pls_gamlss <- function (X, y, a, familyd) 
{
  Xh <- scale(X, center = TRUE, scale = TRUE)
  Xh.res <- Xh
  yh <- y
  datah <- data.frame(yh,Xh)
  n <- nrow(Xh)
  pXo <- ncol(Xh)
  ww <- NULL
  Tx <- NULL
  P <- NULL
  C <- NULL
  W <- NULL
  #if(a > rank(X),"Stop: 'a' must be between 1 and the rank of X","ok") 
  
  for (h in 1:a) {
    if(h<=1) {
      whi <- 0
      for (j in 1:pXo) {
        ww <- gamlss(yh~Xh[,j], family=familyd, data=datah)
        whi[j] <- (coef(ww)[h+1])
      } 
    }  
    else  {
      whi <- 0
      for (j in 1:pXo) {
        ww <- gamlss(yh~Tx+Xh[,j],family=familyd, data=datah)
        whi[j] <- (coef(ww)[h+1])
      }  
    }
    wh <- whi /as.vector(sqrt(t(whi) %*% whi))
    th <- Xh.res %*% wh
    ch <- t(yh)%*% th/as.vector(t(th) %*% th)
    ph <- 0
    for(j in 1:pXo) {
      ph[j] <- t(Xh[,j]) %*% th/as.vector(t(th) %*% th)
    }
    Xh.res <- Xh.res - (th %*% t(ph))
    Tx <- cbind(Tx, th)
    P <- cbind(P, ph); row.names(P) <- colnames(X)
    C <- c(C, ch)
    W <- cbind(W, wh); row.names(W) <- colnames(X)
  }
  datah2 <- data.frame(yh,Tx)
  modPLS <- gamlss(yh ~ Tx, family=familyd,data=datah2)  ##  BETA serian los q
  q <- coef(modPLS)[2:(a+1)]
  sumaryModPLS <- summary(modPLS)
  b.PLS <- W %*% q ; row.names(b.PLS) <- colnames(Xh)
  yPred2 <- fitted(modPLS)
  ress <- sum( (y - yPred2)^2 ) /n
  sst <- sum( (y - mean(y))^2 ) /n
  NumPara <- length(coef(modPLS)) ## numero de parametros del Modelo PLS
  list(P = P, Tx = Tx, W = W,
       AIC= modPLS$aic, 
       BIC= modPLS$sbc,
       resid= modPLS$residuals, 
       b.PLS=b.PLS, q.hat=q, 
       resumen=sumaryModPLS,
       ressH=ress, yPred2= yPred2)
}
