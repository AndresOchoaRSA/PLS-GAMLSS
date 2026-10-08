### se genera la distribución multivariada SEP
### se compara PLS-GAMLSS-SEP vs PLS clásico (Usando AIC y BIC)
### 20 Sept 2026
### propuesta beta= 0.4 y 0.8.



library(gamlss) ## Regresion flexible semi parametrica
library(chemometrics) ## Regresion PLS clasica
library(ggplot2)
source("pls_gamlss.R")
library(ggplot2)
library(gridExtra)

#### SEP


library(mixSPE)

mu1= c(0,0,0,0,0,0,0,0,0,0,0)
cov2= matrix(c(1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
               0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
               0.9,0.9,10,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
               0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,0.9,0.9,
               0.9,0.9,0.9,0.9,10,0.9,0.9,0.9,0.9,0.9,0.9,
               0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,0.9,0.9,0.9,
               0.9,0.9,0.9,0.9,0.9,0.9,10,0.9,0.9,0.9,0.9,
               0.9,0.9,0.9,0.9,0.9,0.9,0.9,10,0.9,0.9,0.9,
               0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1,0.9,0.9,
               0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,10,0.9,
               0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,1), 
             ncol=11, nrow=11, byrow=TRUE)

lambda1 = c(-2,-2,-2,-2,-2,2,2,2,2,2,2)


### SIMULACION

M <- 10
w <- c(50,100,150,200,250,300)

AICmSkewEMSEP <-matrix(0, nrow=M, ncol=length(w))
BICmSkewEMSEP <-matrix(0, nrow=M, ncol=length(w))

AICmClasic <-matrix(0, nrow=M, ncol=length(w))
BICmClasic <-matrix(0, nrow=M, ncol=length(w))


for(j in 1:length(w)){
  for(m in 1:M) {
    
    x1 = rspe(n = w[j], beta = 0.4, location = mu1, scale =
                cov2, psi = lambda1)
    X=x1[,2:11]
    y=x1[,1]
    
    ### Modelo PLS-GAMLSS-SEP
    
    PLS.skewSEP <- try(pls_gamlss(X, y, a=2, familyd=SEP))  ### a vaya hasta p o hasta el rango de X
    
    MODPLS1 <- pls1_nipals(X,y, a=2)
    modxx <- lm(y ~ MODPLS1$T)
    resumenModxx <- summary(modxx)
    
    AICmSkewEMSEP[m,j] <- PLS.skewSEP$AIC
    BICmSkewEMSEP[m,j] <- PLS.skewSEP$BIC
    
    AICmClasic[m,j] <-  AIC(modxx)
    BICmClasic[m,j] <- BIC(modxx)
    
  }
}



write.table(AICmSkewEMSEP, "AICmSkewEMSEP.txt") 
write.table(BICmSkewEMSEP, "BICmSkewEMSEP.txt") 
write.table(AICmClasic, "AICmClasic.txt") 
write.table(BICmClasic, "BICmClasic.txt") 


#### grafico AIC
#### 20 Sept 2026
#### SEP

M2 = M*2

# create a data frame
muestras=c(rep(" 50", each=M2),rep("100", each=M2),
           rep("150", each=M2), rep("200", each=M2),
           rep("250", each=M2), rep("300", each=M2)) 
methods=c(rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M),
            rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M),
            rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M))
AIC2=c(AICmSkewEMSEP[,1], AICmClasic[,1],
       AICmSkewEMSEP[,2], AICmClasic[,2],
       AICmSkewEMSEP[,3], AICmClasic[,3],
       AICmSkewEMSEP[,4], AICmClasic[,4],
       AICmSkewEMSEP[,5], AICmClasic[,5],
       AICmSkewEMSEP[,6], AICmClasic[,6])
data=data.frame(muestras, methods ,  AIC2)

# grouped boxplot
gg7 <- ggplot(data, aes(x=muestras, y=AIC2, fill=methods)) + 
  geom_boxplot()+
  labs(y="AIC", x="sample")


##### BIC SEP

M2 = M*2

# create a data frame
muestras=c(rep(" 50", each=M2),rep("100", each=M2),
           rep("150", each=M2), rep("200", each=M2),
           rep("250", each=M2), rep("300", each=M2)) 
methodPLS=c(rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M),
            rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M),
            rep("PLS-GAMLSS-SEP",each=M),rep("PLS",each=M))
BIC2=c(BICmSkewEMSEP[,1], BICmClasic[,1],
       BICmSkewEMSEP[,2], BICmClasic[,2],
       BICmSkewEMSEP[,3], BICmClasic[,3],
       BICmSkewEMSEP[,4], BICmClasic[,4],
       BICmSkewEMSEP[,5], BICmClasic[,5],
       BICmSkewEMSEP[,6], BICmClasic[,6])
data=data.frame(muestras, methods ,  BIC2)

# grouped boxplot
gg8 <- ggplot(data, aes(x=muestras, y=BIC2, fill=methods)) + 
  geom_boxplot()+
  labs(y="BIC", x="sample")


library(gridExtra)
grid.arrange(gg7,gg8,ncol=2)
