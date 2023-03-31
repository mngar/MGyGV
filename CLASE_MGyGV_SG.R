#clase MGyGV
#GWAS

# 1. INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
# 2. SETEAR UNA SEMILLA ASI TODOS TIENEN LOS MISMOS DATOS
# 3. SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS)
# 4. PONERLE NOMBRE A LOS INDIVIDUOS Y MARCADORES
# 5. RANKEAR MARCADORES ASOCIADOS Y EFECTOS: EXPORTAR COMO TABLA
# 6. REALIZAR ANALISIS GWAS CON DISTINTAS METODOLOGIAS Y COMPARAR RESULTADOS CON NUESTROS "QTLs REALES"

# instalacion de paquetes
#PAQUETE PARA REALIZAR SELECCION GENOMICA MEDIANTE RRBLUP
install.packages("rrBLUP")


#CARGA DE PAQUETES
library(simulMGF)
library(rrBLUP)

#setear la semilla (VAMOS A SIMULAR LOS MISMOS DATOS QUE UTILIZAMOS EN LA PRACTICA DE GWAS)
set.seed(1234)


#SIMULACION DE DATOS  (MISMOS PARAMETROS UTILIZADOS EN PRACTICA DE GWAS)
Nind <- 1000
Nmarkers <- 10000
Nqtl <- 50
Esigma <- .5
Pmean <- 25
Perror <- .25

simulN(Nind, Nmarkers, Nqtl, Esigma, Pmean, Perror)
 str(nsimout)

#los QTLs
  QTL <- cbind(nsimout$QTN, nsimout$Meffects)
  QTL <- cbind(QTL,abs(nsimout$Meffects))
  colnames(QTL) <- c("marker", "effect", "effabs")
  QTL <- as.data.frame(QTL)
  QTL <- QTL[order(-QTL$effabs),]

#matriz de datos genomicos
population <- nsimout$geno
colnames(population) <- c(paste("M", 1:Nmarkers,sep = ""))
rownames(population) <- c(paste("IND", 1:Nind,sep = ""))

#matriz de datos fenotipicos
popfeno <- nsimout$pheno
colnames(popfeno) <- "PHENO"
rownames(popfeno) <- c(paste("IND", 1:Nind,sep = ""))
feno <- data.frame ( IND = rownames(popfeno),
                     Pheno = popfeno)



#mapa
map <- data.frame (SNP = colnames(population),
                  Chromosome = c(rep(1,(Nmarkers/5)),rep(2,(Nmarkers/5)),rep(3,(Nmarkers/5)),rep(4,(Nmarkers/5)),rep(5,(Nmarkers/5))),
                  Position = c(1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5))
                  )






##########
##  SG  ##
##########
#setear el directorio de trabajo
setwd("D:/MGyGV/R/SG")

#esquema de validacion cruzada
train <- round(0.9*Nind, digits = 0)          # porcentaje pob entrenamiento
test <- Nind-train
cv <- 10                                        # cantidad de validaciones cruzadas
if(isTRUE(test+train == Nind) == TRUE) {
xtrain <- rep(1,train)
xtest <- rep(2,test)
xcv <- c(xtrain, xtest)
index <- matrix(nrow = Nind, ncol = cv, NA)
for (i in 1:cv) {
index[,i] <- sample(xcv)
}
} else {
print("train+test es distinto de Nind, adecuar valores")
}







ntest <- test
nind  <- Nind

sets1  <- index
datas2 <- population
datas3 <- feno

#MATRIZ DATOS GENOTIPICOS Y FENOTIPICOS
datas <- cbind(population, feno$PHENO)


for(fold in 1:10){
 correlations <- matrix(NA,1,6)
 datos        <- matrix(NA,1,6)
################################################################################
  itrain <- which(sets1[,fold]==1)
  itest  <- which(sets1[,fold]==2)
  test   <- datas[itest,]
  train  <- datas[itrain,]
  Xtest  <- test[,-ncol(test)]
  Ytest  <- test[,ncol(test)]
  Xtrain <- train[,-ncol(test)]
  Ytrain <- train[,ncol(test)]
################################################################################
indice <- 1
for (columna1 in 1:(ncol(Xtrain)-1)){
    for (columna2 in (columna1+1):ncol(Xtrain)){
        dd <- sum(abs(Xtrain [,columna1]- Xtrain [,columna2]))
    if(dd==0){
        indice <- c(indice,columna1)
}}}
##  RR-BLUP       ##############################################################
  X1     <- rbind(Xtrain,Xtest)
  X1     <- X1[,-indice]
  y1     <- c(Ytrain,Ytest)
  yNa1   <- y1
  train1 <- 1:nrow(Xtrain)
  f      <-  nrow(Xtrain)+1
  pred1  <- f:nrow(X1)
  ans    <- mixed.solve(y=y1[train1],Z=Xtrain) #By default K = I
  intercepto <- rep(ans$beta,ntest)

efectoRR  <- ans$u
efectoRR2 <- abs(efectoRR)
num <- length(efectoRR)
rr  <- cbind(1:num,efectoRR,efectoRR2)
newdata <- rr[order(-efectoRR2),]
m <- matrix(data=NA,nrow=length(efectoRR),ncol=3)
for(i in 1:length(efectoRR)) {
      if(newdata[i,2] > median(efectoRR)){
      m[i,1] <- newdata[i,1]}
}

yTestHat   <- intercepto + Xtest %*% ans$u
accuracyRR <- cor(Ytest,yTestHat)
datos      <- paste(fold,efectoRR,efectoBL,efectoRF,sep = ",")
efectos    <- cbind(efectoRR,efectoBL,efectoRF)

write.table(datos,file="SALIDA_EFECTOS_RR.csv",row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
datos = paste(fold, accuracyRR,sep = ",")
write.table(datos,file="salida_PRECISION_SG.csv",row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")
    }





