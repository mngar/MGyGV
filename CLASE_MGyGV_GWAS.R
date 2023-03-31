#clase MGyGV
#GWAS

# 1. INSTALAR Y CARGAR LOS PAQUETES NECESARIOS
# 2. SETEAR UNA SEMILLA ASI TODOS TIENEN LOS MISMOS DATOS
# 3. SIMULAR MATRIZ DE DATOS GENOMICOS Y FENOTIPICOS (QTLs y EFECTOS)
# 4. PONERLE NOMBRE A LOS INDIVIDUOS Y MARCADORES
# 5. RANKEAR MARCADORES ASOCIADOS Y EFECTOS: EXPORTAR COMO TABLA
# 6. REALIZAR ANALISIS GWAS CON DISTINTAS METODOLOGIAS Y COMPARAR RESULTADOS CON NUESTROS "QTLs REALES"

# instalacion de paquetes
#PAQUETE DE SIMULACION DE DATOS
install.packages("simulMGF")

#PAQUETES PARA CHECKEAR NORMALIDAD DEL FENOTIPO (GRAFICA Y ANALITICAMENTE)
install.packages("ggplot2")
library("nortest")

#PAQUETE DE ANALISIS DE ASOCIACION
install.packages("remotes")
remotes::install_github("jiabowang/GAPIT3")


#SOLO SI LO ANTERIOR NO FUNCIONA PORQUE NO TIENEN CARGADO RTOOLS:
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#CARGA DE PAQUETES
library(simulMGF)
library(ggplot2)
library(ggplot2)
#library(GAPIT)
library(GAPIT3)

#setear la semilla (para que los resultados puedan reproducirse exactamente)
set.seed(1234)


#SIMULACION DE DATOS
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
#armamos la matriz de datos genomicos como la requiere GAPIT3
myGD <- as.data.frame(population))
myGD <- cbind(feno$IND,myGD)
colnames(myGD)[1] <- "IND"
#write.csv(myGD, "myGD.csv")
#myGD <- as.data.frame(read.csv("myGD.csv", header = T))
#myGD <- read.csv("myGD.csv", header = T)
#myGD <- myGD[,-1]


#matriz de datos fenotipicos
popfeno <- nsimout$pheno
colnames(popfeno) <- "PHENO"
rownames(popfeno) <- c(paste("IND", 1:Nind,sep = ""))
feno <- data.frame ( IND = rownames(popfeno),
                     Pheno = popfeno)


#mapa
map <- data.frame (SNP = colnames(population),
                  Chromosome = c(rep(1,(Nmarkers/5)),rep(2,(Nmarkers/5)),rep(3,(Nmarkers/5)),rep(4,(Nmarkers/5)),rep(5,(Nmarkers/5))),
                  Position = c(1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5),1:(Nmarkers/5))                             #Jiabo Wang: GAPIT requiere que el cromosoma 1 tenga >100 marcadores en cromosoma 1.
                  )




#PRUEBAS DE NORMALIDAD
#CHECKEO VISUAL

# HISTOGRAMA Y CUARVA NORMAL: Consiste en representar los datos mediante un histograma
#y superponer la curva de una distribución normal con la misma media y desviación estándar
#que muestran los datos.

ggplot(data = feno, aes(x = PHENO)) +
  geom_histogram(aes(y = ..density.., fill = ..count..)) +
  scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
  stat_function(fun = dnorm, colour = "firebrick",
                args = list(mean = mean(feno$PHENO),
                            sd = sd(feno$PHENO))) +
  ggtitle("Histograma + curva normal teórica") +
  theme_bw()

#Gráfico de cuantiles teóricos (QQplot)
#Consiste en comparar los cuantiles de la distribución observada con los cuantiles
#teóricos de una distribución normal con la misma media y desviación estándar que los datos.
#Cuanto más se aproximen los datos a una normal, más alineados están los puntos entorno a la
#recta.
qqnorm(feno$PHENO, pch = 19, col = "gray50")
qqline(feno$PHENO)


#CHECKEO ANALITICO

#test  Lilliefors: modificación del  Kolmogorov-Smirnov.
#El test Lilliefors asume que la media y varianza son desconocidas, estando especialmente desarrollado
#para contrastar la normalidad. Es la alternativa al test de Shapiro-Wilk cuando el número de
#observaciones es mayor de 50. La función lillie.test() del paquete nortest permite aplicarlo.
library("nortest")
lillie.test(x = feno$PHENO)

#D es el estadístico que hay que mirar que sea <0.05


##########
## GWAS ##
##########
#setear el directorio apropiado para los archivos de salida
setwd("D:/MGyGV/R/GWAS")


#GWAS CON DISTINTOS MODELOS
myGAPIT <- GAPIT(
  Y=feno,
  GD=myGD,
  GM=map,
  model=c("GLM","MLM","SUPER","MLMM","FarmCPU","Blink"),# choose model
  PCA.total=0,                                          # set total PCAs
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T,
  Geno.View.output=FALSE
)


#Comparar marcadores asociados con los QTLs "reales"






