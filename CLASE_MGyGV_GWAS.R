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
install.packages("nortest")

#PAQUETE DE ANALISIS DE ASOCIACION
install.packages("remotes")
remotes::install_github("jiabowang/GAPIT3")


#SOLO SI LO ANTERIOR NO FUNCIONA PORQUE NO TIENEN CARGADO RTOOLS:
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#CARGA DE PAQUETES
library(simulMGF)
library(ggplot2)
library(nortest)
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


#matriz de datos fenotipicos
popfeno <- nsimout$pheno
colnames(popfeno) <- "PHENO"
rownames(popfeno) <- c(paste("IND", 1:Nind,sep = ""))
feno <- data.frame ( IND = rownames(popfeno),
                     Pheno = popfeno)
#armamos la matriz de datos genomicos como la requiere GAPIT3
myGD <- as.data.frame(population)
myGD <- cbind(feno$IND,myGD)
colnames(myGD)[1] <- "IND"
#write.csv(myGD, "myGD.csv")
#myGD <- as.data.frame(read.csv("myGD.csv", header = T))
#myGD <- read.csv("myGD.csv", header = T)
#myGD <- myGD[,-1]


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

QTL$SNP <- paste0("M",QTL$marker)

#blink
blink <- read.csv("GAPIT.Blink.PHENO.GWAS.Results.csv",header = T)
head(blink)
ablink <- blink[,c(1,9)]
head(ablink)
ablink <- ablink[ablink$FDR_Adjusted_P.values <0.05,]
asoc_blink <- merge(QTL, ablink, by = "SNP")

#FarmCPU
FarmCPU <- read.csv("GAPIT.FarmCPU.PHENO.GWAS.Results.csv",header = T)
head(FarmCPU)
aFarmCPU <- FarmCPU[,c(1,9)]
head(aFarmCPU)
aFarmCPU <- aFarmCPU[aFarmCPU$FDR_Adjusted_P.values <0.05,]
asoc_FarmCPU <- merge(QTL, aFarmCPU, by = "SNP")

#MLMM
MLMM <- read.csv("GAPIT.MLMM.PHENO.GWAS.Results.csv",header = T)
head(MLMM)
aMLMM <- MLMM[,c(1,9)]
head(aMLMM)
aMLMM <- aMLMM[aMLMM$FDR_Adjusted_P.values <0.05,]
asoc_MLMM <- merge(QTL, aMLMM, by = "SNP")

#SUPER
SUPER <- read.csv("GAPIT.SUPER.PHENO.GWAS.Results.csv",header = T)
head(SUPER)
aSUPER <- SUPER[,c(1,9)]
head(aSUPER)
aSUPER <- aSUPER[aSUPER$FDR_Adjusted_P.values <0.05,]
asoc_SUPER <- merge(QTL, aSUPER, by = "SNP")

#MLM
MLM <- read.csv("GAPIT.MLM.PHENO.GWAS.Results.csv",header = T)
head(MLM)
aMLM <- MLM[,c(1,9)]
head(aMLM)
aMLM <- aMLM[aMLM$FDR_Adjusted_P.values <0.05,]
asoc_MLM <- merge(QTL, aMLM, by = "SNP")

#GLM
GLM <- read.csv("GAPIT.GLM.PHENO.GWAS.Results.csv",header = T)
head(GLM)
aGLM <- GLM[,c(1,9)]
head(aGLM)
aGLM <- aGLM[aGLM$FDR_Adjusted_P.values <0.05,]
asoc_GLM <- merge(QTL, aGLM, by = "SNP")

#The statistical power of a hypothesis test is the probability of detecting an effect,
#if there is a true effect present to detect.

#Entonces podemos estimar el poder estadístico de cada evaluación en nuestros datos
#simulados de la siguiente manera (ya que conocemos los marcadores con asociación real
#al caracter, imposible en datos reales):
#pow <- Marcadores_Asociados_Reales/Nqtl

pow_blink <- dim(asoc_blink)[1]/Nqtl
pow_FarmCPU <- dim(asoc_FarmCPU)[1]/Nqtl
pow_MLMM <- dim(asoc_MLMM)[1]/Nqtl
pow_SUPER <- dim(asoc_SUPER)[1]/Nqtl
pow_MLM <- dim(asoc_MLM)[1]/Nqtl
pow_GLM <- dim(asoc_GLM)[1]/Nqtl

data <- matrix(c(
dim(ablink)[1]  ,   dim(asoc_blink)[1], dim(ablink)[1]-dim(asoc_blink)[1],     pow_blink,
dim(aFarmCPU)[1], dim(asoc_FarmCPU)[1], dim(aFarmCPU)[1]-dim(asoc_FarmCPU)[1], pow_FarmCPU,
dim(aMLMM)[1]   ,    dim(asoc_MLMM)[1], dim(aMLMM)[1]-dim(asoc_MLMM)[1],       pow_MLMM,
dim(aSUPER)[1]  ,   dim(asoc_SUPER)[1], dim(aSUPER)[1]-dim(asoc_SUPER)[1],     pow_SUPER,
dim(aMLM)[1]    ,     dim(asoc_MLM)[1], dim(aMLM)[1]-dim(asoc_MLM)[1],         pow_MLM,
dim(aGLM)[1]    ,     dim(asoc_GLM)[1], dim(aGLM)[1]-dim(asoc_GLM)[1],         pow_GLM), ncol=4, byrow=TRUE)

# specify the column names and row names of matrix
colnames(data) = c('associatedSNP', 'TRUEassociatedSNP', 'FALSEpositives', 'Power')
rownames(data) <- c('blink', 'FarmCPU', 'MLMM', 'SUPER', 'MLM', 'GLM')

final=as.table(data)
final





