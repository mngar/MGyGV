#Tipo de ensayo	STP (Single Tree Plot)
#Número de sitios	4
#Número de clones	100
#Réplicas clonales por genotipo	3 por sitio
#Número de árboles por clon y sitio	3
#Número total de árboles	1.200
#Variable de interés	Diámetro a la altura del pecho (DAP, cm)
#Unidad experimental	Árbol individual

#Estimaremos el BLUP del efecto genético de cada clon a partir de sus observaciones fenotípicas. 
#También incorporaremos un efecto aleatorio de interacción clon x sitio para representar que el 
#comportamiento relativo de los clones puede variar entre ambientes.

#Modelo estadístico
#	Yijk = mu + Si + Cj + (SC)ij + Eijk

# Donde,
#Yijk: DAP observado en el árbol k, del clon j, en el sitio i.
#mu: media general.
#Si: efecto fijo del sitio \(i\).
#Cj: efecto genético aleatorio del clon \(j\).
#SCij: efecto aleatorio de la interacción clon × sitio.
#Eijk: error residual asociado al árbol individual.

#suponemos que los efectos aleatorios siguen distribuciones normales

#  Cj ~ N (0, VARc)

#  (SC)ij ~ N (0, VARcs)

#  Eijk ~ N (0,VARe)

#El BLUP individual del clon será la predicción de su efecto genético Cj, considerando 
#simultáneamente sus observaciones en los cuatro sitios, las diferencias ambientales y la 
#variabilidad residual.

# ============================================================
# 1. Preparación
# ============================================================

# Instalar una sola vez si es necesario:
# install.packages("lme4")
# install.packages("ggplot2")
# install.packages("dplyr")

library(lme4)
library(ggplot2)
library(dplyr)

set.seed(123)

# ============================================================
# 2. Definición del ensayo
# ============================================================

n_clones <- 100
n_sitios <- 4
n_ramets <- 3

clones <- paste0("C", sprintf("%03d", 1:n_clones))
sitios <- paste0("S", 1:n_sitios)

# Una fila por árbol individual
datos <- expand.grid(
  clon = clones,
  sitio = sitios,
  ramet = 1:n_ramets
)

# Factores
datos$clon <- factor(datos$clon, levels = clones)
datos$sitio <- factor(datos$sitio, levels = sitios)

# Comprobación del número de árboles
nrow(datos)

# 100 clones x 4 sitios x 3 ramets = 1200 árboles

###
### definiremos los componentes de variación antes de generar el DAP.
###

# ============================================================
# 3. Parámetros de simulación
# ============================================================

media_general <- 20

# Efectos fijos de los sitios (cm)
efecto_sitio <- c(
  S1 = -2,
  S2 =  0,
  S3 =  1,
  S4 =  3
)

# Desviaciones estándar de los efectos aleatorios (cm)
sd_clon <- 2.0
sd_interaccion <- 1.2
sd_residual <- 2.5

# ============================================================
# 4. Efectos genéticos de los clones
# ============================================================

efectos_clon <- data.frame(
  clon = factor(clones, levels = clones),
  efecto_genetico = rnorm(
    n_clones,
    mean = 0,
    sd = sd_clon
  )
)

# ============================================================
# 5. Efectos de la interacción clon x sitio
# ============================================================

efectos_cs <- expand.grid(
  clon = clones,
  sitio = sitios
)

efectos_cs$efecto_cs <- rnorm(
  nrow(efectos_cs),
  mean = 0,
  sd = sd_interaccion
)

efectos_cs$clon <- factor(
  efectos_cs$clon,
  levels = clones
)

efectos_cs$sitio <- factor(
  efectos_cs$sitio,
  levels = sitios
)

# ============================================================
# 6. Construcción del fenotipo
# ============================================================

datos <- datos |>
  left_join(efectos_clon, by = "clon") |>
  left_join(efectos_cs, by = c("clon", "sitio"))

datos$efecto_ambiental <- efecto_sitio[
  as.character(datos$sitio)
]

datos$error <- rnorm(
  nrow(datos),
  mean = 0,
  sd = sd_residual
)

datos$DAP <- media_general +
  datos$efecto_ambiental +
  datos$efecto_genetico +
  datos$efecto_cs +
  datos$error

# Inspección de los datos
head(datos)

# Resumen del DAP
summary(datos$DAP)

# Verificar el número de árboles
nrow(datos)



###
### Ajuste del modelo mixto y estimación de los BLUP
###

# ============================================================
# 7. Ajuste del modelo mixto
# ============================================================

modelo <- lmer(
  DAP ~ sitio +
    (1 | clon) +
    (1 | clon:sitio),
  data = datos,
  REML = TRUE
)

# Resultados generales
summary(modelo)

# Componentes de varianza estimados
VarCorr(modelo)

###
### La fórmula del modelo es:
###
### DAP ~ sitio + (1 | clon) + (1 | clon:sitio)

# ============================================================
# 8. Extraer los BLUP individuales de los clones
# ============================================================

blup_clones <- ranef(modelo)$clon |>
  data.frame() |>
  tibble::rownames_to_column("clon")

# Renombrar la columna del efecto aleatorio
names(blup_clones)[2] <- "BLUP"

# Ordenar de mayor a menor valor genético predicho
blup_clones <- blup_clones |>
  arrange(desc(BLUP))

# Mostrar los 10 clones con mayor BLUP
head(blup_clones, 10)



#El BLUP representa la desviación genética clonal predicha respecto de la media genética
# de la población de referencia, expresada en centímetros de DAP.

#Por ejemplo, si un clon tiene un BLUP de +2 cm, su efecto genético clonal predicho es 2 cm 
#superior a la media genética de referencia, una vez considerados los efectos incluidos en el modelo.

#Si otro clon presenta un BLUP de -1 cm, su efecto genético predicho es inferior a esa media.

#Importante: los BLUP no son los valores fenotípicos medios de los clones. 
#Son predicciones de sus efectos aleatorios y están sujetos a contracción hacia cero (shrinkage), 
#en función de la cantidad y calidad de la información disponible.

# ============================================================
# 9. Comparación entre efectos genéticos reales y predichos
# ============================================================

comparacion <- efectos_clon |>
  rename(
    efecto_real = efecto_genetico
  ) |>
  left_join(blup_clones, by = "clon")

# Correlación entre valores reales y predichos
cor(
  comparacion$efecto_real,
  comparacion$BLUP
)

# Gráfico
ggplot(
  comparacion,
  aes(x = efecto_real, y = BLUP)
) +
  geom_point(alpha = 0.7) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "Efectos genéticos reales y BLUP",
    subtitle = "Datos simulados de un ensayo STP multisitio",
    x = "Efecto genético real (cm)",
    y = "BLUP del clon (cm)"
  ) +
  theme_minimal()


#¿Cómo interpretar el gráfico?
#Cada punto representa un clon.
#- Puntos próximos a la diagonal: el modelo predice correctamente el efecto genético del clon.
#- Puntos alejados de la diagonal: existe una diferencia importante entre el efecto real y el predicho.
#- Pendiente aparente menor que uno: los BLUP suelen mostrar menor dispersión que los efectos 
#genéticos reales, debido a la contracción hacia la media.
#La correlación permite evaluar la capacidad del modelo para ordenar los clones, mientras que la 
#proximidad a la diagonal permite evaluar la concordancia entre los valores reales y predichos.

# ============================================================
# 10. Gráfico de los BLUP
# ============================================================

ggplot(
  blup_clones,
  aes(
    x = reorder(clon, BLUP),
    y = BLUP
  )
) +
  geom_point() +
  coord_flip() +
  labs(
    title = "BLUP de los clones",
    x = "Clon",
    y = "Efecto genético clonal predicho (cm)"
  ) +
  theme_minimal()


###
### ¿Qué sucede si queremos evaluar cada clon en cada sitio?
### Hasta ahora estimamos un BLUP por clon, utilizando la información de los cuatro sitios.
### Sin embargo, el modelo también permite obtener las predicciones de la interacción clon X sitio.
###

# ============================================================
# 11. Extraer los efectos de interacción clon x sitio
# ============================================================

blup_interaccion <- ranef(modelo)[["clon:sitio"]] |>
  data.frame() |>
  tibble::rownames_to_column("clon_sitio")

names(blup_interaccion)[2] <- "BLUP_interaccion"

head(blup_interaccion)



# Predicciones que incluyen los efectos fijos y aleatorios
datos$DAP_predicho <- predict(
  modelo,
  newdata = datos,
  re.form = NULL
)

head(
  datos[, c("clon", "sitio", "ramet",
            "DAP", "DAP_predicho")]
)

#Estas predicciones son útiles para comparar el comportamiento esperado de los clones entre sitios.



#El flujo de trabajo es el siguiente:
#1. Organizar los datos: cada fila corresponde a un árbol individual con su clon, sitio y DAP.
#2. Ajustar el modelo mixto: estimar los componentes de varianza y los efectos fijos.
#3. Obtener los BLUP: extraer los efectos aleatorios de clon.
#4. Ordenar los clones: identificar los valores genéticos clonales predichos más altos.
#5. Evaluar las predicciones: comparar los BLUP con los efectos reales en la simulación.
#6. Analizar la interacción clon × sitio: estudiar si los clones presentan respuestas diferenciales entre ambientes.
#La idea central para transmitir a los estudiantes es que el BLUP permite combinar la información 
#fenotípica de los árboles con una estimación de la variación genética y residual para obtener predicciones 
#de los efectos genéticos de los clones, teniendo en cuenta la estructura del ensayo.

#Si posteriormente se quisiera utilizar este ejemplo para un programa de mejoramiento forestal real, 
#el siguiente paso sería adaptar el modelo a la estructura espacial del ensayo y evaluar si los valores 
#genéticos de interés son clonales, aditivos o específicos de cada sitio.