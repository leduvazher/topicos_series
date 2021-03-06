######################################################################
#   An�lisis de intervenci�n - nivel de aprobaci�n de G. Bush        #
#   despu�s del ataque a las Torres Gemelas (9/11).                  #
######################################################################

# Los datos en el archivo BUSHJOB.DTA contienen los niveles de 
# aprobaci�n del G. Bush de enero del 2001 a Febrero del 2005
#
# El objetivo es medir el impacto en el nivel de aprovaci�n de Bush
# despu�s de 9/11
#
# Para leer los datos en formato .DTA se debe instalar el paquete
# "foreign" 
#

# Limpiar el ambiente
#
rm(list=ls())

# Establecer working directory
setwd("C:\\Users\\miguel.villalobos\\OneDrive - Red de Universidades An�huac\\DataSets\\SeriesTiempo")

# Instalar paquetes y librer�as requeridos

if(!require(tseries)) {install.packages("tseries")}
library(tseries)
if(!require(astsa)) {install.packages("astsa")}
library(astsa)
if(!require(xts)) {install.packages("xts")}
library(xts)
if(!require(forecast)) {install.packages("forecast")}
library(forecast)
if(!require(readxl)) {install.packages("readxl")}
library(readxl)
if(!require(ModelMetrics)) {install.packages("ModelMetrics")}
library(ModelMetrics)
if(!require(TSstudio)) {install.packages("TSstudio")}
library(TSstudio)
if(!require(MASS)) install.packages("MASS")
library(MASS)
if(!require(foreign)) install.packages("foreign")
library(foreign)
if(!require(TSA)) install.packages("TSA")
library(TSA)
#
# Funci�n para obtener la gr�fica 4 en de los residuales
#         Esta funci�n toma como entrada los valores residuales
#         del modelo genera la gr�fica Normal de los residuales, 
#         residuales vs tiempo y las funciones de autocorrelaci�n y
#         autocorrelaci�n parcial de los residuales.
#
CuatroEnUno = function(residuales){
  par(mfrow=c(2,2), oma=c(0,0,0,0))
  qqnorm(residuales, xlab="Percentiles te�ricos",
         ylab="Percentiles Observados", 
         main = "Grafica Normal de los Residuales")
  qqline(residuales, distribution=qnorm)
  plot(residuales, type="l", xlab="Tiempo", ylab = "Residuales",
       main = "Residuales vs Tiempo")
  acf(residuales, lag.max=25, type="correlation", 
      main = "FAC de los Residuales")
  acf(residuales, lag.max=25, type="partial", 
      main = "FACP de los Residuales")
  par(mfrow=c(1,1))
}

#funcion para llevar a cabo la prueba box-ljung de los residuales

Box.Ljung.Test = function (z,k) {
  lags = seq(1:k)
  pvalues = rep(0,k)
  for (i in c(1:k)) pvalues[i] = Box.test(z,lag=i)[3]
  plot(x=lags, y=pvalues, ylim=c(0,1), main = "Prueba Box.Ljun")
  abline(h=0.05, col="red")
}

####





#cargar datos y graficar series

bush <- read.dta("BUSHJOB.DTA")
names(bush)
head(bush)

View(bush)

fecha = as.Date(paste(bush$year,bush$month,"28",sep="-"),
                format="%Y-%m-%d")
BApp = ts(bush$approve, start=c(2001,01), freq=12)
tsplot(x=fecha, y=BApp)
abline(v=fecha[9], col="red")

#identificar proceso arima
acf2(BApp) #se fijan pareciera que la funcion de autocorr decae lentamente, hay un lag significativo
d1BApp = diff(BApp, lag=1)
acf2(d1BApp)

# Dadas las caracter�sticas de la ACF y PACF nos lleva a pensar
# en un modelo 
# ARIMA(1,1,0)
# si pensamos en un modelo arima con un parametro autoregresivo
# el modelo con una diferencia es lo mismo que modelar a zt
mod.0 <- sarima(bush$approve,1,1,0)
mod.0
#
# Vemos que el par�metro autoregresivo no es significativo por lo
# por lo que el modelo que implica es # Z(t) = Z(t-1)+a(t), consiguiente,
# podr�amos esperar que un modelo ARIMA(1,0,0) pudiera ser razonable
# e intuitivamente el par�metro phi deber�a ser cercano a 1, por 
# lo que ahora consideramos la estimaci�n de un modelo ARIMA(1,0,0)
# Estimamos modelo arima (1,0,0)

mod.1 <- sarima(bush$approve,1,0,0)
mod.1

# Procedemos ahora a estimar el modelo de intervenci�n para septiembre con 
# la serie integrada
#
# Analicemos la gr�fica para estudiar la forma de la intervenci�n
#
# en la lamina 15 de la presentcion, grafica que sube y baja row 2, column 1

tsplot(x=fecha,y=bush$approve)
abline(v=fecha[9], col="red")

# Claramente el tipo de intervenci�n se puede representar 
# con la ecuaci�n:# (1-delta B)epsilon(I,t)=w P(I,t) con I=9 
# (ver gr�fica (c) en la l�mina 15)
# por lo que epsilon(I,t)=w*delta^(t-I), t >= I, 0 de otra forma

mod.2 <- arimax(bush$approve, order=c(1,0,0), xtransf=bush$s11, 
                transfer=list(c(1,0)), include.mean = TRUE); 
mod.2
#t1 ar1 delta 
#t1 ma0 omega cero

names(mod.2)
CuatroEnUno(mod.2$residuals)

#los residuales se comportan razonablemnte bien, no es sorpresivo porque
#tenemos la intervenci�n

tsplot(x=fecha, y=bush$approve)
lines(x=fecha, y=fitted(mod.2), col="red",type="l") #modelo con arimax, peque�o desfase

# Nota: El segundo par�metro corresponde al numerador. Si es 0, 
# entonces es solo un efecto concurrente. El primer par�metro 
# corresponde al denominador (delta).
# Como vemos, este modelo nos da delta estimada de 0.8984 y 
# w estimada de 27.660, ambos estad�sticamente significativos

# Ahora graficamos el modelo de intervenci�n

# tf = epsilon(I,t) = delta*epsilon(I,t-1)+w0P(I,t)= w*delta^(t-I) t>=I

tf <-stats::filter(1*(seq(1:(length(bush$approve)))==09),
                  filter=0.8984,
                  method='recursive',side=1)*(27.666)
plot(tf)

forecast.arima<-Arima(bush$approve,order=c(1,0,0), 
                      seasonal =list(order=c(0,0,0),period=12),
                      xreg=tf) 
forecast.arima

## cambia porque no se esta ajustando la funcion de transferencia

plot(bush$approve, type="l")
lines(forecast.arima$fitted, col="red")

# Para la predicci�n

#se utiliza la longitud del vector de datos y se agregar 6 periodos hacia adelante

tf.ext = stats::filter(1*(seq(1:(length(bush$approve)+6))==09),
                                filter=0.8984,
                                method='recursive',side=1)*(27.666)
predict(forecast.arima,n.ahead = 6, 
        newxreg=tf.ext[length(bush$approve):length(bush$approve)+6])



