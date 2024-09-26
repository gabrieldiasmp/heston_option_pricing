library(NMOF)
library(dplyr)
library(DEoptim)
library(rio)
library(ggplot2)

iv_surface = rio::import("data/iv_surface.csv") %>% rename("t_ano" = t_dias) %>% filter(t_ano > 0.01) %>% select(spot, strike, t_ano, iv)
iv_surface = rio::import("data/ivolatility.csv") %>% select(stock_price_for_iv, strike, period, iv) %>% mutate(period = period/360)
iv_surface = rio::import("data/spy_options.csv") %>% select(spot, Strike, `Last Price`, `Implied Volatility`, tipo, t_ano) %>% rename("strike" = Strike, "iv" = `Implied Volatility`, "price" = `Last Price`) %>% filter(tipo == "calls", strike > spot) %>% mutate(moneyness = strike/spot) %>% na.omit() #%>% mutate(t_ano = 0.5, moneyness = (spot/strike)-1) #%>% filter(tipo == "calls", strike > spot)
colnames(iv_surface) <- c("spot", "strike", "t_ano", "iv")

ggplot(iv_surface, aes(moneyness, iv)) + geom_point()

# SABR ------------------------------------------------------------

sabrVol <- function(f,k,t,a,b,r,v){
  if(abs(f/k-1) < 1e-04){
    term1 <- a/(f^(1-b))
    term2 <- ((((1-b)^2)/24)*((a^2)/f^(2-2*b)) + r*b*a*v/(4*f^(1-b))
              + (2-3*r^2)*v^2/24)
    y <- term1*(1 + term2*t)
    
  } else{
    fk <- f*k
    z <- v/a*(fk)^((1-b)/2)*log(f/k)
    x <- suppressWarnings(log((sqrt(1-2*r*z + z^2) + z-r)/(1-r)))
    term1 <- a / fk^((1-b)/2) / (1 + (1-b)^2/24*log(f/k)^2 +
                                   (1-b)^4/1920*log(f/k)^4)
    if(is.nan(x)){
      term2 <- 1
    } else{
      term2 <- z / x
    }
    term3 <- 1 + ((1-b)^2/24*a^2/fk^(1-b) +
                    r*b*v*a/(4*fk^((1-b)/2)) + (2-3*r^2)/24*v^2)*t
    y <- term1*term2*term3
  }
  return(y)
}

#df_surface = iv_surface %>% filter(tipo == "C")
#df_surface = iv_surface %>% filter(tipo %in% c("A", "C"))
df_surface = iv_surface

SABR.calibration <- function(t = df_surface$t_ano, f = df_surface$spot, k = df_surface$strike, iv = df_surface$iv){
  
  SABR.prediction <- function(p){
    # objective function for optimization
    # variables are transformed because of satisfing the constraint conditions
    objective <- lapply(1:nrow(df_surface), function(linha) {
      sabr_prediction <- sabrVol(t = t[linha], f = f[linha], k = k[linha], a = p[1], b = p[2], r = p[3], v = p[4])
      
      soma_erro_quadrado <- (iv[linha] - sabr_prediction)^2
      
      soma_erro_quadrado
      
    })
    
    return(do.call(sum, objective))
    
  }
  
  # x <- nlm(SABR.prediction, c(1, 0.5, -0.5, 1))
  # parameter <- x$estimate
  
  parametros_sabr <- DEoptim::DEoptim(SABR.prediction, lower = c(0, 0, -1, 0), upper = c(100, 1, 1, 100), DEoptim.control(itermax = 300))
  parameter <- parametros_sabr$optim$bestmem
  
  
  parameter <- c(parameter[1], parameter[2], parameter[3],parameter[4])
  names(parameter) <- c("Alpha", "Beta", "Rho", "Nu")
  parameter
  
}

bla <- SABR.calibration()
bla

a <- lapply(1:nrow(df_surface), function(n){
  iv_predict <- sabrVol(f = df_surface$spot[n], k = df_surface$strike[n], t = df_surface$t_ano[n], a = bla[1], b = bla[2], r = bla[3], v = bla[4])
  
  iv_predict
  
  # print(c(df_surface$iv[n], iv_predict))
  # 
  # soma_erro_quadrado <- df_surface$iv[n] - iv_predict
  # 
  # print(soma_erro_quadrado)
  
})

a <- do.call(rbind, a)
df_surface <- cbind(df_surface, iv_predict = a)
df_surface$erro <- df_surface$iv - df_surface$Alpha
sum(df_surface$erro^2)
df_surface$moneyness <- (df_surface$strike / df_surface$spot) - 1
with(df_surface, plot(moneyness, Alpha))


# FUNCAO HESTON OTIMIZAVEL ------------------------------------------------

HestonCallClosedForm <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K) {
  PIntegrand <- function(u, lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
    FF <- S0*exp(r*tau)
    x <- log(FF/K)
    a <- lambda * vbar
    
    if (j == 1) {
      b <- lambda - rho* eta
      alpha <- - u^2/2 - u/2 * 1i + 1i * u
      beta <- lambda - rho * eta - rho * eta * 1i * u
    } else {
      b <- lambda
      alpha <- - u^2/2 - u/2 * 1i
      beta <- lambda - rho * eta * 1i * u
    }
    
    gamma <- eta^2/2
    d <- sqrt(beta^2 - 4*alpha*gamma)
    rplus <- (beta + d)/(2*gamma)
    rminus <- (beta - d)/(2*gamma)
    g <- rminus / rplus
    
    D <- rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
    C <- lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
    
    top <- exp(C*vbar + D*v0 + 1i*u*x)
    bottom <- (1i * u)
    Re(top/bottom)
  }
  
  P <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
    value <- integrate(PIntegrand, lower = 0, upper = Inf,
                       lambda, vbar, eta, rho, v0, r, tau,
                       S0, K, j, subdivisions=2000)$value
    0.5 + 1/pi * value
  }
  
  A <- S0*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 1)
  B <- K*exp(-r*tau)*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 0)
  A-B
}

funcao_heston <- function(p){
    
    ## PARAMETERS
    ##
    ## lambda: mean-reversion speed
    ## vbar: long-term average volatility
    ## eta: volatility of vol process
    ## rho: correlation between stock and vol
    ## v0: initial volatility
    ## r: risk-free interest rate
    ## tau: time to maturity
    ## S0: initial share price
    ## K: strike price
  
    df = iv_surface
    df$r = 0.017
    
    erro_quadrado <- lapply(1:nrow(df), function(n){
  
      heston_prediction <- HestonCallClosedForm(K = df$strike[n],
                                                S0 = df$spot[n],
                                                tau = df$t_ano[n],
                                                r = 0.01,
                                                rho = p[1],
                                                v0 = p[2],
                                                vbar = p[3],
                                                eta = p[4],
                                                lambda = p[5])

      # heston_prediction <- NMOF::callHestoncf(S = df$spot[n],
      #                                         X = df$strike[n],
      #                                         tau = df$t_ano[n],
      #                                         r = df$r[n],
      #                                         q = 0,
      #                                         rho = p[1],
      #                                         v0 = p[2],
      #                                         vT = p[3],
      #                                         sigma = p[4],
      #                                         k = p[5],
      #                                         implVol = FALSE)
  
  
      #FUNCAO CUSTO SOMA DOS QUADRADOS DOS ERROS
      soma_erro_quadrado <- (df$price[n] - heston_prediction)^2
      #soma_erro_quadrado <- (df$iv[n] - heston_prediction)^2
      # soma_erro_quadrado <- (df$iv[n] - heston_prediction$impliedVol)^2
      
    })
    
    return(do.call(sum, erro_quadrado))
    
  }

fit_values <- function(df){
  
  a <- lapply(1:nrow(df), function(n){
    
    df$r = 0.017
    
    iv_predict <- NMOF::callHestoncf(S = df$spot[n],
                                       X = df$strike[n],
                                       tau = df$t_ano[n],
                                       r = df$r[n],
                                       q = 0,
                                       #implVol = TRUE,
                                       rho = parameters[1],
                                       v0 = parameters[2],
                                       vT = parameters[3],
                                       sigma = parameters[4],
                                       k = parameters[5])
    iv_predict = iv_predict
    print(iv_predict)
    
  })
  
  a <- do.call(rbind, a)
  a <- cbind(df, a)
  
  return(a)
}


parametros_heston <- DEoptim::DEoptim(funcao_heston, lower = c(-1, 0.01, 0.01, 0.01, 0.01), upper = c(0, 0.5, 0.5, 1, 1), DEoptim.control(itermax = 50))

parametros_heston <- DEoptim::DEoptim(funcao_heston, lower = c(-1, 0.01, 0.01, 0.01, 0.01), upper = c(0, 1, 1, 1, 1), DEoptim.control(itermax = 30))
tst <- NMOF::gridSearch(funcao_heston, lower = c(-1, 0.01, 0.01, 0.001, 0.001), upper = c(0, 5, 5, 0.9^2, 0.9^2))
tst <- 
#parametros_heston <- NMOF::DEopt(funcao_heston, algo = list(min = c(-1, 0.001, 0.001, 0.001, 0.001), max = c(0, 1, 1, 0.9^2, 0.9^2), minmaxConstr = TRUE, nG = 30))
parameters <- parametros_heston$optim$bestmem
names(parameters) <- c("rho", "v0", "vT", "sigma", "k")
parameters

iv_surface <- fit_values(iv_surface)
summary(lm(price ~ a, data = iv_surface))
with(iv_surface, plot(a, price))

with(iv_surface, plot(moneyness, a))
with(iv_surface, plot(a, iv))


# EXEMPLO COM PARAMETROS A MAO --------------------------------------------

df = iv_surface #%>% filter(tipo == "C")

fit_values <- function(df){
  
  a <- lapply(1:nrow(df), function(n){
    
    iv_predict <- NMOF::callHestoncf(S = df$spot[n],
                                     X = df$strike[n],
                                     tau = df$t_ano[n],
                                     r = 0.017,
                                     q = 0,
                                     rho = -0.831705,
                                     v0 = 0.031475,
                                     vT = 0.058971,
                                     sigma = 0.528437,
                                     k = 0.072192,
                                     implVol = FALSE)
    
  })
  
  a <- do.call(rbind, a)
  a <- cbind(df, a)
  
  return(a)
}

tst <- fit_values(df)

teste <- lapply(1:nrow(df), function(n){
  
  # heston_prediction <- NMOF::callHestoncf(S = df$spot[n],
  #                    X = df$strike[n],
  #                    tau = df$t_ano[n],
  #                    r = df$r[n],
  #                    q = 0,
  #                    v0 = 0.03,
  #                    vT = 0.270041,
  #                    rho = -0.623425,
  #                    k = 0.874219,
  #                    sigma = 0.390626,
  #                    implVol = TRUE)
  
  heston_prediction <- HestonCallClosedForm(K = df$strike[n],
                                            S0 = df$spot[n],
                                            tau = df$t_ano[n],
                                            r = 0.01,
                                            rho = -0.854380,
                                            v0 = 0.011831,
                                            vbar = 0.570433,
                                            eta = 0.581208,
                                            lambda = 0.023547)
  
  
  #erro_quadrado = (heston_prediction$impliedVol - df$iv[n])^2
  erro_quadrado = (heston_prediction - df$iv[n])^2
  
  
  
})

do.call(sum, teste)



# OUTRAS MANEIRAS DE CÁLCULO HESTON (FUNCIONANDO) ---------------iv_surface------------------------

HestonCallClosedForm <-
  function(lambda, vbar, eta, rho, v0, r, tau, S0, K) {
    PIntegrand <- function(u, lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
      F <- S0*exp(r*tau)
      x <- log(F/K)
      a <- lambda * vbar
      
      if (j == 1) {
        b <- lambda - rho* eta
        alpha <- - u^2/2 - u/2 * 1i + 1i * u
        beta <- lambda - rho * eta - rho * eta * 1i * u
      } else {
        b <- lambda
        alpha <- - u^2/2 - u/2 * 1i
        beta <- lambda - rho * eta * 1i * u
      }
      
      gamma <- eta^2/2
      d <- sqrt(beta^2 - 4*alpha*gamma)
      rplus <- (beta + d)/(2*gamma)
      rminus <- (beta - d)/(2*gamma)
      g <- rminus / rplus
      
      D <- rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
      C <- lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
      
      top <- exp(C*vbar + D*v0 + 1i*u*x)
      bottom <- (1i * u)
      Re(top/bottom)
    }
    
    P <- function(lambda, vbar, eta, rho, v0, r, tau, S0, K, j) {
      value <- integrate(PIntegrand, lower = 0, upper = Inf,
                         lambda, vbar, eta, rho, v0, r, tau,
                         S0, K, j, subdivisions=1000)$value
      0.5 + 1/pi * value
    }
    
    A <- S0*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 1)
    B <- K*exp(-r*tau)*P(lambda, vbar, eta, rho, v0, r, tau, S0, K, 0)
    A-B
  }

funcao_heston <- function(p){
  
  #print(list(p[1], p[2], p[3], p[4], p[5]))
  
  soma_erro_quadrado = list()
  df = iv_surface %>% filter(tipo == "C")
  
  erro_quadrado <- lapply(1:nrow(df), function(n) {
    heston_prediction <-
      HestonCallClosedForm(K = df$strike[n],
                           S0 = df$spot[n],
                           tau = df$t_ano[n],
                           r = df$r[n],
                           rho = p[1],
                           v0 = p[2],
                           vbar = p[3],
                           eta = p[4],
                           lambda = p[5]
      )
      
    
    soma_erro_quadrado <- (heston_prediction - df$iv[n])^2
    soma_erro_quadrado
    
    # tst <- c(df$iv[n], heston_prediction)
    # print(tst)
  })
  
  return(do.call(sum, erro_quadrado))
  
}

parametros_heston <- DEoptim::DEoptim(funcao_heston, lower = c(-1, 0.01, 0.01, 0.05^2, 0.05^2), upper = c(1, 5, 5, 0.9^2, 0.9^2))
parameters <- parametros_heston$optim$bestmem

df_surface = iv_surface %>% filter(tipo == "C")

a <- lapply(1:nrow(df_surface), function(n){
  
  iv_predict <- HestonCallClosedForm(K = df_surface$strike[n],
                                           S0 = df_surface$spot[n],
                                           tau = df_surface$t_ano[n],
                                           r = df_surface$r[n],
                                           rho = parameters[1],
                                           v0 = parameters[2],
                                           vbar = parameters[3],
                                           eta = parameters[4],
                                           lambda = parameters[5]
  )
  
})

a <- do.call(rbind, a)
df_surface <- cbind(df_surface, a)



# n = 7
# HestonCallClosedForm(K = iv_surface$strike[n],
#                      S0 = iv_surface$spot[n],
#                      tau = iv_surface$t_ano[n],
#                      r = iv_surface$r[n],
#                      rho = -0.75,
#                      v0 = 0.2,
#                      vbar = 0.2,
#                      eta = 0.01,
#                      lambda = 0.01
# )
# 
# iv_surface$iv[7]

# MODELO SABR V2 -------------------------------------------------------------

EPS <- 10^(-8)

# sub function for SABR BS-IV
.x <- function(z, r){log((sqrt(1-2*r*z+z^2)+z-r)/(1-r))}
.z <- function(f, K, a, b, nu){nu/a*(f*K)^(0.5*(1-b))*log(f/K)}

# variable transformation function
.t1  <- function(x){1/(1+exp(x))}
.t2  <- function(x){2/(1+exp(x)) -1}

# Black-Scholes IV apporoximation formula by Hagan(2002)
SABR.BSIV <- function(t, f, K, a, b, r, n)
{
  z <- .z(f, K, a, b, n)
  x <- .x(z, r)
  numerator   <- 1 + ((1-b)^2/24*a^2/(f*K)^(1-b) + 0.25*r*b*n*a/(f*K)^(0.5*(1-b)) + (2-3*r^2)*n^2/24)*t
  denominator <- x*(f*K)^(0.5*(1-b))*(1 + (1-b)^2/24*(log(f/K))^2 + (1-b)^4/1920*(log(f/K))^4)
  ifelse(abs((f-K)/f) < EPS, a*numerator/f^(1-b), z*a*numerator/denominator)
}

#####TESTE SABR COM OTIMIZAÇÃO NEWTON
df_surface = iv_surface %>% filter(tipo == "C")

SABR.calibration <- function(t = df_surface$t_ano, f = df_surface$spot, K = df_surface$strike, iv = df_surface$iv){
  
  SABR.prediction <- function(p){
    # objective function for optimization
    # variables are transformed because of satisfing the constraint conditions
    objective <- lapply(1:nrow(df_surface), function(linha) {
      heston_prediction <- SABR.BSIV(t[linha], f[linha], K[linha], exp(p[1]), .t1(p[2]), .t2(p[3]), exp(p[4]))
      
      soma_erro_quadrado <- (iv[linha] - heston_prediction)^2
      print(soma_erro_quadrado)
      
      soma_erro_quadrado
      
    })
    
    return(do.call(sum, objective))
    
  }
  
  x <- nlm(SABR.prediction, c(0.1, 0.5, 0.0, 0.1))
  # return the optimized parameters
  parameter <- x$estimate
  parameter <- c(exp(parameter[1]), .t1(parameter[2]), .t2(parameter[3]), exp(parameter[4]))
  names(parameter) <- c("Alpha", "Beta", "Rho", "Nu")
  parameter
  
}

bla <- SABR.calibration()


a <- lapply(1:nrow(df_surface), function(n){
  iv_predict <- SABR.BSIV(f = iv_surface$spot[n], K = iv_surface$strike[n], t = iv_surface$t_ano[n], bla[1], bla[2], bla[3], bla[4])
  
  print(c(df_surface$iv[n], iv_predict))
  
  soma_erro_quadrado <- (df_surface$iv[n] - iv_predict)^2
})

a <- do.call(rbind, a)
df_surface <- cbind(df_surface, a)
df_surface[, c("iv", "a")]


#### Parameter calibration function for SABR (ORIGINAL)
SABR.calibration <- function(t, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  objective <- function(x){sum( (iv - SABR.BSIV(t, f, K, exp(x[1]), .t1(x[2]), .t2(x[3]), exp(x[4])))^2) }
  x <- nlm(objective, c(0.1, 0.5, 0.0, 0.1))
  # return the optimized parameters
  parameter <- x$estimate
  parameter <- c(exp(parameter[1]), .t1(parameter[2]), .t2(parameter[3]), exp(parameter[4]))
  names(parameter) <- c("Alpha", "Beta", "Rho", "Nu")
  parameter
}

n = 7
bla <- SABR.calibration(f = iv_surface$spot[n], K = iv_surface$strike[n], t = iv_surface$t_ano[n], iv = iv_surface$iv[n])

for(n in 1:nrow(iv_surface)){
  print(SABR.BSIV(f = iv_surface$spot[n], K = iv_surface$strike[n], t = iv_surface$t_ano[n], a = bla[1], b = bla[2], r = bla[3], n = bla[4]))
}


####FUNCIONAL DEopt

SABR_function <- function(p){ ####CALIBRAÇÃO
  
  soma_erro_quadrado = list()
  df = iv_surface[1:8, ]
  
  erro_quadrado <- lapply(1:nrow(df), function(linhas) {
    heston_prediction <-
      SABR.BSIV(f = df$spot[linhas],
                K = df$strike[linhas],
                t = df$t_ano[linhas],
                a = exp(p[1]),
                b = .t1(p[2]),
                r = .t2(p[3]),
                n = exp(p[4]))
    
    # tst <- c(df$iv[linhas], heston_prediction)
    # print(tst)
    
    soma_erro_quadrado <- (df$iv[linhas] - heston_prediction)^2
    
    soma_erro_quadrado
    
  })
  
  return(do.call(sum, erro_quadrado))
}



DEoptim::DEoptim(SABR_function, lower = c(0, 0, -1, 0), upper = c(1, 1, 1, 1))





# Parameter calibration function for SABR
SABR.calibration <- function(t, f, K, iv){
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  objective <- function(x){sum( (iv - SABR.BSIV(t, f, K, exp(x[1]), .t1(x[2]), .t2(x[3]), exp(x[4])))^2) }
  x <- nlm(objective, c(0.1, 0.5, 0.0, 0.1))
  # return the optimized parameters
  parameter <- x$estimate
  parameter <- c(exp(parameter[1]), .t1(parameter[2]), .t2(parameter[3]), exp(parameter[4]))
  names(parameter) <- c("Alpha", "Beta", "Rho", "Nu")
  parameter
}

SABR.calibration(f = iv_surface$spot[n], K = iv_surface$strike[n], t = iv_surface$t_ano[n], iv = iv_surface$iv[n])

# TESTANDO OUTROS MÉTODOS -------------------------------------------------

SABR.prediction <- function(p){
  t = df_surface$t_ano
  f = df_surface$spot
  k = df_surface$strike
  iv = df_surface$iv
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  objective <- lapply(1:nrow(df_surface), function(linha) {
    sabr_prediction <- sabrVol(t = t[linha], f = f[linha], k = k[linha], a = p[1], b = p[2], r = p[3], v = p[4])

    soma_erro_quadrado <- (iv[linha] - sabr_prediction)^2
    
    soma_erro_quadrado
    
  })
  
  return(do.call(sum, objective))
  
}

#x <- nlm(SABR.prediction, c(0.1, 0.5, 0.0, 0.1))
#parameter <- x$estimate
parametros_sabr <- nlm(SABR.prediction, c(0.5, 0.5, 0.0, 0.5))
parametros_sabr <- parametros_sabr$estimate
parametros_sabr <- DEoptim::DEoptim(SABR.prediction, lower = c(0, 0, -1, 0), upper = c(10000, 1, 1, 10000))
parametros_sabr <- NMOF::gridSearch(SABR.prediction, lower = c(0, 0, -1, 0), upper = c(10000, 1, 1, 10000))
parametros_sabr <- NMOF::GAopt(SABR.prediction, algo = list(nB = 20L))

NMOF::
