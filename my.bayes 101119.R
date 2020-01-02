# Generating Data ---------------------------------------------------------
set.seed(1234)
library(tidyverse)

# B0 = 15, B1 = .3, B2 = .6, sd = .74
mydata<-data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))%>%
  mutate(y = 15 +.3 * x1+.6 * x2 + .15*x3+ rnorm(100, mean = 0, sd = sqrt(1-(.6^2+.3^2+.15^2))))
lmfit  <- lm(y ~ x1 + x2 + x3, mydata)

# sourcing private functions ----------------------------------------------
source(file = "/Users/galensnyder/Documents/R scripts/Nelder.Mead.bounds.R", local = T)
source(file = "/Users/galensnyder/Documents/R scripts/flexlm.R", local = T)

dinvgamma <- function(x, shape, scale, log = F){
  a <- shape
  b <- scale
  log.density <- a * log(b) - lgamma(a) - (a + 1) * log(x) - (b/x)
  if(log == T){
    return(log.density)
  } else {
    return(exp(log.density))
  }
}

autocorrplot <- function(theta, lag, ...){
  num <- vector(mode = "numeric", length = length(lag))
  
  for(i in 1:length(lag)){
    num[i] <- sum((theta[1:(length(theta)-lag[i])]-mean(theta)) * (theta[(1+lag[i]):length(theta)]-mean(theta)))
  }
  
  denom <- sum((theta - mean(theta))^2)
  
  auto <- num / denom
  
  plot(x = lag, y = auto, type = "h", xlab = "Lag", ylab = "Autocorrelation", ...)
}

my.bayes <- function(formula, dat, n = 12500, chains = 4, burn = 2500, adapt = F, acc.star = 0.234, plots = FALSE, x = NULL, y = NULL){
  cols <- c("black", "red", "blue", "green", "orange", "yellow", "violet")
  proposal.sds <- c(100, 100, 100)
  mod <- lm.model.mod(formula, dat)
  p <- mod$p
  thetastart <- runif(p, .01, .99)
  fit <- nelder.mead(thetastart, flexlm.fit, hessian = T, tol = .Machine$double.eps, mod = mod)
  
  x  <- mod$x
  y  <- mod$y
  N  <- mod$N
  p <- mod$p
  df <- mod$N - mod$p
  b  <- fit$parm
  vb <- solve(0.5*fit$hessian)
  rss<- sum((mod$y - mod$x%*%b)^2)
  s2 <- rss / df
  
  theta <- c(b, s2)
  theta <- cbind(theta, matrix(0, nrow = p+1, ncol = n-1))
  # prop.sigma <- c(sqrt(diag(vb)), 1)
  
  old.priors  <- sum(dnorm(theta[1:p, 1], 0, sqrt(1e6), log = T)) + dinvgamma(theta[p+1, 1], .01, .01, log = T)
  old.lik     <- sum(dnorm(y, x%*%theta[1:p, 1], sqrt(theta[p+1, 1]), log = T))
  old.post    <- old.priors + old.lik
  
  rho <- rep(sqrt(2.38), p+1)
  old.rho <- rho
  S <- lapply(1:(p+1), function(k) (rho[[k]]^2)*diag(1, 1))
  alpha <- 0.75
  ar <- rep(.44, p+1)
  ar.hat <- rep(0, p+1)
  b0 <- 0.8
  gam <- 0
  adapt.iter <- rep(1L, p+1)
  last.adapt <- rep(burn, p+1)
  acc <- lapply(1:(p+1), function (k) vector(mode = "numeric", length = n))
  accept <- lapply(1:(p+1), function(i) vector(mode = "numeric", length = n))
  # prop <- theta[,1]
  for(i in 2:n){
    prop      <- theta[, i-1]
    for(j in 1:(p+1)){
      prop[j] <- prop[j] + rnorm(1, 0, rho[j])
      if(j == p+1){
        while(prop[j] < 0){
          prop[j] <- theta[j, i-1] + rnorm(1, 0, rho[j])
        }
        priors  <- sum(dnorm(theta[1:p, i-1], 0, sqrt(1e6), log = T)) + dinvgamma(prop[j], .01, .01, log = T)
        lik     <- sum(dnorm(y, x%*%theta[1:p, i-1], sqrt(prop[j]), log = T))
      } else {
        priors  <- dnorm(prop[j], 0, sqrt(1e6), log = T) + sum(dnorm(theta[-c(j,p+1), i-1], 0, sqrt(1e6), log = T)) + dinvgamma(theta[p+1, i-1], .01, .01, log = T)
        lik     <- sum(dnorm(y, prop[j]*x[,j] + as.matrix(x[, -j])%*%prop[-c(j, p+1)], sqrt(theta[p+1, i-1]), log = T))
      }
      # priors  <- sum(dnorm(prop[1:p], 0, sqrt(1e6), log = T)) + dinvgamma(prop[p+1], .01, .01, log = T)
      # lik     <- sum(dnorm(y, x%*%prop[1:p], sqrt(prop[p+1]), log = T))
      post    <- priors + lik
      acc.r   <- exp(post - old.post)
      u       <- runif(1)
      if(u > min(1, acc.r)){
        prop[j] <- theta[j, i-1]
        acc[[j]][i] <- u
        # ar[j]      <- u
      } else {
        accept[[j]][i] <- 1
        acc[[j]][i] <- acc.r
        # ar[j]        <- acc.r
      }
      
      # prop[p+1] <- prop[p+1] + rnorm(1, 0, 1)
      # while(prop[p+1] < 0){
      #   prop[p+1] <- theta[p+1, i-1] + rnorm(1, 0, 1)
      # }
      # priors  <- sum(dnorm(prop[1:p], 0, 1e6, log = T)) + dinvgamma(prop[p+1], .01, .01, log = T)
      # lik     <- sum(dnorm(y, x%*%prop[1:p], sqrt(prop[p+1]), log = T))
      # post    <- priors + lik
      # acc.r   <- exp(post - old.post)
      # u       <- runif(1)
      # if(u > min(1, acc.r)){
      #   prop[p+1] <- theta[p+1, i-1]
      #   # acc <- u
      # } else {
      #   old.post <- post
      # }
      if(adapt == T && i > burn && i%%20 == 0 && mean(accept[[j]][last.adapt[j]:i]) < .25 ){
        ar.hat[j] <- mean(acc[[j]][last.adapt[j]:i])
        ar[j] <- (1 - alpha)*ar[j] + alpha*ar.hat[j]
        # old.rho[j] <- rho[j]
        rho[j] <- rho[j]*exp((b0/adapt.iter[j]^gam) * (1/pnorm(ar[j]/2) - 1/pnorm(.22)))
        # S[[j]] <- (rho[j]^2)*((1-(b0/adapt.iter[j]^gam))*S[[j]] + ((b0/adapt.iter[j]^gam)*S[[j]]))
        adapt.iter[j] <- adapt.iter[j] + 1
        last.adapt[j] <- i
      }
    }
    old.priors  <- sum(dnorm(prop[1:p], 0, sqrt(1e6), log = T)) + dinvgamma(prop[p+1], .01, .01, log = T)
    old.lik     <- sum(dnorm(y, x%*%prop[1:p], sqrt(prop[p+1]), log = T))
    old.post    <- old.priors + old.lik
    theta[, i]  <- prop
  }
  
  theta.out <- theta[, -(1:burn)]
  theta.names <- c(mod$namesFE, "residual")
  out <- data.frame(theta.out, row.names = theta.names)
  out <- apply(out, 1, function(i){
    c("mean" = mean(i),
      "std. err" = sd(i),
      "t-value" = mean(i)/sd(i),
      "ci lower" = quantile(i, 0.025),
      "ci upper" = quantile(i, 0.975)
    )
  })
  print(t(out))
  
  varcov <- matrix(0, p+1, p+1)
  for(i in 1:(p+1)){
    for(j in 1:(p+1)){
      varcov[i,j] <- var(theta.out[i, ], theta.out[j, ])
    }
  }
  
  print(varcov)
  
  if(plots){
    for(i in 1:(p+1)){
      autocorrplot(theta[i, burn:n], seq(1:100), 
                   main = paste("Autocorrelation of ", theta.names[i], sep = ""))
      
      plot(theta[i, burn:n], type = "l", main = paste("Traceplot of ", theta.names[i], sep = ""),
           xlab = "Iteration", ylab = "Value")
      
      abline(h = mean(theta[i, burn:n]), col = "red")
    }
    message("Use 'previous plot' option in plots & files window to see all plots")
  } 
  invisible(theta)
}

out <- my.bayes(y~x1 + x2 + x3, mydata, n = 12500, adapt = T)
out <- my.bayes(y~x1 + x2 + x3, mydata, n = 12500, adapt = F)





