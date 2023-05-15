rm(list=ls())
library(coda)

#############################################
## Useful functions.
#############################################
rigamma <- function (n, shape, scale = 1){
    return(1/rgamma(n = n, shape = shape, rate = scale))
}

rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
        stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

#############################################
## Generate some data.
#############################################
set.seed(1)

n <- 100
X <- cbind(1, rnorm(n))
p <- ncol(X)

sigma.sq <- 1
beta <- c(-1,10)

## Make some random covariance matrix just for fun.
A <- matrix(rnorm(n^2), n, n)
V.y <- A%*%t(A) 

## Draw the y.
y <- rmvn(1, X%*%beta, sigma.sq*V.y)

#############################################
## Priors.
#############################################
## sigma^2 ~ IG(a, b)
a <- 2
b <- 1

## beta ~ MVN(mu_beta, V_beta)
mu.beta <- rep(0, p)
V.beta <- diag(1000, p)

#############################################
## Fixed quantities and sample stuff
#############################################
V.y.inv <- chol2inv(chol(V.y))
V.beta.inv <- chol2inv(chol(V.beta))
m <- V.beta.inv%*%mu.beta + t(X)%*%V.y.inv%*%y
M <- chol2inv(chol(V.beta.inv + t(X)%*%V.y.inv%*%X))

## Number of exact samples.
N <- 1000
sigma.sq.samps <- rep(0, N)
beta.samps <- matrix(0, N, p)

#############################################
## Draw samples via composition sampling.
#############################################
for(s in 1:N){

    ## Draw sigma^2
    a.str <- a + n/2
    b.str <- b + 1/2*(t(mu.beta)%*%V.beta.inv%*%mu.beta + t(y)%*%V.y.inv%*%y - t(m)%*%M%*%m)
    sigma.sq.s <- rigamma(1, a.str, b.str)

    ## Draw beta
    beta.s <- rmvn(1, M%*%m, sigma.sq.s*M)

    ## Save samples.
    sigma.sq.samps[s] <- sigma.sq.s
    beta.samps[s,] <- beta.s
}

## Take a look.
samps <- cbind(sigma.sq.samps, beta.samps)
colnames(samps) <- c("sigma^2", paste0("beta_",1:p))
plot(mcmc(samps))

## See Banerjee (2021) https://arxiv.org/pdf/2109.04447.pdf for more details.
