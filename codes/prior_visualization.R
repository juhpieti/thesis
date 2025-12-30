### SCRIPT TO VISUALLY EXAMINE PRIOR CHOICES FOR THE MODEL ###

# pdf of inverse gamma
dinv_gamma <- function(y,a,b) {
  return( ((b^a) / gamma(a)) * y^(-(a+1))*exp(-b/y)  )
}

# cdf of inferse gamma
pinvgamma <- function(x, shape, scale) {
  pgamma(scale / x, shape, lower.tail = FALSE)  # Complement of upper incomplete gamma function
}

# pdf of student-t 
dstudent <- function(y,v,u,s) {
  upper <- gamma((v+1)/2)
  lower <- (sqrt(pi*v*s^2))*gamma(v/2)
  const <- upper/lower
  c <- -(v+1)/2
  return(const*((1+((y-u)^2)/(v*s^2))^c))
}

### 1) priors for beta coefficients
### N(0,2)

par(mfrow = c(1,2))
n_rep <- 1000
beta_1 <- rnorm(n_rep,0,sqrt(2))
hist(beta_1, breaks = 40, probability = TRUE)

### 2) covariance function parameters

### S2 (magnitude)
s2_grid <- seq(0.01,10.01,length=500)


### inv_gamma(a,b) vs half-cauchy(0,t) for s2?
### e.g. inv_gamma(2,1) vs half-cauchy(0,1)
plot(NULL, xlim = c(0,10), ylim = c(0,0.8))
lines(s2_grid, 2*dcauchy(s2_grid,0,sqrt(0.001)),col="red")
lines(s2_grid,2*dstudent(s2_grid,2,0,sqrt(0.1)), col = "orange")
lines(s2_grid,2*dstudent(s2_grid,4,0,sqrt(0.1)), col = "purple")
lines(s2_grid,2*dstudent(s2_grid,2,0,1.5), col = "blue")
lines(s2_grid,2*dstudent(s2_grid,4,0,1.5), col =  "darkgreen")
# draw jarnos student-t
lines(s2_grid, 2*dstudent(sqrt(s2_grid),4,0,0.1)/(2*sqrt(s2_grid)), col = "green")

### probability that s2 is under 10
pinvgamma(5,2,1) #0.99
pcauchy(5,0,sqrt(0.001)) - pcauchy(-5,0,sqrt(0.001)) #0.93

### lengthscale?
### the minimum distance is 20km between grid points, and the maximum is 580km
### if l = 865, covariance at distance 600 is 50% of the one at distance 0.
### let's say that we expect that at 80km 5% of the one at distance 0: l \approx 25
### the "realistic" value for l could be e.g. between (0,1000)

### let's visualize jarnos prior for gulf of finland thing
l_grid <- seq(0,400,length=500)
plot(NULL, xlim = c(min(l_grid),max(l_grid)), ylim = c(0,0.02))
lines(l_grid, 2*dstudent(1/l_grid,4,0,1)/l_grid^2, col = "blue")
lines(l_grid,2*dstudent(l_grid,4,10,10), col = "orange")
lines(l_grid,2*dcauchy(l_grid,8,40), col = "red")
lines(l_grid,2*dstudent(l_grid,4,8,40), col = "red")
lines(l_grid,2*dcauchy(l_grid,8,sqrt(50)), col = "purple")
lines(l_grid,dinv_gamma(l_grid,1,1), col = "green")

abline(v = 40)

# calculate probabilities for length-scale being less than ?
pinvgamma(1000,1,1) #0.99
pcauchy(100,8,25) - pcauchy(-100,8,25) #0.95

# calculate expectations
integrate(function(x) (x*2*dstudent(x,4,10,10)),0,500)
integrate(function(x) x*2*dcauchy(x,8,50),0,500)

#### RHO
### since the variance of beta distribution depends inversely on rho, I think it is logical to give it a prior with heavy tails

### half-cauchy(0,sqrt(5)) vs inv_gamma(0.5, 0.5) for rho?
rho_grid <- seq(0,100,length=500)
plot(NULL, xlim = c(min(rho_grid),max(rho_grid)), ylim = c(0,0.2))
lines(rho_grid, 2*dcauchy(rho_grid,0,sqrt(10)), col = "green")
lines(rho_grid, 2*dcauchy(rho_grid,0,5),col="red")
lines(rho_grid, 2*dnorm(rho_grid,0,sqrt(10)), col = "blue")

### probability that rho is under 100
pinvgamma(100,0.1,0.1) #0.92
pcauchy(100,0,sqrt(10)) - pcauchy(-100,0,sqrt(10)) #0.98

### 1/rho
### it could also make sense to give prior for 1/rho 
inv_rho_grid <- seq(0.01,2.01,length=500)

### student-t(df=4,0,?)
plot(NULL, xlim = c(0,2), ylim = c(0,2))
lines(inv_rho_grid,2*dstudent(inv_rho_grid,2,0,0.1), col = "orange")
lines(inv_rho_grid,2*dstudent(inv_rho_grid,4,0,0.1), col = "purple")
abline(v=1, col = "red", lty = 2)

################################################### EXTRA ###########################################
### draw plot of joint prior for l and s2

l_grid <- seq(0,500,length=100)
s2_grid <- seq(0,5,length=100)

### let's put half-cauchy(0,0.01) and half-cauchy(0,50)
ls2_grid <- expand.grid(l = l_grid, s2 = s2_grid)

log_lik <- function(l,s2) {
  return(2*log(2)+dcauchy(l,8,sqrt(100),log=TRUE)+dcauchy(s2,0,sqrt(0.01),log=TRUE))
}

ll <- mapply(log_lik, ls2_grid$l, ls2_grid$s2)

params_ll <- cbind(ls2_grid,ll)

ggplot(params_ll, aes(x = l, y = s2, fill = ll)) +
  geom_tile() +  # Alternative: use geom_raster() for larger datasets
  scale_fill_viridis_c() +  # Better color scale for log-likelihood
  labs(x = "Parameter 1", y = "Parameter 2", fill = "Log-Likelihood") +
  theme_minimal()

### normal vs laplace

#install.packages("VGAM")
library(VGAM)

x_grid <- seq(-5,5,length=200)
plot(x_grid, dnorm(x_grid,0,sqrt(2)), type = "l", col = "blue", ylim = c(0,0.5))
lines(x_grid, dlaplace(x_grid,0,1), col = "red")

var(rnorm(1000,0,sqrt(2)))
var(rlaplace(10000,0,1))
