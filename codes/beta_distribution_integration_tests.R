### tests of integral function
a <- 1
#mu <- mu_i[test_idx]
mu <- 0.25
#rho <- rho_i
rho <- 10

### draw
par(mfrow = c(1,1))
grid <- seq(-a,1,length=200)
plot(grid, dbeta((a+grid)/(a+1),mu*rho,(1-mu)*rho), type = "l")

### compute expectation by function for (Y = max(0,V))

# probability to be under 0
pbeta(a/(a+1),mu*rho,(1-mu)*rho)

calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu*rho,(1-mu)*rho)/(a+1))
calc_density_with_scaling <- function(x) (dbeta((x+a)/(a+1),mu*rho,(1-mu)*rho)/(a+1))
calc_density_without_scaling <- function(x) (dbeta((x+a)/(a+1),mu*rho,(1-mu)*rho))

integrate(calc_density,0,1)$value
integrate(calc_density_with_scaling,0,1)$value
integrate(calc_density_with_scaling,-a,1)$value
integrate(calc_density_without_scaling,-a,1)$value


library(cubature)
adaptIntegrate(calc_density,0,1)

### compute by sampling
sample_01 <- rbeta(100,mu*rho,(1-mu)*rho)
sample_a1 <- (a+1)*sample_01 - a
mean(sapply(sample_a1, FUN = function(x) (max(0,x))))


### problem with Chara beta model :-)
a.sam <- as.matrix(mod_chara.beta, pars = "alpha")
b.sam <- as.matrix(mod_chara.beta, pars = c("beta_1","beta_2"))
rho.sam <- as.matrix(mod_chara.beta, pars = "rho")

pred_mat <- pred_grid_2021_july_df[,2:10]
pred_mat <- scale_covariates(X,pred_mat)
pred_mat <- add_second_order_terms(pred_mat, colnames(pred_mat))
pred_mat <- as.matrix(pred_mat)

for(i in 1:1) {
  alpha_i <- a.sam[i,]
  beta_i <- b.sam[i,]
  rho_i <- rho.sam[i,]
  
  mu_i <- inv_logit(as.vector(alpha_i + pred_mat %*% beta_i))
  test_idx <- 1
  for(j in 1:length(mu_i)) {
    calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu_i[j]*rho_i,(1-mu_i[j])*rho_i)/(a+1))
    integrate(calc_density,0,1)$value
    test_idx <- test_idx + 1
  }
}

grid <- seq(-a,1,length=200)
plot(grid, dbeta((a+grid)/(a+1),mu_i[test_idx]*rho_i,(1-mu_i[test_idx])*rho_i), type = "l")

calc_density <- function(x) (x*dbeta((x+a)/(a+1),mu_i[test_idx]*rho_i,(1-mu_i[test_idx])*rho_i)/(a+1))
integrate(calc_density,0,1)$value
