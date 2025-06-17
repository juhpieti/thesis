### test how the beta distribution actually works
# the interpretation of sample size parameter:
# the larger value, the smaller variance: how concentrated the distribution is around the mean

mu <- 0.5
sample_size <- c(0.1,1,2,5,10,50)
x_grid <- seq(0,1,length=100)

cols <- rainbow(length(sample_size))

plot(NULL, xlim = 0:1, ylim = c(0,6))

for (i in 1:length(sample_size)) {
  ss <- sample_size[i]
  print(ss)
  lines(x_grid,dbeta(x_grid,mu*ss,(1-mu)*ss), type = "l", col = cols[i])
}

legend("topright", legend = paste0("sample size = ", sample_size), lty = 1, col = cols)


### considering the measurement error (5%), think about how large the rho parameter can be
rho <- c(2,4,10,20,50,100,250,500,1000)
mu <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5)
x_seq <- seq(0,1,length=200)
for (mu_i in mu) {
  par(mfrow = c(3,3))
  for (rho_i in rho) {
    prob_within_5 <- pbeta(mu_i+0.025, mu_i*rho_i, (1-mu_i)*rho_i) - pbeta(mu_i-0.025, mu_i*rho_i, (1-mu_i)*rho_i)
    plot(x_seq,dbeta(x_seq,mu_i*rho_i,(1-mu_i)*rho_i), type = "l", main = paste0("mu=",mu_i,"rho=",rho_i))
    legend("topright",legend = prob_within_5)
    abline(v=0.475,col="red",lty=2)
    abline(v=0.525,col="red",lty=2)
    
  }
}

### examine why rho models produce U-shaped responses...?
mod.beta <- readRDS("models/n_500/M1/Cladophora_glomerata.rds")
mod.beta_rho <- readRDS("models/scaled_sigmoid/n_500/M1/Cladophora_glomerata.rds")

so_val <- 8.25
so_val <- 3.4
so_val_scaled <- (so_val - mean(X$so_bottom))/sd(X$so_bottom)

alpha_base <- colMeans(as.matrix(mod.beta, pars = c("alpha")))
b6_base <- colMeans(as.matrix(mod.beta, pars = c("beta_1[6]","beta_2[6]")))
mu_base <- inv_logit(alpha_base + b6_base[1]*so_val_scaled + b6_base[2]*so_val_scaled^2)
rho_base <- colMeans(as.matrix(mod.beta, pars = c("rho")))

x_grid <- seq(-a,1,length=200)
plot(x_grid, dbeta((x_grid + a) /(a+1),mu_base*rho_base,(1-mu_base)*rho_base)/(a+1), type = "l", ylim = 0:1, main = sp_test)
abline(v=0,col="red")

alpha_rho <- colMeans(as.matrix(mod.beta_rho, pars = c("alpha","alpha_rho")))
b6_rho <- colMeans(as.matrix(mod.beta_rho, pars = c("beta_1[6]","beta_2[6]","beta_rho_1[6]","beta_rho_2[6]")))
mu_rho <- inv_logit(alpha_rho[1] + b6_rho[1]*so_val_scaled + b6_rho[2]*so_val_scaled^2)
rho_rho <- 1000*inv_logit(alpha_rho[2] + b6_rho[3]*so_val_scaled + b6_rho[4]*so_val_scaled^2)

plot(x_grid, dbeta((x_grid + a) /(a+1),mu_rho*rho_rho,(1-mu_rho)*rho_rho)/(a+1), type = "l", ylim = 0:1, main = sp_test)
abline(v=0,col="red")




