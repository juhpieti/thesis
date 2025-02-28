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

