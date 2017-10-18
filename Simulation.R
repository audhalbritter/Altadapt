trange <- 1:15
slope <- plyr::rdply(10, {
  slope <- sapply(trange, function(x){
    y <- rnorm(10)
    x <- seq(0, x, length = 10)
    coef(lm(y ~ x))[2]
  })
  data_frame(range = trange, slope = slope)
})

ggplot(slope, aes(x = range, y = abs(slope), group = range)) +
  geom_boxplot() +
  scale_y_log10()

y <- rnorm(10)
x <- seq(0, 15, length = 10)
cor <- cor(x, y)
ggplot(aes(x = trange[1], y = cor))


result <- plyr::rdply(10, {
  corellation <- sapply(trange, function(x){
    y <- rnorm(10)
    x <- seq(0, x, length = 10)
    cor(x, y)
  })
  data_frame(range = trange, corellation = corellation)
})

ggplot(result, aes(x = range, y = corellation, group = range)) +
  geom_boxplot()


