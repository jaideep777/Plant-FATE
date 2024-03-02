Ueq = function(x){
  (1-x)^2/(1+x)^4
}

U0 = function(x){
  (1-x)^2/(1+x)^4 + (1-x)/(1+x)^3
}

abm = read.delim("~/codes/libpspm/abm_init.txt", header = T)
h = hist(abm$x.0., plot = F, breaks=20)

x = seq(0,1, length.out=101)
plot(U0(x)~x, type="l")
points(y=h$counts*first(abm$u)/diff(h$breaks), x=h$mids, col="red")
#       ^ superinds * size = inds / breaks = density



dat = read.delim("~/codes/libpspm/abm_grid.txt", sep=" ", header=F)
abm = read.delim("~/codes/libpspm/abm_init_2d.txt", header = T)

p1 = tibble(x.0. = rnorm(n = 20000, mean=2, sd=sqrt(1/20)), 
       x.1. = rnorm(n = 20000, mean=4, sd=sqrt(1/20))) %>% ggplot(aes(x=x.0., y=x.1.)) +
  stat_bin2d(bins = 50) + 
  xlim(c(1,3)) + ylim(c(3,5)) + 
  ggtitle("expected")

p2 = abm %>% ggplot(aes(x=x.0., y=x.1.)) +
  stat_bin2d(bins = 50) + 
  xlim(c(1,3)) + ylim(c(3,5)) + 
  ggtitle("abm sample")
    
p3 = dat %>% ggplot(aes(x=V7, y=V8, fill=V11)) + 
  geom_tile() + 
  xlim(c(1,3)) + ylim(c(3,5)) + 
  ggtitle("grid")

cowplot::plot_grid(p1, p3, p2)


