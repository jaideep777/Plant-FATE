library(tidyverse)

dat = read.delim("~/codes/Plant-FATE/evolution.txt")
traits = read.delim("~/codes/Plant-FATE/evol_traits.txt")
# traits$iter = traits$t

# dat1 = dat %>% filter(iter==1)
# traits1 = traits %>% filter(iter==1)
matplot(x=dat$t, y=dat[,-c(1,2)], type="l", lty=1, col=rainbow(length(which(traits$resident==1))/1000), lwd=1+2*traits$resident, ylim=c(0, 1))
# matplot(x=traits$t[traits$resident==1], y=traits$K[traits$resident==1], type="l", lty=1, col=scales::muted(rainbow(nrow(traits))))

# 
# traits %>% 
#   filter(resident==T) %>% 
#   ggplot(aes(y=y-dfy*0.01, x=x-dfx*0.01))+
#   theme_classic(base_size = 18)+
#   geom_point(aes(col=iter), alpha=0.7)+
#   scale_color_viridis_c(direction = -1)+
#   xlim(c(-1,1))+
#   ylim(c(-1,1))+
#   geom_point(data=traits %>% filter(iter==0), col="grey80", size=4, alpha=0.6)+
#   geom_point(data=traits %>% filter(iter==max(iter)), col="magenta4", size=3, alpha=1)+
#   geom_segment(data=traits %>% filter(iter==0), aes(xend=x, yend=y), arrow = arrow(length = unit(0.2, "cm")))+
#   ylab("y")+
#   xlab("x")
#   

traits %>% 
  filter(resident==T) %>% 
  ggplot(aes(y=y, x=x))+
  theme_classic(base_size = 18)+
  geom_point(aes(col=t), alpha=0.7)+
  scale_color_viridis_c(direction = -1)+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  geom_point(data=traits %>% filter(t==0 & resident==T), col="grey80", size=4, alpha=0.6)+
  geom_point(data=traits %>% filter(t==max(t) & resident==T), col="magenta4", size=3, alpha=1)+
  #geom_segment(data=traits %>% filter(t==0), aes(xend=x, yend=y), arrow = arrow(length = unit(0.2, "cm")))+
  ylab("y")+
  xlab("x")

