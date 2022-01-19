setwd("/home/jaideep/codes/tmodel_cpp/tests/data/")
dat = read.csv("trait_100.csv")

library(dplyr)

dat = dat %>% mutate(P50..Mpa. = purrr::map_dbl(P50..Mpa., 
                                      .f = function(x) {
                                              ifelse(is.na(x), 
                                                  no = x,
                                                  yes = runif(1, -3.5, -0.5))
                                           }))

write.csv(dat, file = "trait_100_filled.csv", row.names = F)

# CWM traits

dat %>% select(Leaf.LMA..g.m2., Number.of.individuals) %>% na.omit %>% summarize(meanLMA = sum(Leaf.LMA..g.m2.*Number.of.individuals)/sum(Number.of.individuals))
dat %>% select(meanWoodDensity..g.cm3., Number.of.individuals) %>% na.omit %>% summarize(meanWD = 1000*sum(meanWoodDensity..g.cm3.*Number.of.individuals)/sum(Number.of.individuals))
dat %>% select(P50..Mpa., Number.of.individuals) %>% na.omit %>% summarize(meanP50 = sum(P50..Mpa.*Number.of.individuals)/sum(Number.of.individuals))

par(mfrow = c(2,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
h = hist(dat$X2r.diameter, breaks=20)
plot(h$counts~h$mids, log="xy")
