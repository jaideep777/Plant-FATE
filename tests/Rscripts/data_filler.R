setwd("~/codes/Plant-FATE/tests/data/")
dat = read.csv("trait.csv")

library(dplyr)

dat = dat %>% mutate(P50..Mpa. = purrr::map_dbl(P50..Mpa., 
                                      .f = function(x) {
                                              ifelse(is.na(x), 
                                                  no = x,
                                                  yes = runif(1, -5.5, -2.5))
                                           }))

write.csv(dat, file = "Amz_trait_filled_HD.csv", row.names = F)

# CWM traits

dat %>% select(Leaf.LMA..g.m2., Number.of.individuals) %>% na.omit %>% summarize(meanLMA = sum(Leaf.LMA..g.m2.*Number.of.individuals)/sum(Number.of.individuals))
dat %>% select(meanWoodDensity..g.cm3., Number.of.individuals) %>% na.omit %>% summarize(meanWD = 1000*sum(meanWoodDensity..g.cm3.*Number.of.individuals)/sum(Number.of.individuals))
dat %>% select(P50..Mpa., Number.of.individuals) %>% na.omit %>% summarize(meanP50 = sum(P50..Mpa.*Number.of.individuals)/sum(Number.of.individuals))

par(mfrow = c(2,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
h = hist(dat$X2r.diameter, breaks=20)
plot(h$counts~h$mids, log="xy")

# Reduce diversity

dat = read.csv("Amz_trait_filled_HD.csv")
dat = dat[1:100,]
qlma = quantile(dat$Leaf.LMA..g.m2., na.rm=T)
qwd =  quantile(dat$meanWoodDensity..g.cm3., na.rm=T)
qhmat = quantile(dat$Height_Max.m., na.rm=T)
qp50 = quantile(dat$P50..Mpa., na.rm=T)

dat_reduced = dat %>% 
  filter((Leaf.LMA..g.m2. < qlma[4] & Leaf.LMA..g.m2. > qlma[2]) | is.na(Leaf.LMA..g.m2.)) %>% 
  filter((meanWoodDensity..g.cm3. < qwd[4] & meanWoodDensity..g.cm3. > qwd[2]) | is.na(meanWoodDensity..g.cm3.) ) %>% 
  filter((Height_Max.m. < qhmat[4] & Height_Max.m. > qhmat[2]) | is.na(Height_Max.m.) ) %>% 
  filter((P50..Mpa. < qp50[4] & P50..Mpa. > qp50[2]) | is.na(P50..Mpa.) )  

write.csv(dat_reduced, file = "Amz_trait_filled_LD.csv", row.names = F)
