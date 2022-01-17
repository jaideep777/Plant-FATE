setwd("/home/jaideep/codes/tmodel_cpp/tests/data/")
dat = read.csv("trait.csv")

library(dplyr)

dat = dat %>% mutate(P50..Mpa. = purrr::map_dbl(P50..Mpa., 
                                      .f = function(x) {
                                              ifelse(is.na(x), 
                                                  no = x,
                                                  yes = runif(1, -3.5, -0.5))
                                           }))

write.csv(dat, file = "trait_filled.csv")

