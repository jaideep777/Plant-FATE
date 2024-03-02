library(ggplot2)

data <- read.csv("cohort_vector.txt", header = FALSE, col.names = c("birthtime", "X0", "X1", "Z", "spName"))
data <- data[1:100,]

ggplot(data, aes(X0, X1, fill= Z)) + 
  geom_tile() +
  scale_fill_continuous(trans = "log") + 
  scale_x_continuous(trans = "log")
  
