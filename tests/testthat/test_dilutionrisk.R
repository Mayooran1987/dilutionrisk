library(dilutionrisk)

set.seed(1234)

n <- 100
meanlog <-  0
sdlog <-  1
a <- 0
b <- 500
rtrunpoilog(n, meanlog, sdlog, a, b)


lambda <- 20
a <- 0
b <- 300
FDF <- 1000
USL <- 1000
n_sim <- 10000
prob_detection_homogeneous(lambda, a, b, FDF, USL, n_sim)


meanlog <- 20
sdlog <- 0.8
a <- 0
b <- 300
FDF <- 1000
USL <- 1000
n_sim <- 10000
prob_detection_heterogeneous(meanlog, sdlog, a, b, FDF, USL, n_sim)

c <- 0
lambda <- 10
a <- 0
b <- 300
FDF <- 1000
USL <- 1000
n <- 5
n_sim <- 10000
prob_acceptance_homogeneous(c, lambda, a, b, FDF, USL, n, n_sim)

