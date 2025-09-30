# Testing output syncing #

library(data.table)

n <- 100000
dt <- data.table(runif(n), rnorm(n))

fwrite(dt, "./outputs/testing_rsync.csv")