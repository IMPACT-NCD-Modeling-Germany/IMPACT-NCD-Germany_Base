# Testing output syncing #

library(data.table)

n <- 100000
dt <- data.table(runif(n), rnorm(n), rgamma(n, 0.5))

fwrite(dt, "./outputs/testing_rsync.csv")


## Additional test analysis
dt <- data.table(runif(n), rnorm(n), rgamma(n, 0.9))



#####testing backsync#######

####more testing### testiiiing


####testtititiitititing