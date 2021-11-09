install.packages('ratesci')
library(ratesci)

myout <- scaspci(x = 1:3, n = 9333, level = 0.999)
save(myout, file = paste0(".github/workflows/data_", make.names(Sys.time()), ".Rda"))
