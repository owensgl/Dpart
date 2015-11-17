inputfile <- commandArgs(trailingOnly = TRUE)
library(bootstrap)
IN <- read.table ("/home/owens/exil/2015/Dpart/cenann.calann.G100.G115_out_block.txt",header=T)
IN <- read.table(inputfile[1], header=T)

#Column numbers for data
d1num = 3
d1den = 4
d2num = 5
d2den = 6
d12num = 7
d12den = 8
#Jackknife for D1 statistic
D1_num <- function(x, xdata){sum(xdata[x,d1num])}
D1_denom <- function(x, xdata){sum(xdata[x,d1den])}
D1.num.results <- jackknife(1:nrow(IN), D1_num, IN)
D1.denom.results <- jackknife(1:nrow(IN), D1_denom, IN)
D1.jacked <- (((D1.num.results$jack.values/D1.denom.results$jack.values) -
                  sum(IN[,d1num])/sum(IN[,d1den]))^2)
StdErr.D1 <- sqrt(((length(D1.jacked)-1)/length(D1.jacked))*sum(D1.jacked))
D1.total <- sum(IN[,d1num])/sum(IN[,d1den])
D1.Zscore <- D1.total/StdErr.D1
D1.pvalue <- 2*pnorm(-abs(D1.Zscore))

#Jackknife for D2 statistic
D2_num <- function(x, xdata){sum(xdata[x,d2num])}
D2_denom <- function(x, xdata){sum(xdata[x,d2den])}
D2.num.results <- jackknife(1:nrow(IN), D2_num, IN)
D2.denom.results <- jackknife(1:nrow(IN), D2_denom, IN)
D2.jacked <- (((D2.num.results$jack.values/D2.denom.results$jack.values) -
                 sum(IN[,d2num])/sum(IN[,d2den]))^2)
StdErr.D2 <- sqrt(((length(D2.jacked)-1)/length(D2.jacked))*sum(D2.jacked))
D2.total <- sum(IN[,d2num])/sum(IN[,d2den])
D2.Zscore <- D2.total/StdErr.D2
D2.pvalue <- 2*pnorm(-abs(D2.Zscore))

#Jackknife for D12 statistic
D12_num <- function(x, xdata){sum(xdata[x,d12num])}
D12_denom <- function(x, xdata){sum(xdata[x,d12den])}
D12.num.results <- jackknife(1:nrow(IN), D12_num, IN)
D12.denom.results <- jackknife(1:nrow(IN), D12_denom, IN)
D12.jacked <- (((D12.num.results$jack.values/D12.denom.results$jack.values) -
                 sum(IN[,d12num])/sum(IN[,d12den]))^2)
StdErr.D12 <- sqrt(((length(D12.jacked)-1)/length(D12.jacked))*sum(D12.jacked))
D12.total <- sum(IN[,d12num])/sum(IN[,d12den])
D12.Zscore <- D12.total/StdErr.D12
D12.pvalue <- 2*pnorm(-abs(D12.Zscore))

test <- c("D1", "D2", "D12")
value <- c(D1.total, D2.total, D12.total)
stderr <- c(StdErr.D1, StdErr.D2, StdErr.D12)
zscore <- c(D1.Zscore, D2.Zscore, D12.Zscore)
pvalue <- c(D1.pvalue, D2.pvalue, D12.pvalue)
output <- data.frame(test, value, stderr, zscore, pvalue)

write.table(output, file=inputfile[2], quote=F, sep="\t", row.names=F, col.names=T)
