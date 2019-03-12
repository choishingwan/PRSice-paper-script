#!/usr/bin/env Rscript
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
files=list.files(pattern="result$")
res <- NULL
for(i in files){
	tmp <- read.table(i, header=T, sep=" ")
	tmp$Valid.P <- as.numeric(as.character(tmp$Valid.P))
	tmp$Target.P <- as.numeric(as.character(tmp$Target.P))
	res <- rbind(res, tmp)
}
res %>% 
	filter(!(Heritability==0  & Program=="LDPred")) %>%
	filter(!(Heritability==0 & Program=="lassosum")) -> res
print(table(res$Program))
write.csv(res, paste0("../",args[1],".result.csv"), quote=F, row.names=F)
