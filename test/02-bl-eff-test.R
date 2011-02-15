#!/usr/bin/env Rscript

# Redirect output:
sink("02-bl-eff-test.log")
pdf("02-bl-eff-test.pdf")
# Load required libraries:
library(ape)
suppressMessages(require(methods,quietly=TRUE));
system("cd ..;make cat");

# Load package source:
source("../PCRcoalSource.R");

nr.cycles 	<- 30
nr.effs		<- 1000
nr.replications	<- 10000

expected.bl<-function(eff, nr.cycles ){
	nr.cycles * ( 1 - 1/(mean(eff) + 1) )
}

mean.bl <- numeric()
exp.bl	<- numeric()

for(i in 1:nr.effs){

	efficiencies<-runif(nr.cycles, min=0.5,max=1)	

	x<-PCRcoal(
        	initial.size    =1,
        	sample.size     =1,
        	nr.cycles       =nr.cycles,
		efficiencies    =efficiencies
	);

	branch.lengths<-numeric()

	for (i in 1:nr.replications){
		branch.lengths[i] <- sample.tree(x)$edge.length[1]
	}

	mean.bl <- c( mean.bl, mean(branch.lengths) )
	exp.bl <- expected.bl(efficiencies, nr.cycles)
}


plot(exp.bl, mean.bl)

reg<-lm(exp.bl ~ eff.seq)

summary(reg)

plot(reg)




