args <- commandArgs(TRUE)
root <- if(length(args)) args[1] else "."
source(file.path(root,"chen_model.R"))
set.seed(20261013)
n <- 10000L
male <- rbinom(n,1,.45)
pc <- rnorm(n,mean=male*.25)
g <- ifelse(male==1,2*rbinom(n,1,.18),rbinom(n,2,.18))
C <- cbind(intercept=1,sex=male,PC1=pc,sex_PC1=pc*male)
eta <- -2.6+.25*male+.3*pc+.2*pc*male+.22*g+.2*(g==1 & male==0)-.12*g*male
y <- rbinom(n,1,plogis(eta))
fits <- lapply(c(FALSE,TRUE),function(rev)lapply(c("02","01"),function(code)
 chen_test(y,C,g,male,code,rev)))
stats <- unlist(lapply(fits,function(x)vapply(x,`[[`,0.0,"stat")))
stopifnot(max(stats)-min(stats)<1e-6)
# Independent formula fit verifies the intended nested model and reference df.
d <- data.frame(y,male,pc,g,gd=as.numeric(g==1 & male==0))
null <- glm(y~male*pc,family=binomial(),data=d)
full <- glm(y~male*pc+g+gd+g:male,family=binomial(),data=d)
reference_stat <- null$deviance-full$deviance
stopifnot(abs(stats[1]-reference_stat)<1e-6)
# Structural inability to estimate dominance must fail, not silently become 2 df.
bad_g <- 2*rbinom(n,1,.2)
stopifnot(inherits(try(chen_test(y,C,bad_g,male),silent=TRUE),"try-error"))
receipt <- data.frame(check=c("male coding invariance","allele reversal invariance",
 "independent formula comparison","rank deficient dominance rejected"),passed=TRUE)
write.csv(receipt,file.path(root,"synthetic_validation.csv"),row.names=FALSE)
writeLines(c("Synthetic software validation only; no All of Us data analyzed.",
 sprintf("Maximum statistic discrepancy: %.12g",max(stats)-min(stats)),
 capture.output(sessionInfo())),file.path(root,"synthetic_validation.txt"))
print(receipt)
