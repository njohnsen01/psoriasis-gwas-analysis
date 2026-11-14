# Chen et al. (2021), Genet Epidemiol 45:577-592, doi:10.1002/gepi.22422.
# Three-term X omnibus: additive, female heterozygote dominance, additive x sex.
# ML likelihood-ratio implementation for a bounded validation stage.
# Not a mixed model; does not claim control of residual relatedness or rare-cell bias.
chen_test <- function(y, nuisance, g02, male, coding="02", reverse=FALSE) {
  stopifnot(length(y)==length(g02), length(y)==length(male),
            nrow(nuisance)==length(y), all(male %in% c(0,1)),
            all(y %in% c(0,1)))
  keep <- is.finite(g02) & complete.cases(nuisance)
  y <- y[keep]; s <- male[keep]; g <- g02[keep]
  C <- as.matrix(nuisance[keep,,drop=FALSE])
  storage.mode(C) <- "double"
  stopifnot(all(is.finite(C)), all(g %in% 0:2),
            all(g[s==1] %in% c(0,2)), length(unique(s))==2,
            length(unique(y))==2, qr(C)$rank==ncol(C),
            qr(cbind(C,1,s))$rank==ncol(C))
  if(reverse) g <- 2-g
  gd <- as.numeric(s==0 & g==1)
  if(coding=="01") g[s==1] <- g[s==1]/2
  stopifnot(coding %in% c("01","02"))
  X <- cbind(C,ga=g,gd=gd,ga_sex=g*s)
  if(qr(X)$rank-ncol(C)!=3) stop("Genetic design has fewer than 3 estimable terms")
  warnings <- character()
  fit <- function(M) withCallingHandlers(
    glm.fit(M,y,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100)),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
  f0 <- fit(C); f1 <- fit(X)
  if(!f0$converged || !f1$converged || any(!is.finite(f1$coefficients)) ||
     length(warnings)>0) stop(paste("Convergence/warning gate failed:",paste(warnings,collapse="; ")))
  stat <- f0$deviance-f1$deviance
  stopifnot(is.finite(stat),stat>= -1e-7,f1$rank-f0$rank==3)
  score <- max(abs(crossprod(X,y-f1$fitted.values)))
  if(score>1e-5) stop("Unpenalized score convergence gate failed")
  # Cell counts are retained internally. They must not be exported without review.
  cells <- table(sex=s,genotype=g02[keep],case=y)
  list(stat=max(0,stat),df=3,p=pchisq(max(0,stat),3,lower.tail=FALSE),
       n=length(y),cases=sum(y),score_max=score,
       min_observed_cell=min(cells[cells>0]),
       zero_cells=sum(cells==0), cells=cells,
       coefficients=f1$coefficients,iterations=f1$iter)
}
