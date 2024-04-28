#########################################################
## Code for the example in Sec. 5.1                    ##
## Holzmann, Klar (2024) Lancester correlation - a new ##
## dependence measure linked to maximum correlation    ##
#########################################################
oldw = getOption("warn")
options(warn = -1)
library(energy); library(TauStar); library(dHSIC); 
library(Ball); library(carData)
options(warn = oldw)

#########################################################################
#########################################################################
# function for computing different correlation measures
cor.est = function(xx) { #xx: nx2 matrix
  x = xx[,1]
  y = xx[,2]
  pearson = cor(x,y,method="pearson")
  spearman = cor(x,y,method="spearman")
  lc.lin = lcor(x, y, type = "linear")
  lc.rank = lcor(x, y, type = "rank")
  a = ace(x,y) 
  ace.cor =  cor(a$tx, a$ty)
  dc = dcor(x,y)
  tau = tauStarTest(x, y)$tStar
  xi = xicor(x,y)
  bcor = bcor(x, y)
  result = c(pearson, spearman, lc.lin, lc.rank, ace.cor, dc, tau, xi, bcor)
  names(result) = c("pearson", "spearman", "lc.lin", "lc.rank", "acecor", 
                    "dcor", "tau", "xi", "bcor")
  return( result )
}

# function for computing different independence tests
cor.tests = function(xx, nperm=999) { # xx: nx2-Matrix
  x = xx[, 1]
  y = xx[,2]
  pearson = perm.relation(x, y, method = "pearson", R = nperm)$p.value
  spearman = perm.relation(x, y, method = "spearman", R = nperm)$p.value
  lc.lin = lcor.test(x, y, type = "linear", method = "permutation", nperm = nperm)$pval
  lc.rank1 = lcor.test(x, y, type = "rank", method = "asymptotic")$pval
  lc.rank2 = lcor.test(x, y, type = "rank", method = "permutation", nperm = nperm)$pval
  ace = ace.test(x,y, nperm=nperm)$pval
  dc = dcor.test(x, y, R = nperm)$p.value
  tau = tauStarTest(x, y, mode = "permutation", resamples = nperm)$pVal
  xi = xicor(x, y, pvalue = TRUE, method = "permutation", nperm = nperm)$pval
  hsic = dhsic.test(x, y, method="permutation", kernel="gaussian", B=nperm)$p.value
  bcov = bcov.test(x, y, method="permutation", num.permutations=nperm)$p.value
  result = c(pearson, spearman, lc.lin, lc.rank1, lc.rank2, ace, dc, tau, xi, hsic, bcov)
  names(result) = c("pear", "spear", "lc.lin", "lc.rank1", "lc.rank2", "ace", 
                    "dcor", "tau", "xi", "hsic", "bcov")
  return( result )
}

############################################################################
# Fox J. and Weisberg, S. (2019) An R Companion to Applied Regression, 
# Third Edition, Sage, Section 10.7.1
data(Salaries, package="carData")
head(Salaries)
summary(Salaries)
ftable( xtabs( ~ discipline + rank + sex, data=Salaries) )

data = Salaries[Salaries$rank == "Prof" & Salaries$sex == "Male", ]
(n = nrow(data))
x = data$yrs.service
y = data$salary
lcor.comp(x, y, type = "rank", plot = TRUE)
xx = cbind( data$salary, data$yrs.service)
round( cor.est(xx), 2)
set.seed(24)
round( cor.tests(xx, nperm = 9999), 3)

data = Salaries[Salaries$rank == "Prof" & Salaries$sex == "Male" & Salaries$discipline == "A", ]
(n = dim(data)[1])
x = data$yrs.service
y = data$salary
xx = cbind( data$salary, data$yrs.service)
round( cor.est(xx), 2)
set.seed(24)
round( cor.tests(xx, nperm = 9999), 3)

data = Salaries[Salaries$rank == "Prof" & Salaries$sex == "Male" & Salaries$discipline == "B", ]
(n = dim(data)[1])
x = data$yrs.service
y = data$salary
xx = cbind( data$salary, data$yrs.service)
round( cor.est(xx), 2)
set.seed(24)
round( cor.tests(xx, nperm = 9999), 3)
###########################################################################