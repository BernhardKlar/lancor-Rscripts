#########################################################
## Code for computing Table 1 in                       ##
## Holzmann, Klar (2024) Lancester correlation - a new ##
## dependence measure linked to maximum correlation    ##
#########################################################
oldw = getOption("warn")
options(warn = -1)
library(acepack); library(energy); library(XICOR)
library(TauStar); library(Ball)
library(mvtnorm); library(TSA); library(xtable)
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

#########################################################################
#########################################################################
set.seed(21)
n = 1e4
#########################################################################
#########################################################################
#Bivariate standard normal mixture
comp = rbinom(n, 1, 0.5) 
n1 = length(comp[comp==1])
sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
x1 = rmvnorm( n1, sigma=sigma1)
sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
x2 = rmvnorm( n-n1, sigma=sigma2)
xx = rbind(x1,x2)
xx = xx[sample(n),]
res1 = cor.est(xx) 
#########################################################################
# Bivariate standard normal mixture
comp = rbinom(n, 1, 1/3) 
n1 = length(comp[comp==1])
sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
x1 = rmvnorm( n1, sigma=sigma1)
sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
x2 = rmvnorm( n-n1, sigma=sigma2)
xx = rbind(x1,x2)
xx = xx[sample(n),]
res2 = cor.est(xx)  
#########################################################################
#########################################################################
# t distribution examples
xx = rmvt(n, sigma = diag(2), df = 5)
res3 = round(cor.est(xx), 2)
xx = rmvt(n, sigma = diag(2), df = 1)
res4 = cor.est(xx) 
############################################################################
############################################################################
# Uniform distribution on the unit disc
rho = sqrt(runif(n))
theta = runif(n, 0, 2*pi)
x = rho * cos(theta); y = rho * sin(theta)
theta = seq(0,2*pi,0.001); lines( cos(theta), sin(theta), lwd=2)
res5 = cor.est(cbind(x,y)) 
############################################################################
############################################################################
#returns from GARCH
xdata = garch.sim(alpha=c(.01,.6), beta=0.2, n=n)
x = xdata[1:(n-1)]
y = xdata[2:n]
res6 = cor.est(cbind(x,y))
############################################################################
############################################################################
## Bivariate normal
set.seed(21)
sigma = matrix(c(1,0.1,0.1,1), ncol=2)
xx = rmvnorm( n, sigma=sigma)
res0a = cor.est(xx) #theoretical: cor=maxcor=0
sigma = matrix(c(1,0.3,0.3,1), ncol=2)
xx = rmvnorm( n, sigma=sigma)
res0b = cor.est(xx) #theoretical: cor=maxcor=0
sigma = matrix(c(1,0.95,0.95,1), ncol=2)
xx = rmvnorm( n, sigma=sigma)
res0c = cor.est(xx) #theoretical: cor=maxcor=0
#########################################################################
#########################################################################
result.tab = rbind(res0a, res0b, res1, res2, res3, res4, res5, res6)
#print( xtable( result.tab, digits=3), include.rownames=FALSE)

