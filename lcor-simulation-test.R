############################################################
## Code for simulation of independence tests in Sec. 4.2  ##
## Holzmann, Klar (2024) Lancester correlation - a new    ##
## dependence measure linked to maximum correlation       ##
############################################################
oldw = getOption("warn")
options(warn = -1)
library(energy); library(TauStar); library(dHSIC); library(Ball)
library(mvtnorm); library(TSA); library(xtable)
options(warn = oldw)
start_time = Sys.time()

#########################################################################
#########################################################################
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

#########################################################################
#########################################################################
generate.sample = function(n, distr) {
  flag = TRUE  #valid distribution
  if (distr=="bvnorm1") { # Bivariate normal 1
    sigma = matrix(c(1,0,0,4), ncol=2) 
    xx = rmvnorm( n, mean=c(2,7), sigma=sigma)
  }
  else if (distr=="bvnorm2") { # Bivariate normal 2
    sigma = matrix(c(1,0.5*2,0.5*2,4), ncol=2) 
    xx = rmvnorm( n, mean=c(2,7), sigma=sigma)
  }
  else if (distr=="bvnorm3") { # Bivariate normal 3
    sigma = matrix(c(1,0.95*2,0.95*2,4), ncol=2) 
    xx = rmvnorm( n, mean=c(2,7), sigma=sigma)
  }
  else if (distr=="stnormix1") { # Bivariate standard normal mixture 1
    comp = rbinom(n, 1, 0.5) 
    n1 = length(comp[comp==1])
    if (n1 == 0) n1 = 1
    if (n1 == n) n1 = n-1
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="stnormix2") { # Bivariate standard normal mixture2
    comp = rbinom(n, 1, 1/3) 
    n1 = length(comp[comp==1])
    if (n1 == 0) n1 = 1
    if (n1 == n) n1 = n-1
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="stnormix3") { # Bivariate standard normal mixture 3
    p = 1/4
    comp = rbinom(n, 1, p) 
    n1 = length(comp[comp==1])
    if (n1 == 0) n1 = 1
    if (n1 == n) n1 = n-1
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="normix1") { # Bivariate normal mixture 1
    comp = rbinom(n, 1, 0.5)
    n1 = length(comp[comp==1])
    x1 = rmvnorm(n1, mean=c(0,0))
    x2 = rmvnorm(n-n1, mean=c(5,5))
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="normix2") { # Bivariate normal mixture 2
    comp = sample(1:4, n, replace = TRUE)
    x1 = x2 = x3 = x4 = numeric()
    if (length(comp[comp==1])>0) x1 = rmvnorm(length(comp[comp==1]), mean=c(0,0))
    if (length(comp[comp==2])>0) x2 = rmvnorm(length(comp[comp==2]), mean=c(0,5))
    if (length(comp[comp==3])>0) x3 = rmvnorm(length(comp[comp==3]), mean=c(5,0))
    if (length(comp[comp==4])>0) x4 = rmvnorm(length(comp[comp==4]), mean=c(5,5))
    xx = rbind(x1,x2,x3,x4)
    xx = xx[sample(n),]
    }
  else if (distr=="tdist1") { # Bivariate t distribution 1
    xx = rmvt(n, sigma = diag(2), df = 5)
    }
  else if (distr=="tdist2") { # Bivariate t distribution 2
    xx = rmvt(n, sigma = diag(2), df = 2)
  }
  else if (distr=="tdist3") { # Bivariate t distribution 3
    xx = rmvt(n, sigma = diag(2), df = 1)
  }
  else if (distr=="tdist4") { # Bivariate t distribution 4
    xx <- rmvt(n, sigma = matrix( c(1,0.2,0.2,1), nrow=2), df = 5)
  }
  else if (distr=="tdist5") { # Bivariate t distribution 5
    xx <- rmvt(n, sigma = matrix( c(1,0.2,0.2,1), nrow=2), df = 2)
  }
  else if (distr=="tdist6") { # Bivariate t distribution 6
    xx <- rmvt(n, sigma = matrix( c(1,0.2,0.2,1), nrow=2), df = 1)
  }
  else if (distr=="tdist7") { # Bivariate t distribution 7
    xx <- rmvt(n, sigma = matrix( c(1,0.5,0.5,1), nrow=2), df = 5)
  }
  else if (distr=="tdist8") { # Bivariate t distribution 8
    xx <- rmvt(n, sigma = matrix( c(1,0.5,0.5,1), nrow=2), df = 2)
  }
  else if (distr=="tdist9") { # Bivariate t distribution 9
    xx <- rmvt(n, sigma = matrix( c(1,0.5,0.5,1), nrow=2), df = 1)
  }
  else if (distr=="disc.unif") { # Uniform distribution on the unit disc
    rho = sqrt(runif(n))
    theta = runif(n, 0, 2*pi)
    x = rho * cos(theta); y = rho * sin(theta)
    xx = cbind(x,y)
  }
  else if (distr=="rhomb.unif") { # Uniform distribution on rhomb
    u = runif(n); v = runif(n)
    x = (u-v+1)/2; y = (u+v)/2
    xx = cbind(x,y)
  }
  else if (distr=="triangle.unif") { # Uniform distribution on Triangle
    u = runif(n); v = runif(n); 
    x = pmin(u,v); y = 1-pmax(u,v)
    xx = cbind(x,y)
  }
  else if (distr=="garch") { # garch acf
    xdata = garch.sim(alpha=c(.01,.6), beta=0.2, n=n+1)
    x = xdata[1:(n-1)]
    y = xdata[2:n]
    xx = cbind(x,y)
  }
  else if (distr=="reg.lin1") { # linear regression 1
    x = runif(n,0,1)
    #y = x
    y = x + rnorm(n,0,0.15)
    xx = cbind(x,y)
  }
  else if (distr=="reg.lin2") { # linear regression 2
    x = runif(n,0,1)
    y = x + rnorm(n,0,0.3)
    xx = cbind(x,y)
  }
  else if (distr=="reg.lin3") { # linear regression 3
    x = runif(n,0,1)
    y = x + rnorm(n,0,0.45)
    xx = cbind(x,y)
  }
  else if (distr=="reg.quad1") { # quadratic regression 1
    x = runif(n,-1,1)
    y = x^2
    xx = cbind(x,y)
  }
  else if (distr=="reg.quad2") { # quadratic regression 2
    x = runif(n,-1,1)
    y = x^2 + rnorm(n,0,0.15)
    xx = cbind(x,y)
  }
  else if (distr=="reg.quad3") { # quadratic regression 3
    x = runif(n,-1,1)
    y = x^2 + rnorm(n,0,0.3)
    xx = cbind(x,y)
  }
  else if (distr=="reg.trig1") { # trigonometric regression 1
    x = runif(n,0,4*pi)
    y = (sin(x)+1)/2
    xx = cbind(x,y)
  }
  else if (distr=="reg.trig2") { # trigonometric regression 2
    x = runif(n,0,4*pi)
    y = (sin(x)+1)/2 + rnorm(n,0,0.15)
    xx = cbind(x,y)
  }
  else if (distr=="reg.trig3") { # trigonometric regression 3
    x = runif(n,0,4*pi)
    y = (sin(x)+1)/2 + rnorm(n,0,0.3)
    xx = cbind(x,y)
  }
  else flag = FALSE
  if (flag == TRUE) return(xx) else stop("distribution not defined")
}
############################################################################
############################################################################

# simulation function
cor.sim = function(n, B, distr, alpha = 0.05) {
  counter = rep(0, 10)
  for (i in 1:B) {
    xx = generate.sample(n, distr)
    res = cor.tests(xx)
    counter = ifelse(res <= alpha, counter + 1, counter)
  }
  return( counter / B )
}

############################################################################
############################################################################
set.seed(28)
n = 100
B = 1000
#distr.vec = "stnormix3"
distr.vec = c(
  "bvnorm1", "bvnorm2", "bvnorm3", "stnormix1", "stnormix2", "stnormix3", "normix2",
  "tdist1", "tdist2", "tdist3", "tdist4", "tdist5", "tdist6", 
  "disc.unif", "rhomb.unif", "triangle.unif", "garch", 
  "reg.lin2", "reg.lin3", "reg.quad2", "reg.quad3", "reg.trig2", "reg.trig3")

n.distr = length(distr.vec)
ncol = 11 #number of tests
result.tab =  matrix(0, nrow=n.distr, ncol=ncol+1)
colnames(result.tab) = c("distribution", "pear", "spear", "mcl", "mcn1", "mcn2", "ace",
                         "dcor", "tau", "xi", "hsic", "bcov")
result.tab = as.data.frame(result.tab)

i = 1  
for (distr in distr.vec) {
  result.tab[i,1] = distr
  result.tab[i,2:(ncol+1)] = round( cor.sim(n, B, distr), 3)
  print(i)
  i = i+1
}

result.tab
#result.tab[,1] = "MN3"
result.tab[,1] = c("BVN(0)", "BVN(0.5)",  "BVN(0.95)", "MN1", "MN2", "MN3", "MN",
                   "BVT5(0)", "BVT2(0)", "BVT1(0)", "BVT5(0.2)", "BVT2(0.2)", "BVT1(0.2)",
                   "UnifDisc", "UnifDrhomb", "UnifTriangle", "GARCH(2,1)", 
                   "RegLin1", "RegLin2", "RegQuad1", "RegQuad2", "RegTrig1", "RegTrig2")
print( xtable(result.tab), include.rownames=FALSE)

end_time = Sys.time()
print(end_time - start_time)
#########################################################################
#########################################################################