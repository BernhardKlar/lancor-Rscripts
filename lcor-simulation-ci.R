#############################################################
## Code for simulation of confidence intervals in Sec. 4.1 ##
## Holzmann, Klar (2024) Lancester correlation - a new     ##
## dependence measure linked to maximum correlation        ##
#############################################################
oldw = getOption("warn")
options(warn = -1)
library(mvtnorm); library(TSA); library(xtable)
options(warn = oldw)
start_time <- Sys.time()

#########################################################################
#########################################################################

# function for computing confidence intervals
cor.ci = function(xx, conf.level = 0.95) { # xx: nx2-Matrix
  ci.lin.p1 = round( lcor.ci(xx, type = "linear", method="plugin", con = FALSE), 4)
  ci.lin.p2 = round( lcor.ci(xx, type = "linear", method="plugin", con = TRUE), 4)
  ci.lin.b1 = round( lcor.ci(xx, type = "linear",method="boot", con = FALSE), 4)
  ci.lin.b2 = round( lcor.ci(xx, type = "linear",method="boot", con = TRUE), 4)
  ci.rank1 = round( lcor.ci(xx, type = "rank", con = FALSE), 4)
  ci.rank2 = round( lcor.ci(xx, type = "rank", con = TRUE), 4)
  #                     
  return(list( ci.lin.p1, ci.lin.p2, ci.lin.b1, ci.lin.b2, ci.rank1, ci.rank2))
}

#xx = generate.sample(200, "bvnorm2")
#cor.ci(xx)

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
    x1 = x2 = numeric()
    comp = rbinom(n, 1, 0.5) 
    n1 = length(comp[comp==1])
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    if (n1>0) x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    if (n-n1>0) x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="stnormix2") { # Bivariate standard normal mixture2
    x1 = x2 = numeric()
    comp = rbinom(n, 1, 1/3) 
    n1 = length(comp[comp==1])
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    if (n1>0) x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    if (n-n1>0) x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="stnormix3") { # Bivariate standard normal mixture 3
    x1 = x2 = numeric()
    p = 1/4
    comp = rbinom(n, 1, p) 
    n1 = length(comp[comp==1])
    sigma1 = matrix(c(1,-0.5,-0.5,1), ncol=2)
    if (n1>0) x1 = rmvnorm( n1, sigma=sigma1)
    sigma2 = matrix(c(1,0.5,0.5,1), ncol=2)
    if (n-n1>0) x2 = rmvnorm( n-n1, sigma=sigma2)
    xx = rbind(x1,x2)
    xx = xx[sample(n),]
  }
  else if (distr=="normix1") { # Bivariate normal mixture 1
    x1 = x2 = numeric()
    comp = rbinom(n, 1, 0.5)
    n1 = length(comp[comp==1])
    if (n1>0) x1 = rmvnorm(n1, mean=c(0,0))
    if (n-n1>0) x2 = rmvnorm(n-n1, mean=c(5,5))
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
    x = xdata[1:n]
    y = xdata[2:(n+1)]
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

value.rho = function(distr) { #based on sample of size 1e7
  flag = TRUE  #valid distribution
  if (distr=="bvnorm1") {lc.lin = 0.0; lc.rank = 0.0} else
  if (distr=="bvnorm2") {lc.lin = 0.5; lc.rank = 0.5} else
  if (distr=="bvnorm3") {lc.lin = 0.95; lc.rank = 0.95} else
  if (distr=="stnormix1") {lc.lin = 0.25; lc.rank = 0.25} else
  if (distr=="stnormix2") {lc.lin = 0.25; lc.rank = 0.25} else
  if (distr=="stnormix3") {lc.lin = 0.25; lc.rank = 0.25} else
  if (distr=="normix1") {lc.lin = 0.86; lc.rank = 0.64} else
  if (distr=="normix2") {lc.lin = 0.0; lc.rank = 0.0} else
  if (distr=="tdist1") {lc.lin = 1; lc.rank = 0.219} else #t5
  if (distr=="tdist2") {lc.lin = 1; lc.rank = 0.488} else
  if (distr=="tdist3") {lc.lin = 1; lc.rank = 0.723} else
  if (distr=="tdist4") {lc.lin = 1; lc.rank = 0.247} else
  if (distr=="tdist5") {lc.lin = 1; lc.rank = 0.504} else
  if (distr=="tdist6") {lc.lin = 1; lc.rank = 0.730} else
  if (distr=="tdist7") {lc.lin = 1; lc.rank = 0.495} else
  if (distr=="tdist8") {lc.lin = 1; lc.rank = 0.591} else
  if (distr=="tdist9") {lc.lin = 1; lc.rank = 0.771} else
  if (distr=="disc.unif") {lc.lin = 0.333; lc.rank = 0.267} else
  if (distr=="rhomb.unif") {lc.lin = 0.428; lc.rank = 0.335} else
  if (distr=="triangle.unif") {lc.lin = 0.500; lc.rank = 0.476} else
  if (distr=="garch") {lc.lin = 0.52; lc.rank = 0.52} else
  if (distr=="reg.lin1") {lc.lin = 0.887; lc.rank = 0.858} else
  if (distr=="reg.lin2") {lc.lin = 0.693; lc.rank = 0.675} else
  if (distr=="reg.lin3") {lc.lin = 0.540; lc.rank = 0.527} else
  if (distr=="reg.quad1") {lc.lin = 0.598; lc.rank = 0.529} else
  if (distr=="reg.quad2") {lc.lin = 0.423; lc.rank = 0.376} else
  if (distr=="reg.quad3") {lc.lin = 0.238; lc.rank = 0.255} else
  if (distr=="reg.trig1") {lc.lin = 0.390; lc.rank = 0.346} else
  if (distr=="reg.trig2") {lc.lin = 0.359; lc.rank = 0.338} else
  if (distr=="reg.trig3") {lc.lin = 0.297; lc.rank = 0.296} else flag = FALSE
  if (flag == TRUE) return( c(lc.lin, lc.rank) ) else 
    stop("distribution not defined")
}


############################################################################
############################################################################

# simulation function
cor.sim = function(n, B, distr, conf.level = 0.95) {
  counter = rep(0, 6)
  mean.length = rep(0, 6)
  lc.lin = value.rho(distr)[1]
  lc.rank = value.rho(distr)[2]
  
  for (i in 1:B) {
    xx = generate.sample(n, distr)
    res = cor.ci(xx, conf.level)
    if (res[[1]][1] <= lc.lin & lc.lin <= res[[1]][2]) counter[1] = counter[1] + 1
    if (res[[2]][1] <= lc.lin & lc.lin <= res[[2]][2]) counter[2] = counter[2] + 1
    if (res[[3]][1] <= lc.lin & lc.lin <= res[[3]][2]) counter[3] = counter[3] + 1
    if (res[[4]][1] <= lc.lin & lc.lin <= res[[4]][2]) counter[4] = counter[4] + 1
    if (res[[5]][1] <= lc.rank & lc.rank <= res[[5]][2]) counter[5] = counter[5] + 1
    if (res[[6]][1] <= lc.rank & lc.rank <= res[[6]][2]) counter[6] = counter[6] + 1
    
    mean.length[1] = mean.length[1] + res[[1]][2] - res[[1]][1]
    mean.length[2] = mean.length[2] + res[[2]][2] - res[[2]][1]
    mean.length[3] = mean.length[3] + res[[3]][2] - res[[3]][1]
    mean.length[4] = mean.length[4] + res[[4]][2] - res[[4]][1]
    mean.length[5] = mean.length[5] + res[[5]][2] - res[[5]][1]
    mean.length[6] = mean.length[6] + res[[6]][2] - res[[6]][1]
  }
  
  return( list(counter/B, mean.length/B) )
}

#cor.sim(200, 10, "bvnorm1", conf.level = 0.95)

############################################################################
############################################################################
set.seed(28)
n = 200
B = 1000 
#distr.vec = c("bvnorm1")
distr.vec = c(
  "bvnorm1", "bvnorm2", "bvnorm3", "stnormix1", "stnormix2", "stnormix3", 
  "normix2", "tdist1", "tdist2", "tdist3", "tdist4", "tdist5", "tdist6", 
  "disc.unif", "rhomb.unif", "triangle.unif", "reg.lin2", "reg.lin3", 
  "reg.quad2", "reg.quad3", "reg.trig2", "reg.trig3")

n.distr = length(distr.vec)
ncol = 6
result.tab = matrix(0, nrow=n.distr, ncol=ncol+1)
colnames(result.tab) = c("distribution", "lin.plugin", "lin.plugin.con", 
                       "lin.boot", "lin.boot.con", "rank.boot", "rank.boot.con")
result.tab = as.data.frame(result.tab)
mean.tab = result.tab
  
i = 1  
for (distr in distr.vec) {
  print(i)
  result.tab[i,1] = mean.tab[i,1] = distr
  result = cor.sim(n, B, distr)
  result.tab[i,2:(ncol+1)] = round( result[[1]], 3)
  mean.tab[i,2:(ncol+1)] = round( result[[2]], 3)
  i = i+1
}

result.tab
#result.tab[,1] = "BVN(0)"
result.tab[,1] = c("BVN(0)", "BVN(0.5)",  "BVN(0.95)", "MN1", "MN2", "MN3", "MN",
        "BVT5(0)", "BVT2(0)", "BVT1(0)", "BVT5(0.2)", "BVT2(0.2)", "BVT1(0.2)",
        "UnifDisc", "UnifRhomb", "UnifTriangle", "RegLin1", "RegLin2", 
        "RegQuad1", "RegQuad2", "RegTrig1", "RegTrig2")
#print( xtable(result.tab), include.rownames=FALSE)

mean.tab
mean.tab[,1] = result.tab[,1]
#print( xtable(mean.tab), include.rownames=FALSE)


end_time = Sys.time()
end_time - start_time
#########################################################################
#########################################################################





