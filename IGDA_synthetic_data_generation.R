# Generation of dataset synthetic
# IDEAS - Synthetic classification data
# (1) Perfect separation class in training data with testing corrupted data.
# (2) Linear separating decision boundary
# (3) Non-linear separating decision boundary with:
#  (3.1) Polinomial separating decision boundary
# (4) Two class each one with probability half of chances
# (5) Equal and different covariance matrix (homoscedastic and heteroscedastic)
library(MASS)
library(plotly)
library(ggplot2)
# (1) same covariance matrix 
# generate test noise sphere data
gtest_data <- function(ntest=50, nlabels=2, p=2, mus, Sigmas, axes, 
                       eigvals, noise=1, tau=Inf, thick=0.5){
  # generate noise datatest 
  datatest <- matrix(0, ncol=p+1, nrow=nlabels*ntest)
  
  # generate noise data for root class center (0, 0)
  gnoise_data <- function(label, osign, mu, Sigma){
    mv.shift <- NULL
    for(i in 1:p) {
      min.eig <-  (noise)*sqrt(eigvals[label, i])
      max.eig <- (noise+thick)*sqrt(eigvals[label, i])
      mv.shift <- cbind(mv.shift, 0.5*runif(ntest,min.eig, max.eig))
    }
    tmus <- t(replicate(ntest, mu)) + osign*mv.shift
    S <- if(tau == Inf) Sigma else rWishart(n=1, df=tau*p, Sigma=Sigma)[,,1]
    z <- t(apply(tmus, MARGIN=1, FUN=mvrnorm, n=1, Sigma=S))
    cbind(z, label)
  }
  naxes <- (1:nlabels)[-1]
  if(length(naxes) == 1){
    osign <- axes[replicate(ntest, naxes), ]    
  }else{
    osign <- axes[sample(naxes, ntest, replace = T), ]
  }
  datatest[1:ntest, ] <- gnoise_data(1, osign, mu=mus[1, ], Sigma=Sigmas[, , 1]) # root label simulation
  
  # generate noise data for other class amongts center class
  for(i in 2:nlabels){
    translate.mu <- mus[1, ] - mus[i,]
    osign <- t(replicate(ntest, translate.mu/abs(translate.mu)))
    datatest[(ntest*(i-1)+1):(ntest*i), ] <- gnoise_data(i, osign, mus[i, ], Sigmas[, , i])
  }
  datatest
}

# generate perfect separting sphere data
gperf_separating_sphere_data <- function(ntrain = 200, eigv.large = 1e-1, dseparating=4.1,
                                         p = 2, nlabels=2, homoscestic=T, sphere=T){
  g_covariance <- function(){
    H <- qr.Q(qr(matrix(rnorm(p^2), p)))
    eigvalues <- if(sphere) rep(runif(1, max=eigv.large), p) else runif(p, max=eigv.large)
    Sigma <- crossprod(H, H*eigvalues)  
    list('eigvals'=eigvalues, "Sigma"=Sigma)
  }
  
  # first dataset configuration
  gcov <- g_covariance()
  eigvals.orig <- gcov$eigvals 
  Sigma <- gcov$Sigma
  mu <- matrix(0, ncol=p, nrow=nlabels, byrow = T)
  
  # save information
  Sigmas <- array(Sigma, dim=c(p, p, nlabels))
  eigvals <- matrix(eigvals.orig, ncol = p, nrow=nlabels, byrow = T)
  axes <- as.matrix(expand.grid(rep(list(c(-1,1)), p)))
  position <- sample(1:nrow(axes), nlabels-1)
  lbaxes <- matrix(0, ncol=p, nrow=nrow(axes))
  
  # create dataset and save first root label
  x <- mvrnorm(ntrain, mu=mu[1, ], Sigma = Sigma) 
  dataset <- matrix(0, ncol=p+1, nrow=nlabels*ntrain)
  dataset[1:ntrain, ] <- cbind(x, 1)
  for(i in 1:(nlabels-1)){
     osign <- axes[position[i],]
     lbaxes[i+1, ] <- osign
     if(!homoscestic) {
       gcov <- g_covariance()
       eigvals <- gcov$eigvals 
       Sigma <- gcov$Sigma
       Sigmas[, , i+1] <- Sigma
       mv.shift <- (dseparating/2)*sqrt(eigvals.orig) + (dseparating/2)*sqrt(eigvals)
     }else{
       mv.shift <- dseparating*sqrt(eigvals.orig)
     }
     tmu <- mu[1,] + osign*mv.shift
     mu[i+1, ] <- t(tmu)
     z <- mvrnorm(ntrain, mu=tmu, Sigma = Sigma)
     dataset[(ntrain*i+1):(ntrain*(i+1)), ] <- cbind(z, i+1)
  }
  return(list("data"=dataset, "mus"=mu, "Sigmas"=Sigmas, "axes"=lbaxes, "eigvals"=eigvals))
}

# plot data in 2-dimension or 3-dimension 
plot_classification <- function(data, nbylabel, nlabels, p){
  g <- NULL
  if(p==2){
    colnames(data) <- c("x1", "x2", "y")
    data <- as.data.frame(data)
    data$y <- as.factor(data$y)
    g <- ggplot(data, aes(x=x1, y=x2)) + 
      geom_point(aes(shape=y, color=y, size = 1.5)) +
      labs( title="", x = 'X1', y = 'X2', color="Labels") +
      guides(shape = FALSE, size=FALSE) + 
      scale_shape_manual(values = c(1, 2, 3, 8, 13)) +
      theme_bw() +
      theme(legend.position= c(0.1, 0.15), 
            legend.text=element_text(size=rel(1.2)), 
            legend.background = element_blank(),
            axis.title.y=element_text(size=rel(1.1)),
            axis.title.x=element_text(size=rel(1.2)))
  }else if(p==3){
    colnames(data) <- c("x1", "x2", "x3", "y")
    data <- as.data.frame(data)
    data$y <- as.factor(data$y)
    g <- plot_ly(x=data$x1, y=data$x2, z=data$x3, color=data$y, symbol=data$y,
            symbols = c('circle-open-dot','square-open', 'diamond-open', 
                        'circle-open', 'x'), type="scatter3d", mode="markers")
    #newdata <- rbind(data, datatest)
    #plot_ly(x=newdata[, 1], y=newdata[, 2], z=newdata[, 3], color=newdata[, 4])
  }
  return(g)
}

# plot testing data noise information
plot_testing_data2d <- function(gplot, data){
  colnames(data) <- c("x1", "x2", "y")
  data <- as.data.frame(data)
  data$y <- as.factor(data$y)
  g <- gplot + geom_point(data=data, aes(x=x1, y=x2, shape=y, size = 1.5), color="black")
  g
}

##################Synthectic Datasets#################
# generation 1-dataset - SEPARATING MIXING
set.seed(708)
dataset.one <- gperf_separating_sphere_data(p=2, nlabels = 3, dseparating=2, ntrain = 40)
g <- plot_classification(dataset.one$data, nbylabel=Ntrain, nlabels=3, p=2)
# generation testing data with corruption data (gamma)
datatest <- list()
for(i in 1:4){
  set.seed(3020)
  .testing <- gtest_data(ntest=10, nlabels=3, p=2, mus=dataset.one$mus, 
                         Sigmas=dataset.one$Sigmas, axes=dataset.one$axes, 
                         eigvals=dataset.one$eigvals, noise=0)
  datatest[[i]] <- .testing
}
plot_testing_data2d(g, datatest[[1]])
plot_testing_data2d(g, datatest[[2]])
plot_testing_data2d(g, datatest[[3]])
plot_testing_data2d(g, datatest[[4]])

save(dataset.one, datatest, file="conf_data_synthetic_01.RData")

# generation 1-dataset - SEPARATING PERFECT 



# generation testing data with corruption data (tau)
set.seed(307)
datatest.one <- gtest_data(nlabels=3, p=2, mus=dataset.one$mus, Sigmas=dataset.one$Sigmas, 
                           axes=dataset.one$axes, eigvals=dataset.one$eigvals, noise=0.5, tau=10)
plot_testing_data2d(g, datatest.one)


# generation 2-dataset
set.seed(3279865)
dataset.two <- gperf_separating_sphere_data(p=2, nlabels = 3, eigv.large = 100, dseparating=2,
                                            ntrain = Ntrain, homoscestic = F)
plot_classification(dataset.two$data, nbylabel=Ntrain, nlabels=3, p=2)
# generation 3-dataset
set.seed(9269014)
dataset.three <- gperf_separating_sphere_data(p=3, nlabels = 5, eigv.large = 100, dseparating=2,
                                          ntrain = Ntrain, homoscestic = F)
plot_classification(dataset.two$data, nbylabel=Ntrain, nlabels=3, p=2)


# model 
data.train <- data.frame(x1=data[,1], x2=data[,2], y=data[,3])
lda.model <- lda(y~., data = data.train)
test.predict <- predict(lda.model, newdata=data.frame(x1=testx1[,1], x2=testx1[,2]))


# Empirical experiments for imprecise distribution marginal 
n = 10
n1 = 0.5
n2 = 0.2
c = seq(0, 2, length.out = 20)
maxx <- matrix(-1, nrow = 20, ncol = 4)
maxx[,2] <- -n2 - c/n
maxx[,3] <- n1 - c/n - 1
maxx[,4] <- n1 - c/n - n2 - c/n
plot(c, maxx[,1], type='l', ylim=c(min(maxx), max(maxx)))
lines(c, maxx[,2])
lines(c, maxx[,3], col="blue")
lines(c, maxx[,4], col='red')

p1 <- 0.5
p2 <- 0.5
maxx <- matrix(-p2, nrow = 20, ncol = 4)
maxx[,2] <- p2*(-n2 - c/n)
maxx[,3] <- p1*(n1 - c/n) - p2
maxx[,4] <- p1*(n1 - c/n) - p2*(n2 + c/n)
plot(c, maxx[,1], type='l', ylim=c(min(maxx), max(maxx)))
lines(c, maxx[,2])
lines(c, maxx[,3], col="blue")
lines(c, maxx[,4], col='red')

