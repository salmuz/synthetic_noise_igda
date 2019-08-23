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
gtest_data <- function(nlabels=2, p=2, data, mus, Sigmas, axes, 
                       eigvals, noise=1, tau=Inf, thick=0.2){
  # generate noise datatest 
  ntest.by.label <- as.integer(nrow(data)/nlabels)
  datatest <- matrix(0, ncol=p+1, nrow=nlabels*ntest.by.label)
  
  # max size confidence intervall ellipse 
  zelp <- qchisq(1-1e-2, p)
  col_idx_label <- p+1
  
  # generate noise data for root class center (0, 0)
  gnoise_data <- function(label, osign, mu, Sigma){
    # Generating testdata without any perturbation
    z <- subset(data, data[, col_idx_label] == label)[, -col_idx_label]
    # shift distance to closer or faraway
    mv.shift <- NULL
    for(i in 1:p) {
      .eigval <- sqrt(eigvals[label, i]*zelp) # middle of axis value
      min.eig <- if(noise==0) 0 else (noise)*.eigval
      max.eig <- if(noise==0) 0 else (noise+thick)*.eigval
      mv.shift <- cbind(mv.shift, rep(runif(1, min.eig, max.eig), ntest.by.label))
    }
    # noise dispersion wishart random covariance
    if(tau != Inf){
      S <- rWishart(n=1, df=tau*p, Sigma=Sigma)[,,1]
      mv.shift <- t(apply(mv.shift, MARGIN=1, FUN=mvrnorm, n=1, Sigma=S))
    }
    # puttin noise information to testing data
    z <- z + osign*mv.shift
    cbind(z, label)
  }
  naxes <- (1:nlabels)[-1]
  if(length(naxes) == 1){
    osign <- axes[replicate(ntest.by.label, naxes), ]    
  }else{
    osign <- axes[sample(naxes, ntest.by.label, replace = T), ]
  }
  # root label simulation
  datatest[1:ntest.by.label, ] <- gnoise_data(1, osign, mu=mus[1, ], Sigma=Sigmas[, , 1]) 
  
  # generate noise data for other class amongts center class
  for(i in 2:nlabels){
    translate.mu <- mus[1, ] - mus[i,]
    osign <- t(replicate(ntest.by.label, translate.mu/abs(translate.mu)))
    datatest[(ntest.by.label*(i-1)+1):(ntest.by.label*i), ] <- 
      gnoise_data(i, osign, mus[i, ], Sigmas[, , i])
  }
  datatest
}

# generate perfect separting sphere data
# d_sep_level distance separating level
# ntrain number training data by label
gperf_separating_sphere_data <- function(ndata.by.label = 200, p = 2, nlabels=2, 
                                         eigv.large = 1e-1, 
                                         d_sep_level=1-1e-2,
                                         homoscestic=T, 
                                         sphere=T,
                                         pct_test=0.2){
  g_covariance <- function(){
    H <- qr.Q(qr(matrix(rnorm(p^2), p)))
    eigvalues <- if(sphere) rep(runif(1, max=eigv.large), p) else runif(p, max=eigv.large)
    Sigma <- crossprod(H, H*eigvalues)  
    list('eigvals'=eigvalues, "Sigma"=Sigma)
  }
  
  # size confidence intervall ellipse 
  zelp <- qchisq(d_sep_level, p)
  
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
  x <- mvrnorm(ndata.by.label, mu=mu[1, ], Sigma = Sigma) 
  dataset <- matrix(0, ncol=p+1, nrow=nlabels*ndata.by.label)
  dataset[1:ndata.by.label, ] <- cbind(x, 1)
  for(i in 1:(nlabels-1)){
     osign <- axes[position[i],]
     lbaxes[i+1, ] <- osign
     if(!homoscestic) {
       gcov <- g_covariance()
       eigvals[i, ] <- gcov$eigvals 
       Sigma <- gcov$Sigma
       Sigmas[, , i+1] <- Sigma
       mv.shift <- 0.5*sqrt(eigvals.orig*zelp) + 0.5*sqrt(eigvals[i, ]*zelp)
     }else{
       mv.shift <- 2*sqrt(eigvals.orig*zelp) # 2-times diameter axes ellipse
     }
     tmu <- mu[1,] + osign*mv.shift
     mu[i+1, ] <- t(tmu)
     z <- mvrnorm(ndata.by.label, mu=tmu, Sigma = Sigma)
     dataset[(ndata.by.label*i+1):(ndata.by.label*(i+1)), ] <- cbind(z, i+1)
  }
  
  # split training data and testing data (randomly)
  col_idx_label <- p+1
  datatest <- matrix(0, ncol=p+1, nrow=nlabels*ndata.by.label*pct_test)
  ntest.by.label <- ndata.by.label*pct_test
  .test.idxs  <- NULL
  for(label in 1:nlabels){
   .label_idx <- which(dataset[,col_idx_label]==label)
   .ind_idx <- sample(.label_idx, ntest.by.label)
   datatest[(ntest.by.label*(label-1)+1):(ntest.by.label*label), ] <- dataset[.ind_idx,]
   .test.idxs <- c(.test.idxs, .ind_idx)
  }
  datatrain <- dataset[-.test.idxs,]
  
  # create and putting feacture and label names to datasets (training and testing)
  name_feactures <- paste('X', 1:p, sep='')
  name_label <- "y"
  col_names <- c(name_feactures, name_label)
  colnames(dataset) <- col_names
  colnames(datatest) <- col_names
  colnames(datatrain) <- col_names
  
  return(
    list(
      "data" = dataset,
      "train" = datatrain, 
      "test" = datatest,
      "mus" = mu,
      "Sigmas" = Sigmas,
      "axes" = lbaxes,
      "eigvals" = eigvals
    )
  )
}

# plot data in 2-dimension or 3-dimension 
plot_classification <- function(data, p){
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


# calculate eigen sample values
compute_eigenvalue_from_sample <- function(data){
  labels <- unique(data[, ncol(data)])
  for(label in labels){
    .temp <- which(data[, ncol(data)] == label)
    print(paste("Label(", label, "): ", 
                eigen(cov(data[.temp, -ncol(data)]))$values, sep = ""))
  }
}

# create file csv 
create_file_csv <- function(data, name){
  write.table(data, 
              sep=',',
              file=name, 
              col.names=F, 
              row.names = FALSE)
}

# get optimal ell parameter in average
get_ell_mean_optimal <- function(datafile='results_iris.csv') {
  read_data <- function(file_path){
    read.csv(file_path, sep=",", 
             header=FALSE, 
             colClasses = c('numeric', 'numeric', 'numeric', 'numeric'))
  }
  data <- read_data(datafile)
  idx_remove <- which(data[,1]== -9999)
  if(length(idx_remove) >0) data <- data[-idx_remove,]
  k_fold <- sum(data[,1]== -999)
  data_ells <- data[-which(data[,1]== -999), -2]
  data_ells[, 1] <- round(data_ells[, 1], 2)
  keys <- unique(data_ells[,1])
  .ellmeans <- unlist(lapply(keys, function(key) colMeans(data_ells[data_ells[,1]==key, ])))
  .ellsd <- lapply(keys, function(key) apply(data_ells[data_ells[,1]==key, -1], 2, sd))
  sd_data <- matrix(unlist(.ellsd), ncol = 2, byrow=T)
  mean_data <- matrix(.ellmeans, ncol=3, byrow = T)
  ell65_mean_opt <- mean_data[which(max(mean_data[, 2])==mean_data[, 2]), 1]
  ell80_mean_opt <- mean_data[which(max(mean_data[, 3])==mean_data[, 3]), 1]
  print(paste("ELL_65_MEAN_OPTIMAL:", ell65_mean_opt, sep=""))
  print(paste("ELL_80_MEAN_OPTIMAL:", ell80_mean_opt, sep=""))
  return(as.data.frame(mean_data))
}

