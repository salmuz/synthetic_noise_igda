# Generation of dataset synthetic
# IDEAS - Synthetic classification data
# (1) Perfect separation class in training data with testing corrupted data.
# (2) Linear/non-linear separating decision boundary
# (3) Equal and different covariance matrix (homoscedastic and heteroscedastic)
library(MASS)
library(plotly)
library(ggplot2)
library(gridExtra)

# generate new test noise from a pertubation of original distribution
gtest_data <- function(nlabels=2, p=2, nb.test, mus, Sigmas, epsilon_mean=0, 
                       epsilon_sigma=0, sigma_noise=NULL){
  if(is.null(sigma_noise))
    sigma_noise = array(data = 0, dim = c(p, p, nlabels))
  # generate noise datatest 
  nb.test.by.label <- as.integer(nb.test/nlabels)
  datatest <- matrix(0, ncol=p+1, nrow=nlabels*nb.test.by.label)
  
  # modify the means of each sub-population
  muc <- colMeans(mus)
  mus_update <- (1-epsilon_mean)*mus + 
    epsilon_mean*matrix(rep(muc, nrow(mus)), ncol=p, byrow = T) 
    
  # generation de new sampling 
  for(i in 1:nlabels){
    Sigmas[,,i]  <- (1-epsilon_sigma)*Sigmas[,,i] + epsilon_sigma*sigma_noise[,,i]
    x.new <- mvrnorm(nb.test.by.label, mu = mus_update[i, ], Sigma = Sigmas[,,i])
    datatest[(nb.test.by.label*(i-1)+1):(nb.test.by.label*i), ] <- cbind(x.new, i)
  }
  datatest
}

# generate new test noise from a pertubation of original testing data
gperturbation_test_data <- function(nlabels=2, p=2, data, mus, Sigmas, axes, 
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
      S <- rWishart(n=1, df=p + tau, Sigma=Sigma)[,,1]
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
gperf_separating_sphere_data <- function(ndata.by.label = 200, 
                                         p = 2, nlabels=2, 
                                         eigv.large = 1e-1, 
                                         d_sep_level=1-1e-2,
                                         homoscestic=T, 
                                         sphere=T,
                                         pct_test=0.2,
                                         mus=NULL, 
                                         Sigmas=NULL, 
                                         eigvals=NULL,
                                         lbaxes=NULL){
  
  if(xor(is.null(mus), is.null(Sigmas)) || 
     xor(is.null(eigvals), is.null(Sigmas)) ||
     xor(is.null(eigvals), is.null(mus))){
    stop("Means and Matrix-Covariance should be supplied together.")  
  }
  
  g_covariance <- function(){
    H <- qr.Q(qr(matrix(rnorm(p^2), p)))
    eigvalues <- if(sphere) rep(runif(1, max=eigv.large), p) 
        else sort(runif(p, max=eigv.large), decreasing = T)
    Sigma <- crossprod(H, H*eigvalues)  
    list('eigvals'=eigvalues, "Sigma"=Sigma)
  }
  
  # size confidence intervall ellipse if covariance is known 
  # H_0: mu=mu_0 vs H_1: mu/=mu_0 
  # W = n(bx - mu_0)'SigmaInv(bx - mu_0) ~ ChisQ 
  zelp <- qchisq(d_sep_level, p)
  
  # first dataset configuration
  if(is.null(Sigmas)){
    gcov <- g_covariance()
    eigvals.orig <- gcov$eigvals 
    Sigma <- gcov$Sigma
  }else{
    eigvals.orig <- eigvals[1, ]
  }
  
  # save information
  is.rec.means.cov <- FALSE
  if(is.null(mus) && is.null(Sigmas)){
    mus <- matrix(0, ncol=p, nrow=nlabels, byrow = T)
    Sigmas <- array(Sigma, dim=c(p, p, nlabels))
    eigvals <- matrix(eigvals.orig, ncol = p, nrow=nlabels, byrow = T)
    axes <- if(p <= 24) as.matrix(expand.grid(rep(list(c(-1,1)), p))) 
    else matrix(replicate(nlabels, replicate(p, sample(c(-1, 1), 1))), ncol=p, byrow = T)
    position <- sample(1:nrow(axes), nlabels-1)
    lbaxes <- matrix(0, ncol=p, nrow=nrow(axes)+1)
  }else{
    is.rec.means.cov <- TRUE  
  }
  
  # create dataset and save first root label
  x <- mvrnorm(ndata.by.label, mu=mus[1, ], Sigma = Sigmas[,,1]) 
  dataset <- matrix(0, ncol=p+1, nrow=nlabels*ndata.by.label)
  dataset[1:ndata.by.label, ] <- cbind(x, 1)
  for(i in 1:(nlabels - 1)) {
    if (!homoscestic) {
      if (!is.rec.means.cov) {
        gcov <- g_covariance()
        eigvals[i + 1,] <- gcov$eigvals
        Sigma <- gcov$Sigma
        Sigmas[, , i + 1] <- Sigma
      }
      mv.shift <- 0.8 * sqrt(eigvals.orig * zelp) + 0.8 * sqrt(eigvals[i + 1,] * zelp)
    } else{
      mv.shift <- 2 * sqrt(eigvals.orig * zelp) # 2-times diameter axes sphere
    }
    # recovery the position in the p-dimension
    if (!is.rec.means.cov) {
      osign <- axes[position[i], ]
      lbaxes[i + 1,] <- osign
    }else{
      osign <- lbaxes[i + 1,]
    }
    # distance shift of main population
    tmu <- mus[1, ] + osign * mv.shift
    if (!is.rec.means.cov) mus[i + 1,] <- t(tmu)
    # generation of sub-population
    z <- mvrnorm(ndata.by.label, mu = tmu, Sigma = Sigmas[, , i])
    dataset[(ndata.by.label * i + 1):(ndata.by.label * (i + 1)),] <- cbind(z, i + 1)
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
  if(length(.test.idxs) > 0)
    datatrain <- dataset[-.test.idxs,]
  else
    datatrain <- dataset
  
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
      "mus" = mus,
      "Sigmas" = Sigmas,
      "axes" = lbaxes,
      "eigvals" = eigvals
    )
  )
}

# plot data in 2-dimension or 3-dimension 
plot_classification <- function(data, p=NULL, 
                                level=0.8, 
                                covariances=NULL,
                                means=NULL,
                                ggpellipse=F, 
                                def_size=4, 
                                ylegend=0.2){
  g <- NULL
  p <- ifelse(is.null(p), ncol(data)-1, p)
  if(p==2){
    colnames(data) <- c("x1", "x2", "y")
    data <- as.data.frame(data)
    data$y <- as.factor(data$y)
    g <- ggplot(data, aes(x=x1, y=x2)) + 
      geom_point(aes(shape=y, color=y, size=1.5), size = def_size) +
      labs(title="", x = 'X1', y = 'X2', color="Labels") +
      guides(shape = FALSE, size=FALSE) + 
      scale_shape_manual(values = c(49, 50, 51, 52, 53, 54, 55)) +
      theme_bw() +
      theme(legend.position= c(0.1, ylegend), 
            legend.text=element_text(size=rel(1.5)), 
            legend.background = element_blank(),
            axis.title.y=element_text(size=rel(1.5)),
            axis.title.x=element_text(size=rel(1.5)))
    if(ggpellipse)
      g <- g + stat_ellipse(aes(color=y), type = "norm", linetype = 2, level=level)
    # Create confidence region with true covariance values 
    if(!is.null(covariances) && !is.null(means)){
      ellipse_data <- NULL
      for(i in 1:length(levels(data$y))){
        .t <- ellipse(covariances[, ,i], level=level, centre = means[i,])
        ellipse_data <- rbind(ellipse_data, cbind(.t, i))
      }
      ellipse_data <- as.data.frame(ellipse_data)
      colnames(ellipse_data) <- c("xellip", "yellip", "y")
      ellipse_data$y <- as.factor(ellipse_data$y)
      g <- g + geom_path(data=ellipse_data, aes(x=xellip, y=yellip, color=y), linetype = 2)
    }
  }else if(p==3){
    colnames(data) <- c("x1", "x2", "x3", "y")
    data <- as.data.frame(data)
    data$y <- as.factor(data$y)
    g <- plot_ly(x=data$x1, y=data$x2, z=data$x3, color=data$y, symbol=data$y,
            symbols = c('circle-open-dot','square-open', 'diamond-open', 
                        'circle-open', 'x'), type="scatter3d", mode="markers")
  }
  return(g)
}

# plot testing data noise information
plot_testing_data2d <- function(gplot, data, def_size = 4){
  colnames(data) <- c("x1", "x2", "y")
  data <- as.data.frame(data)
  data$y <- as.factor(data$y)
  g <- gplot + geom_point(data=data, aes(x=x1, y=x2, shape=y, size = 1.5),
                          size = def_size, color="black")
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
  keys <- sort(unique(data_ells[,1]))
  .ellmeans <- unlist(lapply(keys, function(key) colMeans(data_ells[data_ells[,1]==key, ])))
  .ellsd <- lapply(keys, function(key) apply(data_ells[data_ells[,1]==key, -1], 2, sd))
  sd_data <- matrix(unlist(.ellsd), ncol = 2, byrow=T)
  mean_data <- matrix(.ellmeans, ncol=3, byrow = T)
  mean_data[mean_data[,2] > 1, 2] <- 1
  mean_data[mean_data[,3] > 1, 3] <- 1
  ell65_mean_opt <- mean_data[which(max(mean_data[, 2])==mean_data[, 2]), 1]
  ell80_mean_opt <- mean_data[which(max(mean_data[, 3])==mean_data[, 3]), 1]
  print(paste("ELL_65_MEAN_OPTIMAL:", ell65_mean_opt, sep=""))
  print(paste("ELL_80_MEAN_OPTIMAL:", ell80_mean_opt, sep=""))
  return(as.data.frame(mean_data))
}

# plot component of principal component analysis
plot_component_pca <- function(data, component=c(1,2), label=-1){
  label <- ifelse(label == -1, ncol(data), label)
  pp <- prcomp(data[, -label])
  print(paste("Info. components:", cumsum(pp$sdev^2 / sum(pp$sdev^2)), sep=""))
  plot(pp$x[, component], col=data[, label])
}

# plot boxplot with 10-fold accuracy and percentage testing data  
plot_pct_testing_data_comparative <- function(root, model, leg.pos.def=c(0.2, 0.1)){
  model_pre <- model
  model_imp <- paste("i", model_pre, sep = "")
  root_imprecise <- paste(root, "/", model_imp, "/", sep="")
  root_precise <- paste(root, "_precise/", model_pre, "/", sep="")
  
  paths_imprecise <- dir(root_imprecise)
  paths_precise <- dir(root_precise)
  
  read_data <- function(root_path, file_name, ncol=3){
    read.csv(paste(root_path,file_name,sep=""), sep=",", 
             header=FALSE, colClasses = rep('numeric', ncol))
  }
  
  .get_mean_u65_u80 <- function(root, in_file, type="ilda", pct){
    .d <- read_data(root, in_file, 4)
    .d <- .d[which(.d[,1] == -999), c(3, 4)]
    .d <- cbind(type, as.integer(pct), .d)
    .d
  }
  
  allpct <- NULL
  for(i in 1:length(paths_imprecise)){
    .len <- nchar(paths_imprecise[i]) - 8
    pct <- substr(paths_imprecise[i], .len, .len+1)
    path_imprecise <- paths_imprecise[grep(pct, paths_imprecise)]
    path_precise <- paths_precise[grep(pct, paths_precise)]
    print(paste(path_imprecise, path_precise, sep=" vs "))
    .di <- .get_mean_u65_u80(root_imprecise, path_imprecise, model_imp, pct)
    .dp <- .get_mean_u65_u80(root_precise, path_precise, model_pre, pct)
    allpct <- rbind(allpct, .di)
    allpct <- rbind(allpct, .dp)
  }
  colnames(allpct) <- c("model", "pct", "u65", "u80")
  allpct[allpct$u65 > 1, 3] <- 1
  allpct[allpct$u80 > 1, 4] <- 1
  allpct$pct <- as.factor(allpct$pct)
  
  ggplot(as.data.frame(allpct), aes(x=pct, y=u80, fill=model)) + 
    geom_boxplot(alpha=0.7, outlier.colour="red", outlier.size=2.5)+
    labs(title="", x = 'Percentage of testing data', 
         y = 'Utility-discounted accuracy', fill="Models") +
    theme_bw() +
    theme(legend.position=leg.pos.def, legend.direction="horizontal",
          legend.text =element_text(size=rel(2.5)), 
          legend.title = element_text(size=rel(2.5)), 
          legend.background = element_blank(),
          axis.text.x = element_text(size=rel(2)),
          axis.text.y = element_text(size=rel(2)),
          axis.title.y=element_text(size=rel(2)),
          axis.title.x=element_text(size=rel(2)))
}

# plotting noise evolution with nrow data 40
plot_noise_evolution <- function(data, u80=TRUE, method='QDA'){
  default.mar <- par("mar")
  idx_u <- ifelse(u80, 4, 3)
  if(nrow(data) != 40) stop("Not 40 values for plotting")
  imethod <- paste('I', method, sep="")
  noise.gamma  <- seq(0, 1, length.out = 20)
  noise.tau  <- 1:20
  par(mar=c(4, 4, 1, 1))
  plot(noise.gamma, data[1:20, idx_u], type='b', col='red', lty=1, pch=2, 
       yaxt="n", xaxt="n", mgp = c(2.3, 1, 0),
       xlab='Noise parameter gamma', 
       ylab='discount-utility', cex.lab=2)
  lines(noise.gamma, data[1:20, 5], type='b', col='black', lty=3, pch=4)
  legend("topright", inset=.02, legend = c(imethod, method), 
         col=c("red", "black"), lty=c(1,3), pch=c(2,4), cex=2.5)
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5)
  
  plot(noise.tau, data[21:40, idx_u], type='b', col='red', yaxt="n", xaxt="n", mgp = c(2.3, 1, 0), cex.lab=2,
       xlab='Noise parameter tau', ylab='discount-utility', lty=1)
  lines(noise.tau, data[21:40, 5], type='b', col='black', lty=3, pch=4)
  legend("topright", inset=.02, legend = c(imethod, method), 
         col=c("red", "black"), lty=c(1,3), pch=c(2,4), cex=2.5)
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5)
  par(mar=default.mar)
}

# Computing the theoric risk bayes 
calculate_theoric_risk_bayes <- function(data, mus, Sigmas){
  if(!require("mvtnorm")) stop("If necessary to install mvtnorm package.")
  p <- ncol(data) - 1
  X <- data[, 1:p]
  y.obs <- data[, p+1]
  nb.labels <- nrow(mus)
  nb.inst <- nrow(X)
  post.prob <- matrix(0, nrow=nb.inst , ncol=nb.labels)
  marg.prob <- c(table(y.obs)/nb.inst) # pi_i
  for (idxlabel in 1:nb.labels){
    post.prob[, idxlabel] <- marg.prob[idxlabel] * 
      dmvnorm(X, mean=mus[idxlabel,], sigma = Sigmas[,,idxlabel]) 
  }
  y.pred <- max.col(post.prob)
  risk.bayes <- 1 - sum(diag(table(y.obs, y.pred)))/nb.inst
  return(list(
    "risk.bayes" = risk.bayes,
    "posterior.prob" = post.prob,
    "y.predition" = y.pred 
  ))
}

# print paper report
recovery_subdata  <- function(data, size=50, start_idx=17){
  list_name <- unlist(data[,1])
  idx_sub_data <- c(NULL)
  for(i in 1:length(list_name)){
    name <- list_name[i]
    if (as.integer(unlist(strsplit(
      substr(name, start_idx, nchar(as.character(name)) - 4), "_"))[2]) <= size){
      idx_sub_data <- c(idx_sub_data, i)
    }
  }
  data[idx_sub_data,]
}

plot_gg <- function(data, .title="", epsilons, start_idx=19){
  cdata <- ddply(data, "V2",  summarise,
                 u65_m = mean(V4),
                 u80_m = mean(V5),
                 set_m = mean(V6),
                 acc_m = mean(V7),
                 N = length(V4),
                 u65_se = sd(V4) / sqrt(length(V4)),
                 u80_se = sd(V5) / sqrt(length(V5)),
                 set_se = sd(V6) / sqrt(length(V6)),
                 acc_se = sd(V7) / sqrt(length(V7)),
                 u65_inf = u65_m - 1.96*u65_se,
                 u65_sup = u65_m + 1.96*u65_se,
                 u80_inf = u80_m - 1.96*u80_se,
                 u80_sup = u80_m + 1.96*u80_se,
                 set_inf = set_m - 1.96*set_se,
                 set_sup = set_m + 1.96*set_se,
                 acc_inf = acc_m - 1.96*acc_se,
                 acc_sup = acc_m + 1.96*acc_se
  )
  cdata[,1] <- as.integer(sapply(cdata[,1], 
                                 function(x) substr(x, start_idx, nchar(as.character(x))-4)))
  cdata <- cdata[ order(cdata[,1]),]
  cdata[, 1] <- epsilons[cdata[, 1]]
  g <- ggplot(cdata, aes(x=V2)) +
    geom_line(aes(y = u65_m, linetype="u65"), size=1.5) +
    geom_line(aes(y = u80_m, linetype="u80"), size=1.5) +
    geom_line(aes(y = acc_m, linetype="acc"), size=1.5) +
    geom_line(aes(y = set_m, linetype="set"), size=1.5) +
    geom_ribbon(aes(ymin=u65_inf, ymax=u65_sup, fill="u65"), alpha=0.3) +
    geom_ribbon(aes(ymin=u80_inf, ymax=u80_sup, fill="u80"), alpha=0.3) +
    geom_ribbon(aes(ymin=acc_inf, ymax=acc_sup, fill="acc"), alpha=0.3) +
    geom_ribbon(aes(ymin=set_inf, ymax=set_sup, fill="set"), alpha=0.3) +
    labs(fill="", linetype="", 
         axis.text=element_text(size=14), 
         axis.title=element_text(size=14,face="bold"),
         x = 'noise parameter', y = 'utility accuracy')+
    theme_bw() +
    scale_fill_discrete(guide=FALSE) +
    theme(legend.position = c(.95, .85), 
          legend.background = element_blank(), 
          legend.direction="vertical",
          legend.justification = c("right", "top"),
          legend.text=element_text(size=rel(3)),
          legend.title=element_text(size=rel(3)),
          axis.title.x=element_text(size=rel(3)),
          axis.title.y=element_text(size=rel(3)),
          axis.text=element_text(size=rel(2.5), color = "black"),
          legend.key.width = unit(3, "cm"),
          panel.border = element_rect(colour = "black", size = 2, linetype = "solid"))
  g
}

plot_all_data <- function(data, size, epsilons, start_idx=19){
  zepsilons <- length(epsilons)
  g1 <- plot_gg(data[1:(size*zepsilons), ], 
                .title="10", epsilons, start_idx)
  g2 <- plot_gg(data[(size*zepsilons+1):(size*zepsilons*2), ], 
                .title="25", epsilons, start_idx)
  g3 <- plot_gg(data[(size*zepsilons*2+1):(size*zepsilons*3), ], 
                .title="50", epsilons, start_idx)
  g4 <- plot_gg(data[(size*zepsilons*3+1):(size*zepsilons*4), ], 
                .title="100", epsilons, start_idx)
  print(g1)
  print(g2)
  print(g3)
  print(g4)
}
