rm(list = ls())
source(file='IGDA_synthetic_data_generation.R')
##################Generation TRAINING DATA#################
is_save_csv <- FALSE
dataset <- list()
ndata.by.label <- c(10, 50, 100, 25)
mus <- NULL
Sigmas <- NULL
eigvals <- NULL
lbaxes <- NULL
set.seed(3279865)
NRepetitions <- 100
nlabels <- 3
p <- 2
eigv.large <- 1
# generation 1-dataset - SEPARATING PERFECT to 99.9% confidence interval 
for(k in 1:length(ndata.by.label)){
  idx.ndbl <- paste0(ndata.by.label[k])
  dataset[[idx.ndbl]] <- list()
  for(i in 1:NRepetitions){
    dataset[[idx.ndbl]][[i]] <- 
      gperf_separating_sphere_data(ndata.by.label = ndata.by.label[k], 
                                   p=2, nlabels = 3, 
                                   d_sep_level=.80,
                                   eigv.large = eigv.large,
                                   pct_test=0,
                                   mus=mus, 
                                   Sigmas=Sigmas, 
                                   eigvals=eigvals,
                                   lbaxes=lbaxes)
    if(is_save_csv) 
      create_file_csv(dataset[[idx.ndbl]][[i]]$data, 
                    paste0("dsyn01_all_data_", idx.ndbl, "_", i,".csv"))
    
    if(is.null(mus)){
      mus <- dataset[[idx.ndbl]][[i]]$mus
      Sigmas <- dataset[[idx.ndbl]][[i]]$Sigmas
      eigvals <- dataset[[idx.ndbl]][[i]]$eigvals
      lbaxes <- dataset[[idx.ndbl]][[i]]$axes
    }
  }
}
gplot <- plot_classification(dataset[["25"]][[2]]$train, p=2)
compute_eigenvalue_from_sample(dataset[["50"]][[2]]$train)
size_obs <-"10"
grid.arrange(plot_classification(dataset[[size_obs]][[1]]$train), 
             plot_classification(dataset[[size_obs]][[2]]$train),
             plot_classification(dataset[[size_obs]][[3]]$train),
             ncol=3)

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (mean)
epsilons <- seq(0, 1, by = 0.02)
corrupt.test.means <- list()
nb.test <- 10000
set.seed(51760)
for(epsilon in epsilons){
  idx.epsilon <- paste0(length(corrupt.test.means)+1)
  corrupt.test.means[[idx.epsilon]]  <- list()
  .testing = gtest_data(nlabels=3, p=2, 
                        nb.test=nb.test,
                        mus=mus, 
                        Sigmas=Sigmas, 
                        epsilon_mean=epsilon)
  corrupt.test.means[[idx.epsilon]] <- .testing
}
# testing correcting testing
plot_testing_data2d(gplot, corrupt.test.means[["1"]])
plot_testing_data2d(gplot, corrupt.test.means[["20"]])
plot_testing_data2d(gplot, corrupt.test.means[["40"]])
plot_testing_data2d(gplot, corrupt.test.means[["51"]])

# save testing data in files csv
for(idx.epsilon in names(corrupt.test.means)){
  datatest <- corrupt.test.means[[idx.epsilon]]
  create_file_csv(datatest, paste0('dsyn01_test_mean_', idx.epsilon, '.csv'))
}

# (2) generation testing data with corruption data (covariance matrix)
psis <- seq(0, 1, by = 0.02)
corrupt.test.sigmas <- list()
nb.test <- 10000
set.seed(51760)
# calculating inertia total
sigma_noise <- array(data = 0, dim = c(p, p, nlabels))
for(idx in 1:nlabels){
  sigma_noise[, , idx] <- rWishart(1, p, eigv.large*diag(p))
}
# calculating noise data sigma
for(psi in psis){
  idx.psi <- paste0(length(corrupt.test.sigmas)+1)
  .testing = gtest_data(nlabels=nlabels, p=p, 
                        nb.test=nb.test,
                        mus=mus, 
                        Sigmas=Sigmas, 
                        epsilon_sigma=psi,
                        sigma_noise=sigma_noise)
  corrupt.test.sigmas[[idx.psi]] <- .testing
}

# testing correcting testing
plot_testing_data2d(gplot, corrupt.test.sigmas[["1"]])
plot_testing_data2d(gplot, corrupt.test.sigmas[["20"]])
plot_testing_data2d(gplot, corrupt.test.sigmas[["40"]])
plot_testing_data2d(gplot, corrupt.test.sigmas[["51"]])

# save testing data in files csv
for(idx.psi in names(corrupt.test.sigmas)){
  datatest <- corrupt.test.sigmas[[idx.psi]]
  create_file_csv(datatest, paste0('dsyn01_test_sigma_', idx.psi, '.csv'))
}

################## Computing theoric risk bayes ##############
for(size in names(dataset)){
  rbayes <- 0
  for(idataset in dataset[[size]])
    rbayes <- rbayes + calculate_theoric_risk_bayes(idataset$data,
                                                    idataset$mus, 
                                                    idataset$Sigmas)$risk.bayes  
  print(paste("Risque bayes training data:", round(rbayes/100, 3)))
}


for(iepsilon in names(corrupt.test.means)){
  idata <- corrupt.test.means[[iepsilon]] 
  rbayes <- calculate_theoric_risk_bayes(idata, mus, Sigmas)$risk.bayes 
  print(paste0("Risque bayes testing-mean data (",
               epsilons[as.integer(iepsilon)], "): ",
               round(rbayes,3)))
}

for(ipsi in names(corrupt.test.sigmas)){
  idata <- corrupt.test.sigmas[[ipsi]] 
  rbayes <- calculate_theoric_risk_bayes(idata, mus, Sigmas)$risk.bayes 
  print(paste0("Risque bayes testing-sigma data (",
               psis[as.integer(ipsi)], "): ",
               round(rbayes,3)))
}

##################Generation configuration#################
save(dataset, corrupt.test.means, corrupt.test.sigmas, file="dsyn01_conf_2x3_ps.RData")

