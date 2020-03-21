rm(list = ls())
source("IGDA_synthetic_data_generation.R")
##################Generation TRAINING DATA#################
is_save_csv <- FALSE
p <- 6
nlabels <- 8
# generation 4-dataset - MIXING data to 1% confidence interval 
dataset <- list()
ndata.by.label <- c(10, 25, 50)
mus <- NULL
Sigmas <- NULL
eigvals <- NULL
lbaxes <- NULL
set.seed(80397)
NRepetitions <- 100
for(k in 1:length(ndata.by.label)){
  idx.ndbl <- paste0(ndata.by.label[k])
  dataset[[idx.ndbl]] <- list()
  for(i in 1:NRepetitions){
    dataset[[idx.ndbl]][[i]] <- 
      gperf_separating_sphere_data(ndata.by.label = ndata.by.label[k], 
                                   p=6, nlabels = 8, 
                                   d_sep_level=0.01,
                                   sphere=F,
                                   homoscestic=F,
                                   eigv.large = 1e-1,
                                   pct_test=0,
                                   mus=mus, 
                                   Sigmas=Sigmas,
                                   eigvals=eigvals,
                                   lbaxes=lbaxes)
    if(is_save_csv) 
      create_file_csv(dataset[[idx.ndbl]][[i]]$data, 
                    paste0("dsyn04_all_data_", idx.ndbl, "_", i,".csv"))
    if(is.null(mus)){
      mus <- dataset[[idx.ndbl]][[i]]$mus
      Sigmas <- dataset[[idx.ndbl]][[i]]$Sigmas
      eigvals <- dataset[[idx.ndbl]][[i]]$eigvals
      lbaxes <- dataset[[idx.ndbl]][[i]]$axes
    }
  }
}
plot_component_pca(dataset[["25"]][[1]]$train, component = c(1,2))
plot_component_pca(dataset[["25"]][[20]]$train, component = c(1,2))
plot_component_pca(dataset[["25"]][[30]]$train, component = c(1,2))
plot_component_pca(dataset[["25"]][[40]]$train, component = c(1,2))

# dataset <- gperf_separating_sphere_data(ndata.by.label = 60, 
#                                         p=6, nlabels = 8, 
#                                         d_sep_level=0.1,
#                                         eigv.large = 1e-1,
#                                         sphere = F,
#                                         homoscestic=F)

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (mean)
epsilons <- seq(0, 1, by = 0.02)
corrupt.test.means <- list()
nb.test <- 10000
set.seed(7879603)
for(epsilon in epsilons){
  idx.epsilon <- paste0(length(corrupt.test.means)+1)
  corrupt.test.means[[idx.epsilon]]  <- list()
  .testing = gtest_data(nlabels=8, p=6, 
                        nb.test=nb.test,
                        mus=mus, 
                        Sigmas=Sigmas, 
                        epsilon_mean=epsilon)
  corrupt.test.means[[idx.epsilon]] <- .testing
}
# testing correcting testing
plot_component_pca(corrupt.test.means[["1"]], component = c(1,2))
plot_component_pca(corrupt.test.means[["20"]], component = c(1,2))
plot_component_pca(corrupt.test.means[["50"]], component = c(1,2))

# save testing data in files csv
for(idx.epsilon in names(corrupt.test.means)){
  datatest <- corrupt.test.means[[idx.epsilon]]
  create_file_csv(datatest, paste0('dsyn04_test_mean_', idx.epsilon, '.csv'))
}


# (2) generation testing data with corruption data (covariance matrix)
psis <- seq(0, 1, by = 0.02)
corrupt.test.sigmas <- list()
nb.test <- 10000
set.seed(7879603)
# calculating inertia total
sigma_noise <- array(data = 0, dim = c(p, p, nlabels))
for(idx in 1:nlabels){
  sigma_noise[, , idx] <- rWishart(1, p, 1e-1*diag(p))
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

# save testing data in files csv
for(idx.psi in names(corrupt.test.sigmas)){
  datatest <- corrupt.test.sigmas[[idx.psi]]
  create_file_csv(datatest, paste0('dsyn04_test_sigma_', idx.psi, '.csv'))
}

# testing correcting testing
plot_component_pca(corrupt.test.sigmas[["1"]], component = c(1,2))
plot_component_pca(corrupt.test.sigmas[["50"]], component = c(1,2))

################## Computing theoric risk bayes ##############
for(size in names(dataset)){
  rbayes <- 0
  for(idataset in dataset[[size]])
    rbayes <- rbayes + calculate_theoric_risk_bayes(idataset$data,
                                                    idataset$mus, 
                                                    idataset$Sigmas)$risk.bayes  
  print(paste("Risque bayes training data:", round(rbayes/NRepetitions, 2)))
}


for(iepsilon in names(corrupt.test.means)){
  idata <- corrupt.test.means[[iepsilon]] 
  rbayes <- calculate_theoric_risk_bayes(idata, mus, Sigmas)$risk.bayes 
  print(paste0("Risque bayes testing data (",epsilons[as.integer(iepsilon)], "): ",
               round(rbayes,2)))
}

for(ipsi in names(corrupt.test.sigmas)){
  idata <- corrupt.test.sigmas[[ipsi]] 
  rbayes <- calculate_theoric_risk_bayes(idata, mus, Sigmas)$risk.bayes 
  print(paste0("Risque bayes testing-sigma data (",
               psis[as.integer(ipsi)], "): ",
               round(rbayes,3)))
}

##################Generation configuration#################
save(dataset, corrupt.test.means, corrupt.test.sigmas, file="dsyn04_conf_6x8_mix.RData")
