source("IGDA_synthetic_data_generation.R")
##################Generation TRAINING DATA#################
# generation 3-dataset - MIXING data to 50% confidence interval 
set.seed(29323)
dataset <- gperf_separating_sphere_data(ndata.by.label = 80, 
                                        p=3, nlabels = 5, 
                                        d_sep_level=0.35,
                                        eigv.large = 1e-1,
                                        sphere = F,
                                        homoscestic=F)
plot_classification(dataset$train, p=3)
plot_component_pca(dataset$train, component = c(1,2))
plot_component_pca(dataset$train, component = c(2,3))
plot_component_pca(dataset$train, component = c(1,3))
create_file_csv(dataset$train, "dsyn03_train.csv")
create_file_csv(dataset$data, "dsyn03_all_data.csv")

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (gamma)
n.test.corrupt <- 20
corrupt.test.gamma <- list()
noise.gamma  <- seq(0, 1, length.out = n.test.corrupt)
set.seed(6265)
for(i in 1:n.test.corrupt){
  .testing <- gtest_data(nlabels=5, p=3, 
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma[i])
  corrupt.test.gamma[[i]] <- .testing
}
#corrupt.test.gamma[[20]][, 4] <- 99
#plot_classification(rbind(dataset$train, corrupt.test.gamma[[20]]), p=3)
# create file test
for(i in 1:n.test.corrupt){
  create_file_csv(corrupt.test.gamma[[i]], paste('dsyn03_test_gamma_', i ,'.csv', sep = ""))
}

# (2) generation testing data with corruption data (tau)
corrupt.test.tau <- list()
noise.tau  <- 1:20
noise.gamma.fixed <- noise.gamma[5]
set.seed(9919)
for(i in 1:length(noise.tau)){
  .testing <- gtest_data(nlabels=5, p=3, 
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma.fixed,
                         tau=noise.tau[i])
  corrupt.test.tau[[i]] <- .testing
}
#corrupt.test.tau[[5]][, 4] <- 99
#plot_classification(rbind(dataset$train, corrupt.test.tau[[5]]), p=3)
# create file test
for(i in 1:length(noise.tau)){
  create_file_csv(corrupt.test.tau[[i]], paste('dsyn03_test_tau_', i ,'.csv', sep = ""))
}

##################Generation configuration#################
save(dataset, noise.gamma, corrupt.test.gamma, noise.tau, corrupt.test.tau, 
     file="dsyn03_conf_3x5_mix.RData")
rm(list = ls())





