source("IGDA_synthetic_data_generation.R")
##################Generation TRAINING DATA#################
# generation 5-dataset - MIXING data 
set.seed(5460)
dataset <- gperf_separating_sphere_data(ndata.by.label = 75, 
                                        p=45, nlabels = 25, 
                                        d_sep_level=1e-55, 
                                        eigv.large = 1,
                                        sphere = F,
                                        homoscestic=F)
plot(dataset$train[, c(1,2)], col=dataset$train[, ncol(dataset$train)])

# plotting data 
# plot_component_pca(dataset$train, component = c(1,2))
# plot_component_pca(dataset$train, component = c(4,2))
# plot_component_pca(dataset$train, component = c(2,3))
# plot_component_pca(dataset$train, component = c(1,3))
# library(rgl)
# pp <- prcomp(dataset$train[, -ncol(dataset$train)])
# plot3d(pp$x[,1], pp$x[,2], pp$x[,3], "x", "y", "z", col=dataset$train[,ncol(dataset$train)])
create_file_csv(dataset$train, "dsyn05_train.csv")
create_file_csv(dataset$data, "dsyn05_all_data.csv")

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (gamma)
n.test.corrupt <- 20
corrupt.test.gamma <- list()
noise.gamma  <- seq(0, 1, length.out = n.test.corrupt)
set.seed(42721)
for(i in 1:n.test.corrupt){
  .testing <- gtest_data(nlabels=25, p=45, 
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma[i])
  corrupt.test.gamma[[i]] <- .testing
}

# create file test
for(i in 1:n.test.corrupt){
  create_file_csv(corrupt.test.gamma[[i]], paste('dsyn05_test_gamma_', i ,'.csv', sep = ""))
}

# (2) generation testing data with corruption data (tau)
corrupt.test.tau <- list()
noise.tau  <- 1:20
noise.gamma.fixed <- noise.gamma[5]
set.seed(35973)
for(i in 1:length(noise.tau)){
  .testing <- gtest_data(nlabels=25, p=45,
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma.fixed,
                         tau=noise.tau[i])
  corrupt.test.tau[[i]] <- .testing
}
# create file test
for(i in 1:length(noise.tau)){
  create_file_csv(corrupt.test.tau[[i]], paste('dsyn05_test_tau_', i ,'.csv', sep = ""))
}

################## Computing theoric risk bayes ##############
rs <- calculate_theoric_risk_bayes(dataset$data, dataset$mus, dataset$Sigmas)
print(paste("Risque bayes all data:", round(rs$risk.bayes, 2)))
rs <- calculate_theoric_risk_bayes(dataset$test, dataset$mus, dataset$Sigmas)
print(paste("Risque bayes testing data:", round(rs$risk.bayes, 2)))
rs <- calculate_theoric_risk_bayes(dataset$train, dataset$mus, dataset$Sigmas)
print(paste("Risque bayes training data:", round(rs$risk.bayes, 2)))

##################Generation configuration#################
save(dataset, noise.gamma, corrupt.test.gamma, noise.tau, corrupt.test.tau, 
     file="dsyn05_conf_45x25_mix.RData")
rm(list = ls())






