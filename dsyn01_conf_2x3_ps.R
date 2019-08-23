source(file='IGDA_synthetic_data_generation.R')
##################Generation TRAINING DATA#################
# generation 1-dataset - SEPARATING PERFECT to 99.9% confidence interval 
set.seed(3279865)
dataset <- gperf_separating_sphere_data(ndata.by.label = 50, p=2, nlabels = 3, 
                                            d_sep_level=.80,
                                            eigv.large = 1)
gplot <- plot_classification(dataset$train, p=2)
compute_eigenvalue_from_sample(dataset$train)
print(dataset$eigvals)
plot_testing_data2d(gplot, dataset$test)
create_file_csv(dataset$train, "dsyn01_train.csv")

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (gamma)
n.test.corrupt <- 20
corrupt.test.gamma <- list()
noise.gamma  <- seq(0, 1, length.out = n.test.corrupt)
set.seed(51760)
for(i in 1:n.test.corrupt){
  .testing <- gtest_data(nlabels=3, p=2, 
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma[i])
  corrupt.test.gamma[[i]] <- .testing
}

# verification shift test-dataset
plot_testing_data2d(gplot, corrupt.test.gamma[[1]])
plot_testing_data2d(gplot, corrupt.test.gamma[[5]])
plot_testing_data2d(gplot, corrupt.test.gamma[[10]])
plot_testing_data2d(gplot, corrupt.test.gamma[[20]])
# create file test
for(i in 1:n.test.corrupt){
  create_file_csv(corrupt.test.gamma[[i]], paste('dsyn01_test_gamma_', i ,'.csv', sep = ""))
}

# (2) generation testing data with corruption data (tau)
corrupt.test.tau <- list()
noise.tau  <- 1:10
noise.gamma.fixed <- noise.gamma[5]
set.seed(9919)
for(i in 1:length(noise.tau)){
  .testing <- gtest_data(nlabels=3, p=2, 
                         data=dataset$test, 
                         mus=dataset$mus, 
                         Sigmas=dataset$Sigmas, 
                         axes=dataset$axes, 
                         eigvals=dataset$eigvals, 
                         noise=noise.gamma.fixed,
                         tau=noise.tau[i])
  corrupt.test.tau[[i]] <- .testing
}
# verification shift test-dataset
plot_testing_data2d(gplot, corrupt.test.tau[[1]])
plot_testing_data2d(gplot, corrupt.test.tau[[3]])
plot_testing_data2d(gplot, corrupt.test.tau[[4]])
plot_testing_data2d(gplot, corrupt.test.tau[[7]])
plot_testing_data2d(gplot, corrupt.test.tau[[10]])
# create file test
for(i in 1:length(noise.tau)){
  create_file_csv(corrupt.test.tau[[i]], paste('dsyn01_test_tau_', i ,'.csv', sep = ""))
}

##################Generation configuration#################
save(dataset, noise.gamma, corrupt.test.gamma, noise.tau, corrupt.test.tau, 
     file="dsyn01_conf_2x3_mix.RData")
rm(list = ls())

