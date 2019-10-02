source(file='IGDA_synthetic_data_generation.R')
##################Generation TRAINING DATA#################
# generation 2-dataset - MIXING data to 95% confidence interval 
set.seed(708)
dataset <- gperf_separating_sphere_data(ndata.by.label = 50, p=2, nlabels = 5, 
                                        d_sep_level=0.6,
                                        sphere=F,
                                        homoscestic=F,
                                        eigv.large = 1e-1)
gplot <- plot_classification(dataset$train, p=2)
compute_eigenvalue_from_sample(dataset$train)
print(dataset$eigvals)
plot_testing_data2d(gplot, dataset$test)
create_file_csv(dataset$train, "dsyn02_train.csv")
create_file_csv(dataset$data, "dsyn02_all_data.csv")

##################Generation CORRUPT TEST DATA#################
# (1) generation testing data with corruption data (gamma)
n.test.corrupt <- 20
corrupt.test.gamma <- list()
noise.gamma  <- seq(0, 1, length.out = n.test.corrupt)
set.seed(3020)
for(i in 1:n.test.corrupt){
  .testing <- gtest_data(nlabels=5, p=2, 
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
  create_file_csv(corrupt.test.gamma[[i]], paste('dsyn02_test_gamma_', i ,'.csv', sep = ""))
}

# (2) generation testing data with corruption data (tau)
corrupt.test.tau <- list()
noise.tau  <- 1:20
noise.gamma.fixed <- noise.gamma[5]
set.seed(3420)
for(i in 1:length(noise.tau)){
  .testing <- gtest_data(nlabels=5, p=2, 
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
  create_file_csv(corrupt.test.tau[[i]], paste('dsyn02_test_tau_', i ,'.csv', sep = ""))
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
     file="dsyn02_conf_2x5_mix.RData")
rm(list = ls())
