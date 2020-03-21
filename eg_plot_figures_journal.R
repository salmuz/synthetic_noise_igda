###### Figures for shifting results
root_path <- '.../dsyn02_2x5_mix/noise_mean/'
data <- read.csv(paste0(root_path, "results_ilda_vs_lda_noise_mean_dsyn02_accuracy.csv"), 
                 header = FALSE)
plot_figures_paper(data, is_noise_sigma = F)

###### Figures for noise dispersion results
root_path <- '.../dsyn02_2x5_mix/noise_sigma/'
data <- read.csv(paste0(root_path, "results_iqda_vs_qda_noise_sigma_dsyn02_accuracy.csv"), 
                 header = FALSE)
plot_figures_paper(data, is_noise_sigma = T)
