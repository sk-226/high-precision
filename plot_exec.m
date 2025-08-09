file_name = "results.mat";
load("outputs/"+file_name)
plot_conv_hist(data.convergence, data.metadata.matrix_name)
