file_name = "nos7_dq.mat";
load("outputs_dq/HB/"+file_name)

plot_conv_hist(convergence, metadata.matrix_name)
