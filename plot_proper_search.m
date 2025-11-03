load("outputs/nos3_proper_search.mat")

lag = proper_search.lag;
hist_res_orthogonality = proper_search.hist_res_orthogonality;
hist_search_direction_A_orthogonality = proper_search.hist_search_direction_A_orthogonality;

x_axis = 0:metadata.iterations_performed;

figure;
hold on, grid on;
plot(x_axis,hist_res_orthogonality,'-*','DisplayName',sprintf('|(r_{k}, r_{k-%d})| / ||r_k||_2 ||r_{k-%d}||_2', lag, lag));
plot(x_axis,hist_search_direction_A_orthogonality,'-*','DisplayName',sprintf('|(p_{k}, Ap_{k-%d})| / ||p_k||_2 ||Ap_{k-%d}||_2', lag, lag));
legend, box on;
title(metadata.matrix_name);
xlabel('Number of Iterations');
ylabel('Log_{10} of orthogonality (cos)');
ylim(gca,[1e-19 1e+1]);
set(gca,...
    'FontSize',16,...
    'YScale','log',...
    'YTick',[1e-18 1e-16 1e-14 1e-12 1e-10 1e-08 1e-06 1e-04 1e-02 1e+0],...
    'YTickLabel',{'-18','-16','-14','-12','-10','-8','-6','-4','-2','0'});
hold off;
savefig(strcat('results/',metadata.matrix_name,'_orthogonality_hist.fig'));
saveas(gcf, strcat('results/',metadata.matrix_name,'_orthogonality_hist.pdf'));
