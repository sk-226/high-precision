function plot_conv_hist(results, figure_title, opt)

    arguments
        results
        figure_title string
        opt.show_plot = true
        opt.save_fig = false
        opt.save_fig_filename = ""
        opt.format_type = "pdf"
        opt.output_dir = ""
        opt.y_lim = [1e-13 1e+1]
        opt.y_tick = [1e-12 1e-09 1e-06 1e-3 1]
        opt.y_tick_label = {'-12', '-9', '-6', '-3', '0'}
    end

    % FOR EXPERIMENTS
    if opt.show_plot
        fig = figure();
    else
        fig = figure('Visible', 'off');
    end

    % plot convergence history
    % CAUTION: x_axis starts from 0
    x_axis = 0:results.iter_final;
    hold on, grid on;
    plot(x_axis, results.hist_relres_2, '-*', 'DisplayName', strcat('||r_k||_2/||b||_2'));
    plot(x_axis, results.hist_relerr_2, '-*', 'DisplayName', strcat('||e_k||_2/||x_{true}||_2'));
    plot(x_axis, results.hist_relerr_A, '-*', 'DisplayName', strcat('||e_k||_A/||x_{true}||_A'));
    legend, box on;
    title(figure_title, 'Interpreter', 'none');
    xlabel('Number of Iterations');
    ylabel('Log_{10} of relative norm');
    ylim(gca, opt.y_lim);
    set(gca, ...
        'FontSize', 16, ...
        'YScale', 'log', ...
        'YTick', opt.y_tick, ...
        'YTickLabel', opt.y_tick_label);
    hold off;

    if opt.save_fig
        % Resolve nested directories included in save_fig_filename and ensure directory exists
        filename_norm = strrep(opt.save_fig_filename, '\\', filesep);
        [filename_dir, filename, ~] = fileparts(filename_norm);
        target_dir = opt.output_dir;
        if ~isempty(filename_dir)
            target_dir = fullfile(opt.output_dir, filename_dir);
        end
        if ~exist(target_dir, 'dir')
            mkdir(target_dir);
        end
        % Build complete file path including extension
        file_name = fullfile(target_dir, strcat(filename, ".", opt.format_type));
        saveas(fig, file_name);
    end

    % CLOSE FIGURE IF SHOW_PLOT IS FALSE
    if ~opt.show_plot
        close(fig);
    end

end
