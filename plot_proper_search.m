function plot_proper_search(results, figure_title, opt)

    arguments
        results
        figure_title string
        opt.show_plot = true
        opt.save_fig = false
        opt.save_fig_filename = ""
        opt.format_type = ["pdf", "fig"]
        opt.output_dir = ""
        opt.y_lim = []
        opt.y_tick = []
        opt.y_tick_label = {}
    end

    % FOR EXPERIMENTS
    if opt.show_plot
        fig = figure();
    else
        fig = figure('Visible', 'off');
    end

    % plot proper search history
    % CAUTION: x_axis starts from 0
    lag = results.proper_search.lag;
    x_axis = 0:(numel(results.proper_search.hist_res_orthogonality) - 1);
    hold on, grid on;
    plot(x_axis, results.proper_search.hist_res_orthogonality, '-*', 'DisplayName', sprintf('|(r_{k}, r_{k-%d})| / ||r_k||_2 ||r_{k-%d}||_2', lag, lag));
    plot(x_axis, results.proper_search.hist_search_direction_A_orthogonality, '-*', 'DisplayName', sprintf('|(p_{k}, Ap_{k-%d})| / ||p_k||_2 ||Ap_{k-%d}||_2', lag, lag));
    legend, box on;
    title(figure_title, 'Interpreter', 'none');
    xlabel('Number of Iterations');
    ylabel('Log_{10} of orthogonality (cos)');
    if ~isempty(opt.y_lim)
        ylim(gca, opt.y_lim);
    end
    if ~isempty(opt.y_tick) && ~isempty(opt.y_tick_label)
        set(gca, ...
            'FontSize', 16, ...
            'YScale', 'log', ...
            'YTick', opt.y_tick, ...
            'YTickLabel', opt.y_tick_label);
    else
        set(gca, ...
            'FontSize', 16, ...
            'YScale', 'log');
    end
    
    hold off;

    if opt.save_fig
        % Handle both string array and scalar string
        if isscalar(opt.format_type)
            format_types = {char(opt.format_type)};
        else
            format_types = cell(size(opt.format_type));
            for j = 1:length(opt.format_type)
                format_types{j} = char(opt.format_type(j));
            end
        end
        for i = 1:length(format_types)
            format_type = format_types{i};
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
            file_name = fullfile(target_dir, strcat(filename, ".", format_type));
            saveas(fig, file_name);
        end
    end

    % CLOSE FIGURE IF SHOW_PLOT IS FALSE
    if ~opt.show_plot
        close(fig);
    end

end
