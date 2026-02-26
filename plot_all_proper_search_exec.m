inputs_dir = 'outputs_dd_k_k-2';
output_dir = 'outputs_dd_proper_search_figs';

inputs = dir(fullfile(inputs_dir,'**','*.mat'));
inputs = inputs(~[inputs.isdir]); % inputs.name を使う
n_files = length(inputs);

for i = 1:n_files
    data = load(fullfile(inputs(i).folder, inputs(i).name));
    % Check if proper_search and metadata exist
    if isfield(data, 'proper_search') && isfield(data, 'metadata')
        plot_proper_search(data, data.metadata.matrix_name, ...
            'show_plot', false, ...
            'save_fig', true, ...
            'save_fig_filename', strcat(data.metadata.matrix_name, '_proper_search_lag', num2str(data.proper_search.lag)), ...
            'format_type', ["pdf", "fig"], ...
            'output_dir', output_dir)
    else
        warning('Skipping %s: missing proper_search or metadata field', inputs(i).name);
    end
end
