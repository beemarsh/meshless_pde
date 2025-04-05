function save_matrix(matrix, nodes, shapes, filename)
% Create the results directory if it doesn't exist
if ~exist('results', 'dir')
    mkdir('results');
end

% Create the full Excel filename with path
excel_file = ['results/' filename '.xlsx'];

% Delete the file if it already exists
if exist(excel_file, 'file')
    delete(excel_file);
end

% Current row position tracker
current_row = 1;

% For each node configuration
for i = 1:length(nodes)
    % Prepare header row for nodes
    header_text = {sprintf('Number of Nodes = %d (Grid Size = %d)', (nodes(i)+1)^2, nodes(i)+1)};
    writecell(header_text, excel_file, 'Sheet', 1, 'Range', sprintf('A%d', current_row));

    % Create shape parameter headers
    shape_headers = cell(1, length(shapes) + 1);
    shape_headers{1} = 'Eigenvalue';
    for j = 1:length(shapes)
        shape_headers{j+1} = sprintf('SP = %.3e', shapes(j));
    end

    % Write shape headers
    writecell(shape_headers, excel_file, 'Sheet', 1, ...
        'Range', sprintf('A%d', current_row + 1));

    % Write matrix data for this node
    data_slice = squeeze(matrix(:,:,i));
    % Add eigenvalue numbers to first column
    data_with_index = [((1:size(data_slice,1))'),...
        data_slice];

    % Write the data
    writematrix(data_with_index, excel_file, 'Sheet', 1, ...
        'Range', sprintf('A%d', current_row + 2));

    % Update current row position:
    % 1 row for node header
    % 1 row for shape parameters header
    % size(data_slice,1) rows for data
    % 3 rows for spacing
    current_row = current_row + 2 + size(data_slice,1) + 3;
end
end