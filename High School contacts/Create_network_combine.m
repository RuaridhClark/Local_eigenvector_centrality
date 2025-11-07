addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Local_eigenvector_centrality'

function [A_sparse, A_dense, nodeIDs_nonzero] = csv_to_adjacency(filename)
    %% --- Read Thiers 2011 Contact List ---
    
    % Read the tab-separated file: columns are [t i j Ci Cj]
    contactData = readtable(filename, 'FileType', 'text', ...
                            'Delimiter', '\t', 'ReadVariableNames', false);
    
    % Assign column names for clarity
    contactData.Properties.VariableNames = {'t', 'i', 'j', 'Ci', 'Cj'};
    
    % Extract unique node IDs
    nodeIDs = unique([contactData.i; contactData.j]);
    numNodes = length(nodeIDs);
    fprintf('Number of unique nodes: %d\n', numNodes);
    
    % Map node IDs to indices 1..N
    [~, sources] = ismember(contactData.i, nodeIDs);
    [~, targets] = ismember(contactData.j, nodeIDs);
    
    % You can define weights (e.g., 1 per contact)
    weights = ones(size(sources));
    
    %% --- Build sparse adjacency matrix ---
    fprintf('Building sparse adjacency matrix...\n');
    A_sparse = sparse(sources, targets, weights, numNodes, numNodes);
    
    % Symmetrize (since contacts are undirected)
    A_sparse = A_sparse + A_sparse';
    A_sparse(A_sparse > 0 & A_sparse < 1) = 1;  % make binary if needed
    
    %% --- Extract nodes with degree > 0 ---
    deg = full(sum(A_sparse, 2) + sum(A_sparse, 1)');
    nonzero_idx = find(deg > 0);
    nodeIDs_nonzero = nodeIDs(nonzero_idx);
    
    %% --- Build dense adjacency matrix for connected nodes ---
    fprintf('Creating dense adjacency matrix for %d connected nodes.\n', length(nonzero_idx));
    A_dense = full(A_sparse(nonzero_idx, nonzero_idx));
    
    fprintf('Adjacency matrices created successfully.\n');

end

function [] = plot_classes(filename,positions_out, nodeIDs_out)
    %% Load metadata
    meta = readtable(filename, 'Delimiter','\t','ReadVariableNames',false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    
    %% Assume you already have:
    % nodeIDs_out  -> N×1 string array of node IDs from positions extraction
    % positions_out -> N×3 numeric array of [x, y, z]
    
    % Find class for each node in positions_out
    numNodes = numel(nodeIDs_out);
    nodeClasses = strings(numNodes,1);
    
    for i = 1:numNodes
        idx = find(meta.ID == nodeIDs_out(i), 1);
        if ~isempty(idx)
            nodeClasses(i) = meta.Class(idx);
        else
            nodeClasses(i) = "Unknown";
        end
    end
    
    %% Assign distinct colors to each class
    %% Define the exact order of classes you want
    classOrder = unique(nodeClasses);
    numClasses = numel(classOrder);
    % uniqueClasses = unique(nodeClasses);
    % numClasses = numel(uniqueClasses);
    
    % Custom colormap: 10 visually distinct colors for 10 classes
    % You can tweak RGB values for better differentiation
    colors = [
        230 25 75;    % red
        60 180 75;    % green
        255 225 25;   % yellow
        0 130 200;    % blue
        245 130 48;   % orange
        % 145 30 180;   % purple
        % 70 240 240;   % cyan
        % 240 50 230;   % magenta
        % 210 245 60;   % lime
        % 250 190 190;  % pink
        % 0 0 0;
        ] / 255;      % normalize to [0,1]
    
    % Map node to color
    nodeColors = zeros(numNodes,3);
    for i = 1:numClasses
        nodeColors(nodeClasses==classOrder(i),:) = repmat(colors(i,:), sum(nodeClasses==classOrder(i)), 1);
    end
    
    %% Plot in 3D
    figure;
    
    % Plot actual nodes without adding them to the legend
    scatter(positions_out(:,1), positions_out(:,2), 100, nodeColors, 'filled', 'HandleVisibility','off');
    hold on;
    
    % Add dummy points for legend
    for i = 1:numClasses
        scatter3(NaN, NaN, NaN, 100, colors(i,:), 'filled');
    end
    
    % Add legend
    legend(classOrder, 'Location','bestoutside');
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Nodes Colored by Class');
    grid on; axis equal;
    hold off;

end

function [positions] = plot_centrality_colourvary(filename,A,centrality, X, nodeIDs_out,scaled)

    % --- Inputs for class coloring ---
    % nodeClasses: N×1 string array of class labels corresponding to nodes in A
    % classOrder: string array of unique classes in desired order
    % colors: N_classes×3 RGB array for each class

    G = graph(A);

    %% Load metadata
    meta = readtable(filename, 'Delimiter','\t','ReadVariableNames',false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    
    % Convert IDs to strings to match positions_out nodeIDs
    meta.ID = meta.ID;
    
    %% Assume you already have:
    % nodeIDs_out  -> N×1 string array of node IDs from positions extraction
    % positions_out -> N×3 numeric array of [x, y, z]
    
    % Find class for each node in positions_out
    numNodes = numel(nodeIDs_out);
    nodeClasses = strings(numNodes,1);
    
    for i = 1:numNodes
        idx = find(meta.ID == nodeIDs_out(i), 1);
        if ~isempty(idx)
            nodeClasses(i) = meta.Class(idx);
        else
            nodeClasses(i) = "Unknown";
        end
    end
    
    % Example (user provides):
    classOrder = unique(nodeClasses);
    
    % Custom colormap: 10 visually distinct colors for 10 classes
    % You can tweak RGB values for better differentiation
    colors = [
        230 25 75;    % red
        60 180 75;    % green
        255 225 25;   % yellow
        0 130 200;    % blue
        245 130 48;   % orange
        % 145 30 180;   % purple
        % 70 240 240;   % cyan
        % 240 50 230;   % magenta
        % 210 245 60;   % lime
        % 250 190 190;  % pink
        % 0 0 0;
        ] / 255;      % normalize to [0,1]
    
    % Map node to color based on class
    nodeColors = zeros(length(centrality),3);
    for i = 1:numel(classOrder)
        nodeColors(nodeClasses == classOrder(i), :) = repmat(colors(i,:), sum(nodeClasses == classOrder(i)), 1);
    end
    
    %% Plotting local centrality scaled nodes but colored by class
    figure;
    
    if nargin < 2 || isempty(X)
        p = plot(G, 'Layout','force', ...
            'MarkerSize', (centrality ./ sum(centrality)) * scaled, ... % scaling by centrality
            'NodeColor', nodeColors, ...                              % color by class
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0],'HandleVisibility','off');

    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
            'MarkerSize', (centrality ./ sum(centrality)) * scaled, ...
            'NodeColor', nodeColors, ...
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0],'HandleVisibility','off');
    end
    
    
    % --- Legend ---
    hold on;
    for i = 1:numel(classOrder)
        scatter3(NaN, NaN, NaN, 100, colors(i,:), 'filled'); % dummy for legend
    end
    legend(classOrder, 'Location','bestoutside');
    
    % title('Local Eigenvector Centrality Scaled Nodes Colored by Class');
    axis off; grid on; box on;

    positions = [p.XData',p.YData'];

end

[A_sparse1, A_dense1, nodeIDs_nonzero1] = csv_to_adjacency('thiers_2011.csv');
[centrality] = local_eigenvector_centrality(A_dense1,[],1,5);
[positions] = plot_centrality_colourvary('metadata_HS_2011.txt',A_dense1, centrality,[],nodeIDs_nonzero1,120);
[centrality] = local_eigenvector_centrality(A_dense1,[],0,1);
plot_centrality_colourvary('metadata_HS_2011.txt',A_dense1, centrality,[],nodeIDs_nonzero1,90);
plot_classes('metadata_HS_2011.txt',positions, nodeIDs_nonzero1)

% [A_sparse1, A_dense1, nodeIDs_nonzero1] = csv_to_adjacency('thiers_2012.csv');
% [centrality] = local_eigenvector_centrality(A_dense1,[],1,5);
% [positions] = plot_centrality_colourvary('metadata_HS_2012.txt',A_dense1, centrality,[],nodeIDs_nonzero1,150);
% [centrality] = local_eigenvector_centrality(A_dense1,[],0,1);
% plot_centrality_colourvary('metadata_HS_2012.txt',A_dense1, centrality,[],nodeIDs_nonzero1,90);
% plot_classes('metadata_HS_2012.txt',positions, nodeIDs_nonzero1)