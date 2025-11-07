addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Local_eigenvector_centrality'

function [A_sparse, A_dense, nodeIDs_nonzero] = gexf_to_adjacency(filename)
% gexf_to_adjacency
% Reads a GEXF (XML) graph file and builds:
%   - A sparse adjacency matrix (using all nodes and IDs)
%   - A dense adjacency matrix for nodes with degree > 0
% Edge weight is taken from <attvalue for="2" value="..."/>.

    % --- Parse XML ---
    fprintf('Reading GEXF file: %s\n', filename);
    xDoc = xmlread(filename);

    % --- Extract node IDs ---
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength();
    nodeIDs = strings(numNodes, 1);
    for i = 1:numNodes
        nodeElement = nodeList.item(i-1);
        nodeIDs(i) = string(char(nodeElement.getAttribute('id')));
    end

    % Map node ID to index
    idMap = containers.Map(nodeIDs, 1:numNodes);

    % --- Extract edges ---
    edgeList = xDoc.getElementsByTagName('edge');
    numEdges = edgeList.getLength();
    sources = zeros(numEdges, 1);
    targets = zeros(numEdges, 1);
    weights = ones(numEdges, 1); % default = 1

    for i = 1:numEdges
        edgeElement = edgeList.item(i-1);
        srcID = string(char(edgeElement.getAttribute('source')));
        tgtID = string(char(edgeElement.getAttribute('target')));

        if isKey(idMap, srcID) && isKey(idMap, tgtID)
            sources(i) = idMap(srcID);
            targets(i) = idMap(tgtID);
        else
            warning('Edge %d references unknown node(s): %s -> %s', i, srcID, tgtID);
            continue
        end

        % --- Extract weight from <attvalue for="2"> ---
        attvalues = edgeElement.getElementsByTagName('attvalue');
        for a = 0:attvalues.getLength-1
            att = attvalues.item(a);
            if strcmp(char(att.getAttribute('for')), '2')
                weights(i) = str2double(char(att.getAttribute('value')));
                break
            end
        end
    end

    % --- Build sparse adjacency matrix ---
    fprintf('Building sparse adjacency matrix...\n');
    A_sparse = sparse(sources, targets, weights, numNodes, numNodes);

    % Symmetrize for undirected graph
    A_sparse = A_sparse + A_sparse';
    A_sparse(A_sparse > 0 & A_sparse < 1) = 1; % optional if only binary

    % --- Extract nodes with degree > 0 ---
    deg = full(sum(A_sparse, 2) + sum(A_sparse, 1)');
    nonzero_idx = find(deg > 0);
    nodeIDs_nonzero = nodeIDs(nonzero_idx);

    % --- Build dense adjacency matrix for connected nodes ---
    fprintf('Creating dense adjacency matrix for %d connected nodes.\n', length(nonzero_idx));
    A_dense = full(A_sparse(nonzero_idx, nonzero_idx));

    fprintf('Adjacency matrices created successfully.\n');
end

function [nodeIDs_out, positions_out] = gexf_extract_positions(filename, includeIDs)
% gexf_extract_positions
% Extracts positions of nodes from a GEXF XML file, filtered by a list of node IDs.
%
% Usage:
%   [nodeIDs_out, positions_out] = gexf_extract_positions(filename, includeIDs);
%
% Inputs:
%   filename   - string, path to GEXF file
%   includeIDs - string array or cell array of node IDs to include
%
% Outputs:
%   nodeIDs_out   - M×1 string array of node IDs included
%   positions_out - M×3 numeric array of [x, y, z] positions

    % Convert includeIDs to string array if cell
    if iscell(includeIDs)
        includeIDs = string(includeIDs);
    end

    % Read XML
    xDoc = xmlread(filename);

    % Get all <node> elements
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength;

    % Initialize outputs
    nodeIDs_out = [];
    positions_out = [];

    for i = 0:numNodes-1
        nodeElem = nodeList.item(i);
        nodeID = string(char(nodeElem.getAttribute('id')));

        % Skip node if not in includeIDs
        if ~ismember(nodeID, includeIDs)
            continue;
        end

        % Extract <ns0:position>
        positionElems = nodeElem.getElementsByTagName('ns0:position');
        if positionElems.getLength > 0
            posElem = positionElems.item(0);
            x = str2double(char(posElem.getAttribute('x')));
            y = str2double(char(posElem.getAttribute('y')));
            z = str2double(char(posElem.getAttribute('z')));
        else
            x = NaN; y = NaN; z = NaN;
        end

        % Append to outputs
        nodeIDs_out(end+1,1) = nodeID;
        positions_out(end+1,:) = [x, y, z];
    end

    fprintf('Extracted positions for %d nodes from the input list.\n', numel(nodeIDs_out));
end

function [nodeIDs_out, positions_out] = gexf_create_positions(filename, includeIDs)
% gexf_extract_positions
% Extracts node class information from a GEXF file and assigns new
% 2D positions clustered by year and class.
%
% Inputs:
%   filename   - string, path to GEXF file
%   includeIDs - string array or cell array of node IDs to include
%
% Outputs:
%   nodeIDs_out   - M×1 string array of node IDs included
%   positions_out - M×3 numeric array of [x, y, z] positions

    % Convert includeIDs to string array if cell
    if iscell(includeIDs)
        includeIDs = string(includeIDs);
    end

    % Read XML
    xDoc = xmlread(filename);
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength;

    % Initialize outputs
    nodeIDs_out = [];
    positions_out = [];

    % Define all possible classes
    uniqueClasses = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];

    % Define cluster centers for each class
    numYears = 5;
    ySpacing = 10; % distance between year clusters
    xSpacing = 2;  % distance between A/B classes
    centers = containers.Map;

    for y = 1:numYears
        centers(uniqueClasses((y-1)*2 + 1)) = [-(xSpacing/2), -((y-1)*ySpacing)];
        centers(uniqueClasses((y-1)*2 + 2)) = [(xSpacing/2), -((y-1)*ySpacing)];
    end
    centers("Teachers") = [0, -(numYears*ySpacing + 8)];

    rng(1); % reproducibility

    for i = 0:numNodes-1
        nodeElem = nodeList.item(i);
        nodeID = string(char(nodeElem.getAttribute('id')));

        % Skip if not in includeIDs
        if ~ismember(nodeID, includeIDs)
            continue;
        end

        % --- Extract class info from attvalue with for="0"
        attvalues = nodeElem.getElementsByTagName('attvalue');
        nodeClass = "Teachers"; % default fallback

        for j = 0:attvalues.getLength-1
            attElem = attvalues.item(j);
            attrFor = string(char(attElem.getAttribute('for')));
            if attrFor == "0"
                nodeClass = string(char(attElem.getAttribute('value')));
                break;
            end
        end

        % --- Determine base position
        if isKey(centers, nodeClass)
            base = centers(nodeClass);
        else
            base = [0, 0]; % fallback if unknown
        end

        % --- Random jitter around base position
        jitter = (rand(1,2)-0.5) * 2;  % random jitter in [-1,1]
        pos2D = base; % + 0.8 * jitter;   % cluster tightness
        pos2D(1) = pos2D(1) + 0.8 * jitter(1);
        pos2D(2) = pos2D(2) + 2 * jitter(2);
        z = 0;

        % --- Store outputs
        nodeIDs_out(end+1,1) = nodeID;
        positions_out(end+1,:) = [pos2D, z];
    end

    fprintf('Generated clustered positions for %d nodes.\n', numel(nodeIDs_out));
end


function [] = plot_classes(positions_out, nodeIDs_out)
    %% Load metadata
    meta = readtable('metadata.txt', 'Delimiter','\t','ReadVariableNames',false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    
    % Convert IDs to strings to match positions_out nodeIDs
    meta.ID = string(meta.ID);
    
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
    uniqueClasses = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];
    numClasses = numel(uniqueClasses);
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
        145 30 180;   % purple
        70 240 240;   % cyan
        240 50 230;   % magenta
        210 245 60;   % lime
        250 190 190;  % pink
        0 0 0;
        ] / 255;      % normalize to [0,1]
    
    % Map node to color
    nodeColors = zeros(numNodes,3);
    for i = 1:numClasses
        nodeColors(nodeClasses==uniqueClasses(i),:) = repmat(colors(i,:), sum(nodeClasses==uniqueClasses(i)), 1);
    end
    
    %% Plot in 3D
    figure;
    
    % Plot actual nodes without adding them to the legend
    scatter3(positions_out(:,1), positions_out(:,2), positions_out(:,3), 100, nodeColors, 'filled', 'HandleVisibility','off');
    hold on;
    
    % Add dummy points for legend
    for i = 1:numClasses
        scatter3(NaN, NaN, NaN, 100, colors(i,:), 'filled');
    end
    
    % Add legend
    legend(uniqueClasses, 'Location','bestoutside');
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Nodes Colored by Class');
    grid on; axis equal;
    hold off;

end

function [] = plot_centrality_colourvary(A,centrality, X, nodeIDs_out,scaled)

    % --- Inputs for class coloring ---
    % nodeClasses: N×1 string array of class labels corresponding to nodes in A
    % classOrder: string array of unique classes in desired order
    % colors: N_classes×3 RGB array for each class

    G = graph(A);

    %% Load metadata
    meta = readtable('metadata.txt', 'Delimiter','\t','ReadVariableNames',false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    
    % Convert IDs to strings to match positions_out nodeIDs
    meta.ID = string(meta.ID);
    
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
    classOrder = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];
    
    % Custom colormap: 10 visually distinct colors for 10 classes
    % You can tweak RGB values for better differentiation
    colors = [
        230 25 75;    % red
        60 180 75;    % green
        255 225 25;   % yellow
        0 130 200;    % blue
        245 130 48;   % orange
        145 30 180;   % purple
        70 240 240;   % cyan
        240 50 230;   % magenta
        210 245 60;   % lime
        250 190 190;  % pink
        0 0 0;
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

end

function [] = plot_centrality_colourvary_years(A,centrality, X, nodeIDs_out,scaled)

    % --- Inputs for class coloring ---
    % nodeClasses: N×1 string array of class labels corresponding to nodes in A
    % classOrder: string array of unique classes in desired order
    % colors: N_classes×3 RGB array for each class

    G = graph(A);

    %% Load metadata
    meta = readtable('metadata.txt', 'Delimiter','\t','ReadVariableNames',false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    
    % Convert IDs to strings to match positions_out nodeIDs
    meta.ID = string(meta.ID);
    
    %% Assume you already have:
    % nodeIDs_out  -> N×1 string array of node IDs from positions extraction
    % positions_out -> N×3 numeric array of [x, y, z]
    
    % Find class for each node in positions_out
    numNodes = numel(nodeIDs_out);
    nodeClasses = strings(numNodes,1);
    
    for i = 1:numNodes
        idx = find(meta.ID == nodeIDs_out(i), 1);
        if ~isempty(idx)
            % if meta.Class(idx) contain
            
            nodeClasses(i) = strcat("Year ",meta.Class{idx}(1));
        else
            nodeClasses(i) = "Unknown";
        end
    end
    
    % Example (user provides):
    classOrder = ["Year 1","Year 2","Year 3","Year 4","Year 5","Teachers"];
    
    % Custom colormap: 10 visually distinct colors for 10 classes
    % You can tweak RGB values for better differentiation
    colors = [
        0.094, 0.588, 0.784;   % teal / cyan-blue
        0.902, 0.624, 0.000;   % amber / mustard
        0.376, 0.188, 0.541;   % violet
        0.000, 0.620, 0.451;   % green-turquoise
        0.835, 0.369, 0.000;   % burnt orange
        0.000, 0.000, 0.000;   % black
    ];

    
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

end

[A_sparse1, A_dense1, nodeIDs_nonzero1] = gexf_to_adjacency('sp_data_school_day_1_g.gexf');
[nodeIDs1, positions1] = gexf_create_positions('sp_data_school_day_1_g.gexf',nodeIDs_nonzero1);
plot_classes(positions1, nodeIDs_nonzero1)
[centrality] = local_eigenvector_centrality(A_dense1,positions1,0,1);
plot_centrality_colourvary(A_dense1, centrality,positions1,nodeIDs_nonzero1,800);
[centrality] = local_eigenvector_centrality(A_dense1,positions1,0,5);
plot_centrality_colourvary_years(A_dense1, centrality,positions1,nodeIDs_nonzero1,800);
[centrality] = local_eigenvector_centrality(A_dense1,positions1,0,10);
plot_centrality_colourvary(A_dense1, centrality,positions1,nodeIDs_nonzero1,800);

% [A_sparse2, A_dense2, nodeIDs_nonzero2] = gexf_to_adjacency('sp_data_school_day_2_g.gexf');
% [nodeIDs2, positions2] = gexf_create_positions('sp_data_school_day_2_g.gexf',nodeIDs_nonzero2);
% plot_classes(positions2, nodeIDs_nonzero2)
% [centrality] = local_eigenvector_centrality(A_dense2,positions2,0,10);
% plot_centrality_colourvary(A_dense2, centrality,positions2,nodeIDs_nonzero2,800);
% [centrality] = local_eigenvector_centrality(A_dense2,positions2,0,5);
% plot_centrality_colourvary(A_dense2, centrality,positions2,nodeIDs_nonzero2,800);