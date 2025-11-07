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

function [A_class, nodeIDs_class] = gexf_to_class_adjacency(filename)
% gexf_to_class_adjacency
% Reads a GEXF graph file and builds a separate adjacency matrix
% for each class ID.
%
% Each adjacency matrix includes only nodes that belong to the same class,
% and only edges connecting those nodes.
%
% Outputs:
%   A_class       - struct, with fields named by class ID (e.g. A_class.("1A"))
%                   Each field contains a sparse adjacency matrix.
%   nodeIDs_class - struct, mapping each class to its corresponding node IDs.

    fprintf('Reading GEXF file: %s\n', filename);
    xDoc = xmlread(filename);

    % --- Extract node list and class info ---
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength();
    nodeIDs = strings(numNodes, 1);
    nodeClasses = strings(numNodes, 1);

    for i = 1:numNodes
        nodeElem = nodeList.item(i-1);
        nodeIDs(i) = string(char(nodeElem.getAttribute('id')));

        % Default class
        nodeClass = "Unknown";

        % Extract <attvalue for="0"> which stores class
        attvalues = nodeElem.getElementsByTagName('attvalue');
        for j = 0:attvalues.getLength-1
            attElem = attvalues.item(j);
            if strcmp(char(attElem.getAttribute('for')), '0')
                nodeClass = string(char(attElem.getAttribute('value')));
                break;
            end
        end
        nodeClasses(i) = nodeClass;
    end

    % Map node IDs to indices
    idMap = containers.Map(nodeIDs, 1:numNodes);

    % --- Extract edges ---
    edgeList = xDoc.getElementsByTagName('edge');
    numEdges = edgeList.getLength();
    sources = zeros(numEdges, 1);
    targets = zeros(numEdges, 1);
    weights = ones(numEdges, 1);

    for i = 1:numEdges
        edgeElem = edgeList.item(i-1);
        srcID = string(char(edgeElem.getAttribute('source')));
        tgtID = string(char(edgeElem.getAttribute('target')));

        if isKey(idMap, srcID) && isKey(idMap, tgtID)
            sources(i) = idMap(srcID);
            targets(i) = idMap(tgtID);
        else
            warning('Edge %d references unknown node(s): %s -> %s', i, srcID, tgtID);
            continue;
        end

        % Optional: extract weight if present
        attvalues = edgeElem.getElementsByTagName('attvalue');
        for a = 0:attvalues.getLength-1
            att = attvalues.item(a);
            if strcmp(char(att.getAttribute('for')), '2')
                weights(i) = str2double(char(att.getAttribute('value')));
                break;
            end
        end
    end

    % --- Build full sparse adjacency matrix ---
    A_full = sparse(sources, targets, weights, numNodes, numNodes);
    A_full = A_full + A_full';  % Symmetrize

    % --- Split by class ---
    uniqueClasses = unique(nodeClasses);
    A_class = struct();
    nodeIDs_class = struct();

    for c = 1:numel(uniqueClasses)
        className = uniqueClasses(c);
        classMask = nodeClasses == className;
        classIdx = find(classMask);

        if isempty(classIdx)
            continue;
        end

        % Extract subgraph for this class
        A_sub = A_full(classIdx, classIdx);

        % Remove isolated nodes (optional)
        deg = full(sum(A_sub, 2) + sum(A_sub, 1)');
        nonzeroIdx = find(deg > 0);

        if isempty(nonzeroIdx)
            fprintf('Skipping class %s (no edges).\n', className);
            continue;
        end

        A_sub = A_sub(nonzeroIdx, nonzeroIdx);
        nodeIDs_sub = nodeIDs(classIdx(nonzeroIdx));

        % Store in struct
        classStructname = strcat("class",className);
        A_class.(classStructname) = A_sub;
        nodeIDs_class.(classStructname) = nodeIDs_sub;

        fprintf('Class %s: %d nodes, %d edges retained.\n', ...
            classStructname, numel(nodeIDs_sub), nnz(A_sub)/2);
    end

    fprintf('Created adjacency matrices for %d classes.\n', numel(fieldnames(A_class)));
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
    ySpacing = 2; % distance between year clusters
    xSpacing = 7;  % distance between A/B classes
    centers = containers.Map;

    for y = 1:numYears
        centers(uniqueClasses((y-1)*2 + 1)) = [-(xSpacing/2), -((y-1)*ySpacing)];
        centers(uniqueClasses((y-1)*2 + 2)) = [(xSpacing/2), -((y-1)*ySpacing)];
    end
    centers("Teachers") = [0, -(numYears*ySpacing)];

    rng(1); % reproducibility

    % Initialize output variables
    nodeIDs_out = {};
    positions_out = [];
    classes_out = {}; % To store class information
    
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
        jitter = (rand(1,2)-0.5) * 5;  % random jitter in [-1,1]
        pos2D = base + [0.9 * jitter(1), 0.25 * jitter(2)];
        z = 0;
    
        % --- Store outputs
        nodeIDs_out{end+1,1} = nodeID;
        positions_out(end+1,:) = [pos2D, z];
        classes_out{end+1,1} = nodeClass; % Store the class
    end
    
    % Convert to table for easier filtering
    positions_out = table(nodeIDs_out, positions_out(:,1), positions_out(:,2), positions_out(:,3), classes_out, ...
                       'VariableNames', {'NodeID', 'X', 'Y', 'Z', 'Class'});

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
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0],'HandleVisibility','off', ...
            'NodeLabel', {});

    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
            'MarkerSize', (centrality ./ sum(centrality)) * scaled, ...
            'NodeColor', nodeColors, ...  #[.8,.8,.8], ...                         % color by class
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0],'HandleVisibility','off', ...
            'NodeLabel', {});
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
    % colors = [
    %     0.094, 0.588, 0.784;   % teal / cyan-blue
    %     0.902, 0.624, 0.000;   % amber / mustard
    %     0.376, 0.188, 0.541;   % violet
    %     0.000, 0.620, 0.451;   % green-turquoise
    %     0.835, 0.369, 0.000;   % burnt orange
    %     0.000, 0.000, 0.000;   % black
    % ];
    basecolours = [
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
    ] / 255;  % normalize to [0,1]
    
    % Blend odd-even pairs
    nPairs = floor(size(basecolours,1) / 2);
    colors = zeros(nPairs, 3);

    for i = 1:nPairs
        colr_idx = (i-1)*2 + (1:2);
        colors(i,:) = mean(basecolours(colr_idx,:), 1);
    end
    colors = [colors;0,0,0];

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
            'EdgeAlpha', 0.03, 'EdgeColor', [0 0 0],'HandleVisibility','off', ...
            'NodeLabel', {});

    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
            'MarkerSize', (centrality ./ sum(centrality)) * scaled, ...
            'NodeColor', nodeColors, ...
            'EdgeAlpha', 0.03, 'EdgeColor', [0 0 0],'HandleVisibility','off', ...
            'NodeLabel', {});
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

% Day 1
[A_sparse1, A_dense1, nodeIDs_nonzero1] = gexf_to_adjacency('sp_data_school_day_1_g.gexf');
[nodeIDs1, positions1] = gexf_create_positions('sp_data_school_day_1_g.gexf',nodeIDs_nonzero1);
allPositions = positions1{:, {'X', 'Y', 'Z'}};

[A_class, nodeIDs_class] = gexf_to_class_adjacency('sp_data_school_day_1_g.gexf');

classList = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B"];

stringClass = string(positions1.Class);
TMask = stringClass == "Teachers";
TIdx = find(TMask);
allPositions(TIdx(1),1:2) = [0,-8];
allPositions(TIdx(2),1:2) = [0,-3];
allPositions(TIdx(3),1:2) = [0,-5];
allPositions(TIdx(4),1:2) = [0,-6];
allPositions(TIdx(5),1:2) = [0,-2];
allPositions(TIdx(6),1:2) = [0,0];
allPositions(TIdx(7),1:2) = [0,-4];
allPositions(TIdx(8),1:2) = [0,-8];
allPositions(TIdx(9),1:2) = [0,0];
allPositions(TIdx(10),1:2) = [0,-2];

% Initialize storage
allCentralityData = []; % will store nodeID and centrality pairs

for i = 1:length(classList)
    % Access specific class:
    A_1class = A_class.(strcat("class", classList(i)));
    ids_1class = nodeIDs_class.(strcat("class", classList(i)));

    % Filter positions1 for entries from current class
    stringClass = string(positions1.Class);
    idx = strcmp(stringClass, classList(i));
    filtered_positions1 = positions1{idx, {'X', 'Y', 'Z'}};

    % Compute centrality
    class_centrality = local_eigenvector_centrality(A_1class, filtered_positions1, 0, 1);

    % Store node IDs and centrality values
    classCentralityTable = table(ids_1class(:), class_centrality(:), ...
                                 'VariableNames', {'NodeID', 'Centrality'});
    allCentralityData = [allCentralityData; classCentralityTable];

    % % Visualise
    % plot_centrality_colourvary(A_1class, class_centrality, filtered_positions1, ids_1class, 250);
    % axis equal
end

% plot_classes(positions1, nodeIDs_nonzero1)
% [centrality] = local_eigenvector_centrality(A_dense1,allPositions,0,1);
G = graph(A_dense1);
PRank_centrality = centrality(G,'pagerank','FollowProbability',0.15,'Importance',G.Edges.Weight);
plot_centrality_colourvary(A_dense1, PRank_centrality,allPositions,nodeIDs_nonzero1,1400);
axis equal
[eigcentrality] = local_eigenvector_centrality(A_dense1,allPositions,0,5);
plot_centrality_colourvary_years(A_dense1, eigcentrality,allPositions,nodeIDs_nonzero1,1200);
axis equal
[eigcentrality] = local_eigenvector_centrality(A_dense1,allPositions,0,10);
plot_centrality_colourvary(A_dense1, eigcentrality,allPositions,nodeIDs_nonzero1,1200);
axis equal

% --- Add Local-centrality rows for missing nodes (TIdx) ---
% Create a table with zero values for all TIdx nodes
TeachCentralityTable = table(string(nodeIDs1(TIdx(:))), eigcentrality(TIdx), ...
                            'VariableNames', {'NodeID', 'Centrality'});

% Append the zero-centrality entries
allCentralityData = [allCentralityData; TeachCentralityTable];

[~, sortIdx] = ismember(string(nodeIDs1), allCentralityData.NodeID);

% Reorder the table
orderedCentralityData = allCentralityData(sortIdx, :);

% Extract ordered centrality vector
orderedCentrality = orderedCentralityData.Centrality;

plot_centrality_colourvary(A_dense1, orderedCentrality,allPositions,nodeIDs_nonzero1,1200);
axis equal

%%%%%%% Boxplot comparison %%%%%%%
% Local_norm = (eigcentrality - mean(eigcentrality)) / std(eigcentrality);
% PR_norm = (PRank_centrality - mean(PRank_centrality)) / std(PRank_centrality);
% Class_norm = (orderedCentrality - mean(orderedCentrality)) / std(orderedCentrality);

Local_norm = (eigcentrality - median(eigcentrality)) / mad(eigcentrality, 1);
PR_norm = (PRank_centrality - median(PRank_centrality)) / mad(PRank_centrality,1);
Class_norm = (orderedCentrality - median(orderedCentrality)) / mad(orderedCentrality,1);

% --- Remove missing indices from all vectors ---
validIdx = setdiff(1:length(Local_norm), TIdx);  % keep only valid indices

Local_norm_valid = Local_norm(validIdx);
PR_norm_valid = PR_norm(validIdx);
Class_norm_valid = Class_norm(validIdx);

% --- Compute difference vectors using valid indices only ---
diff_vec_PRank = Local_norm_valid - PR_norm_valid;
diff_vec_class = Local_norm_valid - Class_norm_valid;

% --- Step 3: Plot ---
[~, sortIdx] = sort(Local_norm_valid);  % sort by Local_norm_valid
Local_sorted = Local_norm_valid(sortIdx);
PR_sorted = PR_norm_valid(sortIdx);
Class_sorted = Class_norm_valid(sortIdx);

x = 1:length(Local_sorted);

% Make sure all are row vectors
x = x(:)'; 
Local_sorted = Local_sorted(:)'; 
PR_sorted = PR_sorted(:)'; 
Class_sorted = Class_sorted(:)'; 

figure
subplot(1,9,[1:6]); hold on;

% --- Area between Local and Class ---
X_fill = [x, fliplr(x)];
Y_fill = [Local_sorted, fliplr(Class_sorted)];
fill(X_fill, Y_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% --- Area between Local and PR ---
X_fill = [x, fliplr(x)];
Y_fill = [Local_sorted, fliplr(PR_sorted)];
fill(X_fill, Y_fill, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% --- Overlay curves for clarity ---
% plot(x, Local_sorted, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
% plot(x, Class_sorted, '--', 'Color', [0 0.4470 0.7410]*0.6 + 0.4, 'LineWidth', 1.5); % slightly different shade
% plot(x, PR_sorted, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

legend('Local − Class area','Local − PageRank area','Location', 'northwest');
ylabel('Normalised Centrality');
xlabel('Sorted Index');
grid off;
xlim([0, 227])


% Plot differences as boxplots
% Prepare data for grouped boxplot
all_sets = [diff_vec_class(:); diff_vec_PRank(:)];
grps = [repmat({'Local − Class'}, numel(diff_vec_class), 1);
        repmat({'Local − PageRank'}, numel(diff_vec_PRank), 1)];

%%%%%
subplot(1,9,[7:9]);
hold on;

% Unique groups and numeric positions
uniqueGroups = unique(grps, 'stable');
numGroups = numel(uniqueGroups);
positions = 1:numGroups; % numeric X positions for boxchart

% Assign colors for each group
colors = [0, 0.4470, 0.7410;  % blue
          0.8500, 0.3250, 0.0980]; % orange
alpha = 0.3;  % how much to blend with white (0 = original, 1 = white)
lighterColors = colors + (1 - colors) * alpha;

% Plot boxchart for each group individually
for i = 1:numGroups
    % Extract data for this group
    idx = strcmp(grps, uniqueGroups{i});
    bc = boxchart(positions(i)*ones(sum(idx),1), all_sets(idx), ...
                  'BoxWidth', .7, 'Notch', 'off');
    bc.BoxFaceColor = lighterColors(i,:);
    bc.LineWidth = 1.5;
    bc.MarkerColor = lighterColors(i,:);
end

% Scatter points overlay
[~,~,ic] = unique(grps, 'stable');
for i = 1:numGroups
    Ind = find(ic == i);
    scatter(positions(i) + randn(size(Ind))*0.1, all_sets(Ind), ...
        25, colors(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.3);
end

% Format axes
ylabel('Normalised Centrality Difference');
set(gca, 'XTick', positions, 'XTickLabel', uniqueGroups);
grid off;
%%%%%%%

addpath('..\wasserstein-distance-master')

% %% Create CDF for Wasserstein Distance
% x = 1:length(eigcentrality);
% y1 = eigcentrality;
% y2 = PRank_centrality;
% y3 = orderedCentrality;
% 
% % Compute PDFs
% y1 = y1 ./ trapz(x, y1);
% y2 = y2 ./ trapz(x, y2);
% y3 = y3 ./ trapz(x, y3);
% 
% % Compute CDFs
% F1 = cumtrapz(x, y1);
% F2 = cumtrapz(x, y2);
% F3 = cumtrapz(x, y3);
% 
% wsd_PR = ws_distance(F1, F2, 1)
% wsd_class = ws_distance(F1, F3, 1)

wsd_PR = norm(Local_sorted - PR_sorted, 2)
wsd_class = norm(Local_sorted - Class_sorted, 2)

wsd_PR = ws_distance(Local_sorted, PR_sorted, 1)
wsd_class = ws_distance(Local_sorted, Class_norm, 1)

plot_classes(allPositions, nodeIDs_nonzero1)

% % Day 2
% [A_sparse1, A_dense1, nodeIDs_nonzero1] = gexf_to_adjacency('sp_data_school_day_2_g.gexf');
% [nodeIDs1, positions1] = gexf_create_positions('sp_data_school_day_2_g.gexf',nodeIDs_nonzero1);
% allPositions = positions1{:, {'X', 'Y', 'Z'}};
% 
% [A_class, nodeIDs_class] = gexf_to_class_adjacency('sp_data_school_day_2_g.gexf');
% 
% classList = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B"];
% 
% stringClass = string(positions1.Class);
% TMask = stringClass == "Teachers";
% TIdx = find(TMask);
% allPositions(TIdx(1),1:2) = [0,-8];
% allPositions(TIdx(2),1:2) = [0,-3];
% allPositions(TIdx(3),1:2) = [0,-5];
% allPositions(TIdx(4),1:2) = [0,-6];
% allPositions(TIdx(5),1:2) = [0,-2];
% allPositions(TIdx(6),1:2) = [0,0];
% allPositions(TIdx(7),1:2) = [0,-4];
% allPositions(TIdx(8),1:2) = [0,-8];
% allPositions(TIdx(9),1:2) = [0,0];
% allPositions(TIdx(10),1:2) = [0,-2];
% 
% 
% % plot_classes(positions1, nodeIDs_nonzero1)
% [centrality] = local_eigenvector_centrality(A_dense1,allPositions,0);
% plot_centrality_colourvary(A_dense1, centrality,allPositions,nodeIDs_nonzero1,1200);
% axis equal
% [centrality] = local_eigenvector_centrality(A_dense1,allPositions,0,5);
% plot_centrality_colourvary_years(A_dense1, centrality,allPositions,nodeIDs_nonzero1,1200);
% axis equal
% [centrality] = local_eigenvector_centrality(A_dense1,allPositions,0,10);
% plot_centrality_colourvary(A_dense1, centrality,allPositions,nodeIDs_nonzero1,1200);
% axis equal
% % 
% % plot_classes(allPositions, nodeIDs_nonzero1)