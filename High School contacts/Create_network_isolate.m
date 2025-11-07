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

function [A_class, nodeIDs_class] = csv_to_class_adjacency(filename)
% csv_to_class_adjacency
% Reads a Thiers 2011-style contact CSV (columns: t, i, j, Ci, Cj)
% and builds a separate adjacency matrix for each class ID.
%
% Each adjacency matrix includes only nodes that belong to the same class,
% and only edges connecting nodes of that class.
%
% Outputs:
%   A_class       - struct, fields named by class ID (e.g. A_class.class1A)
%                   Each field contains a sparse adjacency matrix.
%   nodeIDs_class - struct, mapping each class to its corresponding node IDs.

    %% --- Read contact list ---
    fprintf('Reading contact CSV: %s\n', filename);
    contactData = readtable(filename, 'FileType', 'text', ...
                            'Delimiter', '\t', 'ReadVariableNames', false);
    contactData.Properties.VariableNames = {'t', 'i', 'j', 'Ci', 'Cj'};

    % Extract unique node IDs
    nodeIDs = unique([contactData.i; contactData.j]);
    numNodes = numel(nodeIDs);
    fprintf('Number of unique nodes: %d\n', numNodes);

    % Map node IDs to indices
    [~, sources] = ismember(contactData.i, nodeIDs);
    [~, targets] = ismember(contactData.j, nodeIDs);
    weights = ones(size(sources));

    %% --- Build global adjacency matrix ---
    fprintf('Building sparse adjacency matrix...\n');
    A_full = sparse(sources, targets, weights, numNodes, numNodes);
    A_full = A_full + A_full';  % symmetrize (undirected)
    A_full(A_full > 0 & A_full < 1) = 1;  % make binary if needed

    %% --- Determine class of each node ---
    % Since each contact involves Ci and Cj, we can assign class per node.
    classMap = containers.Map('KeyType','double','ValueType','char');

    for r = 1:height(contactData)
        i = contactData.i(r);
        j = contactData.j(r);
        Ci = string(contactData.Ci(r));
        Cj = string(contactData.Cj(r));

        if ~isKey(classMap, i)
            classMap(i) = Ci;
        end
        if ~isKey(classMap, j)
            classMap(j) = Cj;
        end
    end

    % Convert to arrays aligned with nodeIDs
    nodeClasses = strings(numNodes, 1);
    for n = 1:numNodes
        id = nodeIDs(n);
        if isKey(classMap, id)
            nodeClasses(n) = classMap(id);
        else
            nodeClasses(n) = "Unknown";
        end
    end

    %% --- Split by class ---
    uniqueClasses = unique(nodeClasses);
    A_class = struct();
    nodeIDs_class = struct();

    for c = 1:numel(uniqueClasses)
        className = uniqueClasses(c);
        if className == "Unknown"
            continue;
        end

        classMask = nodeClasses == className;
        classIdx = find(classMask);

        if isempty(classIdx)
            continue;
        end

        % Extract subgraph for this class
        A_sub = A_full(classIdx, classIdx);

        % Remove isolated nodes
        deg = full(sum(A_sub, 2) + sum(A_sub, 1)');
        nonzeroIdx = find(deg > 0);

        if isempty(nonzeroIdx)
            fprintf('Skipping class %s (no edges).\n', className);
            continue;
        end

        A_sub = A_sub(nonzeroIdx, nonzeroIdx);
        nodeIDs_sub = nodeIDs(classIdx(nonzeroIdx));

        % Store in struct (use valid field name)
        classStructname = matlab.lang.makeValidName("class" + className);
        A_class.(classStructname) = A_sub;
        nodeIDs_class.(classStructname) = nodeIDs_sub;

        fprintf('Class %s: %d nodes, %d edges retained.\n', ...
            classStructname, numel(nodeIDs_sub), nnz(A_sub)/2);
    end

    fprintf('Created adjacency matrices for %d classes.\n', numel(fieldnames(A_class)));
end

function [nodeIDs_out, positions_out] = csv_create_positions(filename, includeIDs)
% csv_create_positions
% Extracts node class information from a Thiers 2011-style CSV file and assigns
% 2D clustered positions grouped by class (e.g., PC, PC*, PSI*, teacher).
%
% Inputs:
%   filename   - string, path to CSV file with columns [t i j Ci Cj]
%   includeIDs - string array or cell array of node IDs to include
%
% Outputs:
%   nodeIDs_out   - M×1 string array of node IDs included
%   positions_out - table with columns [NodeID, X, Y, Z, Class]

    %% --- Setup ---
    if iscell(includeIDs)
        includeIDs = string(includeIDs);
    end

    % Read contact CSV
    contactData = readtable(filename, 'FileType', 'text', ...
                            'Delimiter', '\t', 'ReadVariableNames', false);
    contactData.Properties.VariableNames = {'t', 'i', 'j', 'Ci', 'Cj'};

    %% --- Determine class of each node ---
    nodeList = unique([contactData.i; contactData.j]);
    classMap = containers.Map('KeyType', 'double', 'ValueType', 'char');

    for r = 1:height(contactData)
        i = contactData.i(r);
        j = contactData.j(r);
        Ci = string(contactData.Ci(r));
        Cj = string(contactData.Cj(r));

        if ~isKey(classMap, i)
            classMap(i) = Ci;
        end
        if ~isKey(classMap, j)
            classMap(j) = Cj;
        end
    end

    %% --- Define class layout ---
    uniqueClasses = unique(string(values(classMap)))';
    numClasses = numel(uniqueClasses);

    % Layout grid spacing
    xSpacing = 7;
    ySpacing = 4;

    % Assign cluster centers (arranged in rows)
    centers = containers.Map;
    for c = 1:numClasses
        if c ~= 4
            xPos = mod(c-1, 3) * xSpacing - xSpacing;  % wrap every 3 classes
            yPos = mod(c-1, 2) * ySpacing;
        else
            xPos = 0;
            yPos = -ySpacing;
        end
        centers(uniqueClasses(c)) = [xPos, yPos];
    end

    rng(1); % reproducibility

    %% --- Generate positions ---
    nodeIDs_out = {};
    positions_out = [];
    classes_out = {};

    for n = 1:numel(nodeList)
        nodeID = nodeList(n);

        % Skip if not in includeIDs (if filtering)
        if ~ismember(nodeID, includeIDs)
            continue;
        end

        % Determine class
        if isKey(classMap, double(nodeList(n)))
            nodeClass = string(classMap(double(nodeList(n))));
        else
            nodeClass = "Unknown";
        end

        % --- Determine base position ---
        if isKey(centers, nodeClass)
            base = centers(nodeClass);
        else
            base = [0, 0];
        end

        % --- Add random jitter ---
        jitter = (rand(1,2) - 0.5) * 5;
        pos2D = base + [0.9 * jitter(1), 0.9 * jitter(2)];
        z = 0;

        % --- Store outputs ---
        nodeIDs_out{end+1,1} = nodeID;
        positions_out(end+1,:) = [pos2D, z];
        classes_out{end+1,1} = nodeClass;
    end

    %% --- Convert to table for output ---
    positions_out = table(nodeIDs_out, positions_out(:,1), positions_out(:,2), positions_out(:,3), classes_out, ...
        'VariableNames', {'NodeID', 'X', 'Y', 'Z', 'Class'});

    fprintf('Generated clustered positions for %d nodes across %d classes.\n', ...
        numel(nodeIDs_out), numClasses);
end



function [] = plot_classes(positions_out, nodeIDs_out)
    %% Load metadata
    meta = readtable('metadata_HS_2011.txt', 'Delimiter','\t','ReadVariableNames',false);
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
    uniqueClasses = ["PC","PC*","PSI*","teacher"];
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
    scatter3(positions_out(:,1), positions_out(:,2), positions_out(:,3), 200, nodeColors, 'filled', 'HandleVisibility','off');
    hold on;
    
    % Add dummy points for legend
    for i = 1:numClasses
        scatter3(NaN, NaN, NaN, 200, colors(i,:), 'filled');
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
    meta = readtable('metadata_HS_2011.txt', 'Delimiter','\t','ReadVariableNames',false);
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
    
    % Example (user provides):
    classOrder = ["PC","PC*","PSI*","teacher"];
    
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
            'EdgeAlpha', 0.03, 'EdgeColor', [0 0 0],'HandleVisibility','off', ...
            'NodeLabel', {});

    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
            'MarkerSize', (centrality ./ sum(centrality)) * scaled, ...
            'NodeColor', [.6,.6,.6], ... %nodeColors, ...  %                          % color by class
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

[A_sparse1, A_dense1, nodeIDs_nonzero1] = csv_to_adjacency('thiers_2011.csv');
[nodeIDs1, positions1] = csv_create_positions('thiers_2011.csv',nodeIDs_nonzero1);
allPositions = positions1{:, {'X', 'Y', 'Z'}};

[A_class, nodeIDs_class] = csv_to_class_adjacency('thiers_2011.csv');

classList = ["PC","PC*","PSI*","teacher"];

% for i = 1:length(classList)-1
%     % Access specific class:
%     className = matlab.lang.makeValidName(strcat("class",classList(i)));
%     A_1class = A_class.(className);
%     ids_1class = nodeIDs_class.(className);
%     % Example: Filter positions1 for entries from class '1A'
%     stringClass = string(positions1.Class);
%     idx = strcmp(stringClass, classList(i));
%     filtered_positions1 = positions1{idx, {'X', 'Y', 'Z'}};
%     [centrality] = local_eigenvector_centrality(A_1class,filtered_positions1,0,1);
%     plot_centrality_colourvary(A_1class,centrality,filtered_positions1,ids_1class,100);
%     axis equal
% end

stringClass = string(positions1.Class);
TMask = stringClass == "teacher";
TIdx = find(TMask);
% allPositions(TIdx(1),1:2) = [-7,-2];
allPositions(TIdx(2),1:2) = [-1,-3];
% allPositions(TIdx(3),1:2) = [-3,-2];
% allPositions(TIdx(4),1:2) = [-1,-2];
% allPositions(TIdx(5),1:2) = [1,-2];
% allPositions(TIdx(6),1:2) = [3,-2];
allPositions(TIdx(7),1:2) = [.25,-2.5];
% allPositions(TIdx(8),1:2) = [7,-2];

% plot_classes(positions1, nodeIDs_nonzero1)
[eigcentrality] = local_eigenvector_centrality(A_dense1,allPositions,1,5);

% page rank
G = graph(A_dense1);
PRank_centrality = centrality(G,'pagerank','FollowProbability',0.85,'Importance',G.Edges.Weight);
plot_centrality_colourvary(A_dense1, PRank_centrality,allPositions,nodeIDs_nonzero1,600);
axis equal

% Optimise wsd_power by varying pow
% Range of exponents to test
pow_list = 0.05:0.05:1;  % you can adjust step/range
wsd_power_values = zeros(size(pow_list));
addpath('..\wasserstein-distance-master')

for i = 1:length(pow_list)
    pow = pow_list(i);

    % Compute "warped" vector using current power
    warped = eigcentrality.^pow;
    warped = warped / sum(warped);  % normalize to sum 1
    
    difference_vector = PRank_centrality - warped;
    wsd_power_values(i) = norm(difference_vector, 2); 
end

% Find optimal power
[wsd_min, idx_min] = min(wsd_power_values);
pow = pow_list(idx_min);

warped = eigcentrality.^pow./sum(eigcentrality.^pow);
plot_centrality_colourvary(A_dense1, warped,allPositions,nodeIDs_nonzero1,700);
axis equal

plot_classes(allPositions, nodeIDs_nonzero1)

%% Boxplot comparison 
% Local_norm = (eigcentrality - mean(eigcentrality)) / std(eigcentrality);
% PR_norm = (PRank_centrality - mean(PRank_centrality)) / std(PRank_centrality);
% warped_norm = (warped - mean(warped)) / std(warped);

Local_norm = (eigcentrality - median(eigcentrality)) / mad(eigcentrality, 1);
PR_norm = (PRank_centrality - median(PRank_centrality)) / mad(PRank_centrality,1);
warped_norm = (warped - median(warped)) / mad(warped,1);

% --- Remove missing indices from all vectors ---
validIdx = setdiff(1:length(Local_norm), TIdx);  % keep only valid indices

Local_norm_valid = Local_norm(validIdx);
PR_norm_valid = PR_norm(validIdx);
warped_norm_valid = warped_norm(validIdx);

% --- Compute difference vectors using valid indices only ---
diff_vec_PRank = Local_norm_valid - PR_norm_valid;
diff_vec_warped = warped_norm_valid - PR_norm_valid;

[~, sortIdx] = sort(Local_norm_valid);  % sort by Local_norm_valid
Local_sorted = Local_norm_valid(sortIdx);
PR_sorted = PR_norm_valid(sortIdx);
warped_sorted = warped_norm_valid(sortIdx);

x = 1:length(Local_sorted);

% Make sure all are row vectors
x = x(:)'; 
Local_sorted = Local_sorted(:)'; 
PR_sorted = PR_sorted(:)'; 
warped_sorted = warped_sorted(:)'; 

%% Euclidean Distance

wsd_local = norm(PR_sorted - Local_sorted, 2)
wsd_power = norm(PR_sorted - warped_sorted, 2)

wsd_local = ws_distance(PR_sorted, Local_sorted, 1)
wsd_power = ws_distance(PR_sorted, warped_sorted, 1)

%% Figure
figure
subplot(1,9,[1:6]); hold on;

% --- Area between Local and PR ---
X_fill = [x, fliplr(x)];
Y_fill = [Local_sorted, fliplr(PR_sorted)];
fill(X_fill, Y_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% --- Area between Warped and PR ---
X_fill = [x, fliplr(x)];
Y_fill = [warped_sorted, fliplr(PR_sorted)];
fill(X_fill, Y_fill, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% % --- Overlay curves for clarity ---
% plot(x, Local_sorted, '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5);
% plot(x, warped_sorted, '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
% plot(x, PR_sorted, ':', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);

legend('Local − PageRank area',sprintf('Local (p=%1.2f) − PageRank area', round(pow,2)),'Location','northwest');
ylabel('Normalised Centrality');
xlabel('Sorted Index');
grid off;
xlim([0, 120])

% Plot differences as boxplots
% Prepare data for grouped boxplot
all_sets = [diff_vec_PRank(:); diff_vec_warped(:)];
grps = [repmat({'Local − PageRank'}, numel(diff_vec_PRank), 1);
        repmat({sprintf('Local (p=%1.2f) − PageRank', round(pow,2))}, numel(diff_vec_warped), 1)];

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

% [A_sparse2, A_dense2, nodeIDs_nonzero2] = gexf_to_adjacency('sp_data_school_day_2_g.gexf');
% [nodeIDs2, positions2] = gexf_create_positions('sp_data_school_day_2_g.gexf',nodeIDs_nonzero2);
% plot_classes(positions2, nodeIDs_nonzero2)
% [centrality] = local_eigenvector_centrality(A_dense2,positions2,0,10);
% plot_centrality_colourvary(A_dense2, centrality,positions2,nodeIDs_nonzero2,800);
% [centrality] = local_eigenvector_centrality(A_dense2,positions2,0,5);
% plot_centrality_colourvary(A_dense2, centrality,positions2,nodeIDs_nonzero2,800);