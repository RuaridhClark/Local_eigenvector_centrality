%% main_analysis.m
% Reproducible analysis script for "A local eigenvector centrality"
%
% Usage: run this script from its folder. It will load the GEXF file,
% compute centralities, and produce plots corresponding with the paper.

%% Configuration
clear; clc; close all;

% Include folder with local_eigenvector_centrality.m
addpath("../")

% Filenames
gexfFile = 'sp_data_school_day_2_g.gexf';
metaFile = 'metadata.txt';

% Check files
assert(isfile(gexfFile), 'GEXF file not found: %s', gexfFile);
assert(isfile(metaFile),  'Metadata file not found: %s', metaFile);

if exist('local_eigenvector_centrality', 'file') ~= 2
    warning(['local_eigenvector_centrality.m was not found on the MATLAB path.\n' ...
             'Please add it with addpath(...) before running the script.']);
end

%% Load adjacency (sparse) and dense subgraph for non-isolated nodes
fprintf('\n--- Reading GEXF and building adjacency matrices ---\n');
[A_sparse, A_dense, nodeIDs_nonzero] = gexf_to_adjacency(gexfFile);

% Create a graph for layout and PageRank
G_dense = graph(A_dense);

%% Generate (clustered) node positions for plotting
fprintf('\n--- Generating node positions (clustered by class/year) ---\n');
[nodeIDs_pos, positions_table] = gexf_create_positions(gexfFile, nodeIDs_nonzero);
Pos = positions_table{:, {'X', 'Y', 'Z'}};    % Nx3 numeric

%% Create class-separated adjacency matrices
fprintf('\n--- Building class adjacency matrices ---\n');
[A_classes, nodeIDs_class] = gexf_to_class_adjacency(gexfFile);

%% Compute PageRank centrality (for comparison)
fprintf('\n--- Computing PageRank centrality ---\n');
if any(G_dense.Edges.Weight)
    PRank_centrality = centrality(G_dense, 'pagerank', 'FollowProbability', 0.85, 'Importance', G_dense.Edges.Weight);
else
    PRank_centrality = centrality(G_dense, 'pagerank', 'FollowProbability', 0.85);
end

%% Compute local eigenvector centrality for whole dense network (example Imax choices)
fprintf('\n--- Computing local eigenvector centralities (Imax = 5 and 10) ---\n');
[local_centrality_5,  ~] = local_eigenvector_centrality(A_dense, Pos, false, 5);
[local_centrality_10, ~] = local_eigenvector_centrality(A_dense, Pos, false, 10);

%% Plot PageRank and Local centrality maps coloured by class
fprintf('\n--- Plotting centrality maps ---\n');
plot_centrality_colourvary(A_dense, PRank_centrality, Pos, nodeIDs_nonzero, 15, metaFile);
title('PageRank'); axis equal;

plot_centrality_colourvary_years(A_dense, local_centrality_5, Pos, nodeIDs_nonzero, 20, metaFile);
title('Local (i=5)'); axis equal;

plot_centrality_colourvary(A_dense, local_centrality_10, Pos, nodeIDs_nonzero, 20, metaFile);
title('Local (i=10)'); axis equal;

%% Compute class-level centralities and assemble ordered vector
fprintf('\n--- Computing class-level centralities (per-class local centrality) ---\n');
classList = string(fieldnames(nodeIDs_class));
allCentralityData = table('Size', [0 2], 'VariableTypes', {'string','double'}, 'VariableNames', {'NodeID','Centrality'});

numClass = numel(classList);
for k = 1:numClass
    classStructname = classList(k);
    ids_1class = nodeIDs_class.(classStructname);
    A_class = A_classes.(classStructname);

    if k < numClass
        % Compute class-level local centrality (Imax = 1 -> principal eigenvector of class subgraph)
        [class_centrality, ~] = local_eigenvector_centrality(A_class, [], false, 1);
    else % Teachers nodes given local_centrality_10 values
        class_centrality = local_centrality_10(ismember(nodeIDs_nonzero,nodeIDs_class.class_Teachers));
    end

    % Store results (NodeID order corresponds to ids_1class)
    T = table(string(ids_1class(:)), class_centrality(:), 'VariableNames', {'NodeID','Centrality'});
    allCentralityData = [allCentralityData; T];
end

% Ensure Teacher nodes present: use local_centrality_10 values for missing teacher indices
% Map nodeIDs_nonzero to their centrality entries (fallback to local_centrality_10)
[commonIdx, ia, ib] = intersect(nodeIDs_nonzero, allCentralityData.NodeID, 'stable');
orderedCentrality = zeros(size(nodeIDs_nonzero));

% Fill known entries
if ~isempty(commonIdx)
    orderedCentrality(ia) = allCentralityData.Centrality(ib);
end

% Final plot of class-based centrality
plot_centrality_colourvary(A_dense, orderedCentrality, Pos, nodeIDs_nonzero, 15, metaFile);
title('Class (per-class centrality)'); axis equal;

%% Boxplot comparisons (Local vs PageRank vs Class)
fprintf('\n--- Boxplot and distance comparison ---\n');
% Normalise using median & MAD (robust)
Local_norm = (local_centrality_10 - median(local_centrality_10)) / mad(local_centrality_10, 1);
PR_norm    = (PRank_centrality - median(PRank_centrality)) / mad(PRank_centrality, 1);
Class_norm = (orderedCentrality - median(orderedCentrality)) / mad(orderedCentrality, 1);

% Remove teacher indices if necessary (find from metadata)
meta = readtable(metaFile, 'Delimiter', '\t', 'ReadVariableNames', false);
meta.Properties.VariableNames = {'ID','Class','Gender'};
meta.ID = string(meta.ID);
Tmask = startsWith(string(positions_table.Class), 'Teachers');
TIdx = find(Tmask);

validIdx = setdiff(1:length(Local_norm), TIdx);
Local_valid = Local_norm(validIdx)';
PR_valid    = PR_norm(validIdx)';
Class_valid = Class_norm(validIdx)';

% Sort by Local for plotting
[Local_sorted, sidx] = sort(Local_valid);
PR_sorted    = PR_valid(sidx);
Class_sorted = Class_valid(sidx);

% Figure layout
figure; %('Units','normalized','Position',[0.05 0.05 0.9 0.7])

t = tiledlayout(1,9, 'TileSpacing', 'compact', 'Padding', 'compact');

% Difference area plots on left
nexttile(1,[1 6]); hold on;
x = 1:numel(Local_sorted);

% Area between Local and Class
X_fill = [x, fliplr(x)];
Y_fill = [Local_sorted, fliplr(Class_sorted)];
fill(X_fill, Y_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Area between Local and PageRank
Y_fill = [Local_sorted, fliplr(PR_sorted)];
fill(X_fill, Y_fill, [0.8500 0.3250 0.0980], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

legend('Local − Class area','Local − PageRank area','Location','northwest');
ylabel('Normalised Centrality'); xlabel('Sorted Index');
box off; axis tight;

% Boxcharts on right
nexttile(7,[1 3]); hold on;
% Differences
diff_vec_class  = Local_valid - Class_valid;
diff_vec_PRank  = Local_valid - PR_valid;
all_sets = [diff_vec_class(:); diff_vec_PRank(:)];
grps = [repmat({'Local − Class'}, numel(diff_vec_class), 1);
        repmat({'Local − PageRank'}, numel(diff_vec_PRank), 1)];

uniqueGroups = unique(grps, 'stable');
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
for i = 1:numel(uniqueGroups)
    idx = strcmp(grps, uniqueGroups{i});
    bc = boxchart(i*ones(sum(idx),1), all_sets(idx), 'BoxWidth', .7);
    bc.BoxFaceColor = colors(i,:);
    bc.LineWidth = 1.2;
end
scatter(ones(size(diff_vec_class))*1 + randn(size(diff_vec_class))*0.05, diff_vec_class, 16, colors(1,:), 'filled', 'MarkerFaceAlpha', 0.2);
scatter(ones(size(diff_vec_PRank))*2 + randn(size(diff_vec_PRank))*0.05, diff_vec_PRank, 16, colors(2,:), 'filled', 'MarkerFaceAlpha', 0.2);
set(gca, 'XTick', [1 2], 'XTickLabel', uniqueGroups);
ylabel('Normalised Centrality Difference'); box off;

% Euclidean norms
EN_PR = norm(Local_sorted - PR_sorted, 2);
EN_class = norm(Local_sorted - Class_sorted, 2);
fprintf('Euclidean norm (Local - PageRank): %.4f\n', EN_PR);
fprintf('Euclidean norm (Local - Class): %.4f\n', EN_class);

fprintf('\n--- Analysis complete. ---\n');

%% ----------------------
% Local helper functions
%% ----------------------

function [A_sparse, A_dense, nodeIDs_nonzero] = gexf_to_adjacency(filename)
% Reads a GEXF file and returns sparse adjacency (all nodes) and dense
% adjacency restricted to nodes with degree>0. Edge weights are read from
% <attvalue for="2" value="..."/> when present; otherwise weight=1.
    fprintf('Reading GEXF file: %s\n', filename);
    xDoc = xmlread(filename);

    % Nodes
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength();
    nodeIDs = string(zeros(numNodes,1));
    for ii = 1:numNodes
        nodeIDs(ii) = string(char(nodeList.item(ii-1).getAttribute('id')));
    end
    idMap = containers.Map(nodeIDs, 1:numNodes);

    % Edges
    edgeList = xDoc.getElementsByTagName('edge');
    numEdges = edgeList.getLength();
    sources = zeros(numEdges,1); targets = zeros(numEdges,1);
    weights = ones(numEdges,1);

    for ii = 1:numEdges
        e = edgeList.item(ii-1);
        srcID = string(char(e.getAttribute('source')));
        tgtID = string(char(e.getAttribute('target')));
        if isKey(idMap, srcID) && isKey(idMap, tgtID)
            sources(ii) = idMap(srcID);
            targets(ii) = idMap(tgtID);
        else
            warning('Edge %d references unknown node(s): %s -> %s', ii, srcID, tgtID);
            continue;
        end
        % look for attvalue for="2"
        attvalues = e.getElementsByTagName('attvalue');
        for a = 0:attvalues.getLength-1
            att = attvalues.item(a);
            if strcmp(char(att.getAttribute('for')), '2')
                weights(ii) = str2double(char(att.getAttribute('value')));
                break;
            end
        end
    end

    % Build sparse adjacency and symmetrise
    A_sparse = sparse(sources, targets, weights, numNodes, numNodes);
    A_sparse = A_sparse + A_sparse';
    A_sparse(A_sparse > 0 & A_sparse < 1) = 1; % optional binary threshold

    % Nodes with degree > 0
    deg = full(sum(A_sparse,2) + sum(A_sparse,1)');
    nonzero_idx = find(deg > 0);
    nodeIDs_nonzero = nodeIDs(nonzero_idx);

    A_dense = full(A_sparse(nonzero_idx, nonzero_idx));
    fprintf('Adjacency matrices created: %d nodes, %d edges (dense).\n', size(A_dense,1), nnz(A_dense)/2);
end

function [A_class, nodeIDs_class] = gexf_to_class_adjacency(filename)
% Splits input graph by node attribute "for=0" (class). Returns struct of
% sparse adjacency matrices for each class and corresponding node ID lists.
    fprintf('Reading GEXF file for class splitting: %s\n', filename);
    xDoc = xmlread(filename);

    % Nodes and classes
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength();
    nodeIDs = string(zeros(numNodes,1));
    nodeClasses = strings(numNodes,1);
    for ii = 1:numNodes
        n = nodeList.item(ii-1);
        nodeIDs(ii) = string(char(n.getAttribute('id')));
        nodeClasses(ii) = "Unknown";
        attvalues = n.getElementsByTagName('attvalue');
        for j = 0:attvalues.getLength-1
            att = attvalues.item(j);
            if strcmp(char(att.getAttribute('for')), '0')
                nodeClasses(ii) = string(char(att.getAttribute('value')));
                break;
            end
        end
    end
    idMap = containers.Map(nodeIDs, 1:numNodes);

    % Edges
    edgeList = xDoc.getElementsByTagName('edge');
    numEdges = edgeList.getLength();
    sources = zeros(numEdges,1); targets = zeros(numEdges,1);
    weights = ones(numEdges,1);
    for ii = 1:numEdges
        e = edgeList.item(ii-1);
        srcID = string(char(e.getAttribute('source')));
        tgtID = string(char(e.getAttribute('target')));
        if isKey(idMap, srcID) && isKey(idMap, tgtID)
            sources(ii) = idMap(srcID);
            targets(ii) = idMap(tgtID);
        else
            warning('Edge %d references unknown node(s): %s -> %s', ii, srcID, tgtID);
            continue;
        end
        attvalues = e.getElementsByTagName('attvalue');
        for a = 0:attvalues.getLength-1
            att = attvalues.item(a);
            if strcmp(char(att.getAttribute('for')), '2')
                weights(ii) = str2double(char(att.getAttribute('value')));
                break;
            end
        end
    end

    A_full = sparse(sources, targets, weights, numNodes, numNodes);
    A_full = A_full + A_full';

    uniqueClasses = unique(nodeClasses);
    A_class = struct(); nodeIDs_class = struct();
    for c = 1:numel(uniqueClasses)
        cname = uniqueClasses(c);
        mask = nodeClasses == cname;
        idx = find(mask);
        if isempty(idx), continue; end
        A_sub = A_full(idx, idx);
        deg = full(sum(A_sub,2) + sum(A_sub,1)');
        nz = find(deg > 0);
        if isempty(nz)
            fprintf('Skipping class %s (no internal edges).\n', cname);
            continue;
        end
        A_sub = A_sub(nz, nz);
        nodeIDs_sub = nodeIDs(idx(nz));
        fname = matlab.lang.makeValidName(strcat('class_', char(cname)));
        A_class.(fname) = A_sub;
        nodeIDs_class.(fname) = nodeIDs_sub;
        fprintf('Class %s: %d nodes, %d edges.\n', fname, numel(nodeIDs_sub), nnz(A_sub)/2);
    end
    fprintf('Created adjacency matrices for %d classes.\n', numel(fieldnames(A_class)));
end

function [nodeIDs_out, positions_out] = gexf_create_positions(filename, includeIDs)
% % gexf_extract_positions
% % Extracts node class information from a GEXF file and assigns new
% % 2D positions clustered by year and class.
% %
% % Inputs:
% %   filename   - string, path to GEXF file
% %   includeIDs - string array or cell array of node IDs to include
% %
% % Outputs:
% %   nodeIDs_out   - M×1 string array of node IDs included
% %   positions_out - M×3 numeric array of [x, y, z] positions
% 
%     % Convert includeIDs to string array if cell
%     if iscell(includeIDs)
%         includeIDs = string(includeIDs);
%     end
% 
%     % Read XML
%     xDoc = xmlread(filename);
%     nodeList = xDoc.getElementsByTagName('node');
%     numNodes = nodeList.getLength;
% 
%     % Initialize outputs
%     nodeIDs_out = [];
%     positions_out = [];
% 
%     % Define all possible classes
%     uniqueClasses = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];
% 
%     % Define cluster centers for each class
%     numYears = 5;
%     ySpacing = 2; % distance between year clusters
%     xSpacing = 7;  % distance between A/B classes
%     centers = containers.Map;
% 
%     for y = 1:numYears
%         centers(uniqueClasses((y-1)*2 + 1)) = [-(xSpacing/2), -((y-1)*ySpacing)];
%         centers(uniqueClasses((y-1)*2 + 2)) = [(xSpacing/2), -((y-1)*ySpacing)];
%     end
%     centers("Teachers") = [0, -(numYears*ySpacing)];
% 
%     rng(1); % reproducibility
% 
%     % Initialize output variables
%     nodeIDs_out = {};
%     positions_out = [];
%     classes_out = {}; % To store class information
% 
%     for i = 0:numNodes-1
%         nodeElem = nodeList.item(i);
%         nodeID = string(char(nodeElem.getAttribute('id')));
% 
%         % Skip if not in includeIDs
%         if ~ismember(nodeID, includeIDs)
%             continue;
%         end
% 
%         % --- Extract class info from attvalue with for="0"
%         attvalues = nodeElem.getElementsByTagName('attvalue');
%         nodeClass = "Teachers"; % default fallback
% 
%         for j = 0:attvalues.getLength-1
%             attElem = attvalues.item(j);
%             attrFor = string(char(attElem.getAttribute('for')));
%             if attrFor == "0"
%                 nodeClass = string(char(attElem.getAttribute('value')));
%                 break;
%             end
%         end
% 
%         % --- Determine base position
%         if isKey(centers, nodeClass)
%             base = centers(nodeClass);
%         else
%             base = [0, 0]; % fallback if unknown
%         end
% 
%         % --- Random jitter around base position
%         jitter = (rand(1,2)-0.5) * 5;  % random jitter in [-1,1]
%         pos2D = base + [0.9 * jitter(1), 0.25 * jitter(2)];
%         z = 0;
% 
%         % --- Store outputs
%         nodeIDs_out{end+1,1} = nodeID;
%         positions_out(end+1,:) = [pos2D, z];
%         classes_out{end+1,1} = nodeClass; % Store the class
% 
%         stringClass = string(positions1.Class);
%         TMask = stringClass == "Teachers";
%         TIdx = find(TMask);
%         Pos(TIdx(1),1:2) = [0,-8];
%         Pos(TIdx(2),1:2) = [0,-3];
%         Pos(TIdx(3),1:2) = [0,-5];
%         Pos(TIdx(4),1:2) = [0,-6];
%         Pos(TIdx(5),1:2) = [0,-2];
%         Pos(TIdx(6),1:2) = [0,0];
%         Pos(TIdx(7),1:2) = [0,-4];
%         Pos(TIdx(8),1:2) = [0,-8];
%         Pos(TIdx(9),1:2) = [0,0];
%         Pos(TIdx(10),1:2) = [0,-2];
%     end
% 
%     % Convert to table for easier filtering
%     positions_out = table(nodeIDs_out, positions_out(:,1), positions_out(:,2), positions_out(:,3), classes_out, ...
%                        'VariableNames', {'NodeID', 'X', 'Y', 'Z', 'Class'});
% 
%     fprintf('Generated clustered positions for %d nodes.\n', numel(nodeIDs_out));
% end
% Extracts node class and assigns clustered 2D positions (with z=0).
% includeIDs: string array of node IDs to include (order matters).
    if iscell(includeIDs), includeIDs = string(includeIDs); end
    xDoc = xmlread(filename);
    nodeList = xDoc.getElementsByTagName('node');
    numNodes = nodeList.getLength();

    % Prepare class centres
    uniqueClasses = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];
    numYears = 5; ySpacing = 2; xSpacing = 7;
    centers = containers.Map;
    for y = 1:numYears
        centers(uniqueClasses((y-1)*2 + 1)) = [-(xSpacing/2), -((y-1)*ySpacing)];
        centers(uniqueClasses((y-1)*2 + 2)) = [(xSpacing/2),  -((y-1)*ySpacing)];
    end
    centers('Teachers') = [0, -(numYears*ySpacing)];
    rng(1);

    nodeIDs_out = {};
    positions = [];
    classes_out = {};

    for ii = 0:numNodes-1
        n = nodeList.item(ii);
        nid = string(char(n.getAttribute('id')));
        if ~ismember(nid, includeIDs), continue; end
        % extract class from attvalue for="0"
        nodeClass = 'Teachers'; % default
        attvalues = n.getElementsByTagName('attvalue');
        for j = 0:attvalues.getLength-1
            att = attvalues.item(j);
            if strcmp(char(att.getAttribute('for')), '0')
                nodeClass = string(char(att.getAttribute('value')));
                break;
            end
        end
        if isKey(centers, nodeClass)
            base = centers(nodeClass);
        else
            base = [0,0];
        end
        jitter = (rand(1,2)-0.5) * 5;
        pos2D = base + [0.9*jitter(1), 0.25*jitter(2)];
        z = 0;
        nodeIDs_out{end+1,1} = nid;
        positions(end+1,:) = [pos2D, z];
        classes_out{end+1,1} = nodeClass;
    end

    classes_out_str = string(classes_out);
    TMask = classes_out_str == "Teachers";
    TIdx = find(TMask);
    positions(TIdx(1),1:2) = [0,-8];
    positions(TIdx(2),1:2) = [0,-3];
    positions(TIdx(3),1:2) = [0,-5];
    positions(TIdx(4),1:2) = [0,-6];
    positions(TIdx(5),1:2) = [0,-2];
    positions(TIdx(6),1:2) = [0,0];
    positions(TIdx(7),1:2) = [0,-4];
    positions(TIdx(8),1:2) = [0,-8];
    positions(TIdx(9),1:2) = [0,0];
    positions(TIdx(10),1:2) = [0,-2];

    positions_out = table(nodeIDs_out, positions(:,1), positions(:,2), positions(:,3), classes_out, ...
        'VariableNames', {'NodeID','X','Y','Z','Class'});
    fprintf('Generated clustered positions for %d nodes.\n', size(positions_out,1));
end

function plot_centrality_colourvary(A, centrality, X, nodeIDs_out, scaled, metaFile)
% Plot graph with node colours based on class (from metadata) and marker
% sizes scaled from centrality. 'scaled' controls maximum marker multiplication.
    if nargin < 6, metaFile = fullfile(dataDir,'metadata.txt'); end

    G = graph(A);
    meta = readtable(metaFile, 'Delimiter','\t', 'ReadVariableNames', false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    meta.ID = string(meta.ID);

    % Find class for each node in nodeIDs_out
    nodeIDs_out = string(nodeIDs_out);
    numNodes = numel(nodeIDs_out);
    nodeClasses = strings(numNodes,1);
    for i = 1:numNodes
        idx = find(meta.ID == nodeIDs_out(i), 1);
        if ~isempty(idx)
            nodeClasses(i) = meta.Class{idx};
        else
            nodeClasses(i) = 'Unknown';
        end
    end

    classOrder = ["1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","Teachers"];
    cmap = lines(numel(classOrder));
    cmap(end,:) = [0,0,0];
    nodeColors = zeros(numNodes,3);
    for i = 1:numel(classOrder)
        nodeColors(nodeClasses==classOrder(i), :) = repmat(cmap(i,:), sum(nodeClasses==classOrder(i)), 1);
    end

    % Marker sizing: consistent scaling across networks
    ms = scale_markers(centrality, 1e-16, 10);
    ms = ms / max(ms) * scaled;

    figure;
    if isempty(X)
        p = plot(G, 'Layout', 'force', 'MarkerSize', ms, 'NodeColor', nodeColors, 'EdgeAlpha', 0.01, 'EdgeColor', [0,0,0],'HandleVisibility','off');
    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), 'MarkerSize', ms, 'NodeColor', nodeColors, 'EdgeAlpha', 0.01, 'EdgeColor', [0,0,0],'HandleVisibility','off');
    end
    axis off; box on; grid on;
    % Legend
    hold on;
    for i = 1:numel(classOrder)
        scatter3(NaN, NaN, NaN, 100, cmap(i,:), 'filled');
    end
    legend(classOrder, 'Location','bestoutside');
end

function plot_centrality_colourvary_years(A, centrality, X, nodeIDs_out, scaled, metaFile)
% Variation that colours nodes by year group rather than class. Uses the
% first character of class label as year.
    if nargin < 6, metaFile = fullfile(dataDir,'metadata.txt'); end
    meta = readtable(metaFile, 'Delimiter','\t', 'ReadVariableNames', false);
    meta.Properties.VariableNames = {'ID','Class','Gender'};
    meta.ID = string(meta.ID);

    % Find class for each node in positions_out
    numNodes = numel(nodeIDs_out);
    nodeYears = strings(numNodes,1);
    
    for i = 1:numNodes
        idx = find(meta.ID == nodeIDs_out(i), 1);
        if ~isempty(idx)
            % if meta.Class(idx) contain
            
            nodeYears(i) = strcat("Year ",meta.Class{idx}(1));
        else
            nodeYears(i) = "Unknown";
        end
    end
    
    % Example (user provides):
    yearOrder = ["Year 1","Year 2","Year 3","Year 4","Year 5","Teachers"];

    cmap = lines(numel(yearOrder));
    cmap(end,:) = [0,0,0];
    nodeColors = zeros(numNodes,3);
    for i = 1:numel(yearOrder)
        nodeColors(nodeYears==yearOrder(i), :) = repmat(cmap(i,:), sum(nodeYears==yearOrder(i)), 1);
    end

    % Marker sizing
    ms = scale_markers(centrality, 1e-16, 10);
    ms = ms / max(ms) * scaled;

    figure;
    if isempty(X)
        plot(graph(A), 'Layout','force', ...
            'MarkerSize', ms, ... % scaling by centrality
            'NodeColor', nodeColors, ...                              % color by class
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0], 'HandleVisibility','off', ...
            'NodeLabel', {});

    else
        plot(graph(A), 'XData', X(:,1), 'YData', X(:,2), ...
            'MarkerSize', ms, ...
            'NodeColor', nodeColors, ...  #[.8,.8,.8], ...                         % color by class
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0], 'HandleVisibility','off', ...
            'NodeLabel', {});
    end
    axis off; box on; grid on;
    hold on;
    for i = 1:numel(yearOrder)
        scatter3(NaN, NaN, NaN, 100, cmap(i,:), 'filled');
    end
    legend(yearOrder, 'Location','bestoutside');
end

function msizes = scale_markers(values, minSize, maxSize)
% SCALE_MARKERS Linearly scales values to a marker size range [minSize,maxSize].
    if nargin < 2, minSize = 0; end
    if nargin < 3, maxSize = 20; end
    v = abs(values(:));
    v = v - min(v);
    if max(v) > 0
        v = v ./ max(v);
    end
    msizes = minSize + v * (maxSize - minSize);
end
