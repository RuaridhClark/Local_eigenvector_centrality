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
gexfFile = 'sp_data_school_day_1_g.gexf';
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

% rowSum = sum(A_dense, 2);
% A_dense = A_dense ./ rowSum;
% A_dense(rowSum == 0, :) = 0;   % handle isolated nodes

% Create a graph for layout and PageRank
G_dense = digraph(A_dense);

%% Generate (clustered) node positions for plotting
fprintf('\n--- Generating node positions (clustered by class/year) ---\n');
[nodeIDs_pos, positions_table, TIdx] = gexf_create_positions(gexfFile, nodeIDs_nonzero);
Pos = positions_table{:, {'X', 'Y', 'Z'}};    % Nx3 numeric

%% Create class-separated adjacency matrices
fprintf('\n--- Building class adjacency matrices ---\n');
[A_classes, nodeIDs_class] = gexf_to_class_adjacency(gexfFile);

%% Compute PageRank centrality (for comparison)
fprintf('\n--- Computing PageRank centrality ---\n');
if any(G_dense.Edges.Weight)
    PRank_centrality = centrality(G_dense, 'pagerank', 'FollowProbability', .85, 'Importance', G_dense.Edges.Weight);
else
    PRank_centrality = centrality(G_dense, 'pagerank', 'FollowProbability', .85);
end

%% Compute local eigenvector centrality for whole dense network (example Imax choices)
fprintf('\n--- Computing local eigenvector centralities (Imax = 5 and 10) ---\n');
[EC,  ~] = local_eigenvector_centrality(A_dense', Pos, true, 1);
close([1,3,4,5,6])
[LEC_5,  ~] = local_eigenvector_centrality(A_dense', Pos, false, 5);
[LEC_10, ~] = local_eigenvector_centrality(A_dense', Pos, false, 10);

%% Plot PageRank and Local centrality maps coloured by class
fprintf('\n--- Plotting centrality maps ---\n');

plot_centrality_colourvary_years(A_dense, EC, Pos, nodeIDs_nonzero, 20, metaFile);
title('Global'); axis equal;

plot_centrality_colourvary_years(A_dense, LEC_5, Pos, nodeIDs_nonzero, 20, metaFile);
title('Local (i=5)'); axis equal;

plot_centrality_colourvary(A_dense, LEC_10, Pos, nodeIDs_nonzero, 20, metaFile);
title('Local (i=10)'); axis equal;

plot_centrality_colourvary(A_dense, PRank_centrality, Pos, nodeIDs_nonzero, 15, metaFile);
title('PageRank'); axis equal;

%% Compute class-level centralities and assemble ordered vector
fprintf('\n--- Computing class-level centralities (per-class local centrality) ---\n');
classList = string(fieldnames(nodeIDs_class));
allCentralityData = table('Size', [0 2], 'VariableTypes', {'string','double'}, 'VariableNames', {'NodeID','Centrality'});

numClass = numel(classList);
classCorrs = zeros(numClass,1);
sumClass = zeros(numClass,1);
for k = 1:numClass
    classStructname = classList(k);
    ids_1class = nodeIDs_class.(classStructname);
    A_class = A_classes.(classStructname);

    if k < numClass
        % Compute class-level local centrality (Imax = 1 -> principal eigenvector of class subgraph)
        [class_centrality, ~] = local_eigenvector_centrality(A_class, Pos(ismember(string(positions_table.NodeID), ids_1class),:), false, 1);
    else % Teachers nodes given LEC_10 values
        class_centrality = LEC_10(ismember(nodeIDs_nonzero,nodeIDs_class.class_Teachers));
    end

    display(classStructname)
    CC = corr([class_centrality,LEC_10(ismember(string(positions_table.NodeID), ids_1class))], 'Rows','pairwise', 'Type', 'Pearson');
    classCorrs(k) = CC(1,2);
    sumCorrs(k) = sum(A_class(:));

    % Store results (NodeID order corresponds to ids_1class)
    T = table(string(ids_1class(:)), class_centrality(:), 'VariableNames', {'NodeID','Centrality'});
    allCentralityData = [allCentralityData; T];
end

corr([classCorrs(1:end-1),full(sumCorrs(1:end-1))'], 'Rows','pairwise', 'Type', 'Pearson')
corr([classCorrs(1:end-1),full(sumCorrs(1:end-1))'], 'Rows','pairwise', 'Type', 'Spearman')

% Ensure Teacher nodes present: use LEC_10 values for missing teacher indices
% Map nodeIDs_nonzero to their centrality entries (fallback to LEC_10)
[commonIdx, ia, ib] = intersect(nodeIDs_nonzero, allCentralityData.NodeID, 'stable');
KCEC = zeros(size(nodeIDs_nonzero));

% Fill known entries
if ~isempty(commonIdx)
    KCEC(ia) = allCentralityData.Centrality(ib);
end

% Final plot of class-based centrality
plot_centrality_colourvary(A_dense, KCEC, Pos, nodeIDs_nonzero, 15, metaFile);
title('Class (per-class centrality)'); axis equal;

%% Area plot
%% --- Normalization & difference vectors ---

normalize_vec = @(x) x(:)./sum(x);

LEC_norm = normalize_vec(LEC_10);
PR_norm = normalize_vec(PRank_centrality);
KCEC_norm = normalize_vec(KCEC);
EC_norm = normalize_vec(EC);

validIdx = setdiff(1:length(LEC_norm), TIdx);

LEC_valid = LEC_norm(validIdx);
PR_valid = PR_norm(validIdx);
KCEC_valid = KCEC_norm(validIdx);

diff_Class = LEC_valid - KCEC_valid;
diff_PRank = LEC_valid - PR_valid;

[~, sortIdx] = sort(LEC_valid);
LEC_sorted = LEC_valid(sortIdx);
PR_sorted = PR_valid(sortIdx);
KCEC_sorted = KCEC_valid(sortIdx);
x = (1:length(LEC_sorted))';

area_plot(LEC_sorted,KCEC_sorted,LEC_sorted,PR_sorted,"Local","Class","Local","PageRank")


%% Nuanced correlation analysis: split by dominance
fprintf('\n--- Computing class-wise correlation matrices with dominance split ---\n');

centralityNames = { ...
    'Global', ...
    'Local (i=5)', ...
    'Local (i=10)', ...
    'PageRank', ...
    'Class-based'};

centralityAll = [ ...
    EC, ...
    LEC_5(:), ...
    LEC_10(:), ...
    PRank_centrality(:), ...
    KCEC(:)];

classNames = string(fieldnames(nodeIDs_class));
numMeasures = size(centralityAll,2);

% Containers for the three conditions
Corr_All   = struct();
Corr_HighC = struct();  % KCEC > LEC_10
Corr_LowC  = struct();  % KCEC <= LEC_10

for k = 1:numel(classNames)

    cname = classNames(k);
    classNodeIDs = nodeIDs_class.(cname);

    % Map to dense indices
    [~, idxDense] = ismember(classNodeIDs, nodeIDs_nonzero);
    idxDense = idxDense(idxDense > 0);

    if numel(idxDense) < 3
        warning('Skipping class %s (too few nodes).', cname);
        continue;
    end

    % Extract data
    C = centralityAll(idxDense, :);
    ordC = KCEC(idxDense);
    locC = LEC_10(idxDense);

    % --- Masks
    maskHigh = ordC > locC;
    maskLow  = ~maskHigh;

    % --- All nodes
    Corr_All.(cname) = corr(C, 'Rows','pairwise', 'Type', 'Pearson');

    % --- High KCEC set
    if sum(maskHigh) > 0
        Corr_HighC.(cname) = corr(C(maskHigh,:), 'Rows','pairwise', 'Type', 'Spearman');
    end

    % --- Low KCEC set
    if sum(maskLow) > 0
        Corr_LowC.(cname) = corr(C(maskLow,:), 'Rows','pairwise', 'Type', 'Spearman');
    end

    fprintf('%s: %d total | %d high-ordered | %d low-ordered\n', ...
        cname, numel(idxDense), sum(maskHigh), sum(maskLow));
end

% Split into High vs Low 
maskHigh = (KCEC_norm > LEC_norm);
maskLow  = (KCEC_norm < LEC_norm);

C = centralityAll;
Corr_1 = corr(C, 'Rows','pairwise', 'Type', 'Pearson');
Corr_1_HighC = corr(C(maskHigh,:), 'Rows','pairwise', 'Type', 'Spearman');
Corr_1_LowC = corr(C(maskLow,:), 'Rows','pairwise', 'Type', 'Spearman');

Corr_Global = struct();
Corr_Global.All = Corr_1 + diag(NaN(size(Corr_1_HighC,1),1));
Corr_Global.Community_centric = Corr_1_HighC + diag(NaN(size(Corr_1_HighC,1),1));
Corr_Global.Global_centric  = Corr_1_LowC + diag(NaN(size(Corr_1_LowC,1),1));
centralityShortlist = { ...
    'Global', ...
    'Local (i=5)', ...
    'Local (i=10)', ...
    'PageRank', ...
    'Class-based'};

plot_corr_allplots( ...
    Corr_Global, ...
    centralityShortlist, ...
    'Global Centrality Correlations');

% create a boxchart
boxplot_comparison(EC,~maskLow,maskLow,"Global")
boxplot_comparison(KCEC,~maskLow,maskLow,"Class")
boxplot_comparison(PRank_centrality,~maskLow,maskLow,"PageRank")

% Correlation analysis 
LEC_diff = LEC_norm - KCEC_norm;
EC_diff = EC_norm - KCEC_norm;
PR_diff = PR_norm - KCEC_norm;

% Assess correlations between diff vectors
EC_corr = corr(LEC_diff,EC_diff,'Type','Spearman');
PR_corr = corr(LEC_diff,PR_diff,'Type','Spearman');

%% ----------------------
% Local helper functions
%% ----------------------
function [] = area_plot(metricA,metricB,metric1,metric2,nameA,nameB,name1,name2)
    x = (1:length(metric1))';
    figure;
    
    % --- Area between Local and PR ---
    t = tiledlayout(1,9, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Difference area plots on left
    nexttile(1,[1 6]); hold on;
    X_fill = [x; flipud(x)];
    Y_fill = [metricA; flipud(metricB)];
    fill(X_fill, Y_fill, [0 0.4470 0.7410], 'FaceAlpha',0.3,'EdgeColor','none');
    
    % --- Area between ---
    Y_fill = [metric1; flipud(metric2)];
    fill(X_fill, Y_fill, [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,'EdgeColor','none');
    
    legend(""+nameA+" − "+nameB+" area", ""+name1+" − "+name2+" area",'Location','northwest');
    ylabel('Normalised Centrality'); xlabel('Sorted Index'); grid off; 
    % xlim([0 120]);
end

function boxplot_comparison(metric,maskHigh,maskLow,nametitle)
    % create a boxchart of Global centrality
    yHigh = metric(maskHigh);
    yLow  = metric(maskLow);
    
    y = [yHigh; yLow];
    g = [repmat("Community-centric", numel(yHigh), 1); repmat("Globally-connected",  numel(yLow),  1)];
    g = categorical(g);    % optional

    % Perform a statistical test to compare the two distributions
    [h, p] = ttest2(yHigh, yLow);
    fprintf('Comparison between %s: t-test p-value = %.3e\n', nametitle, p);
    
    figure;
    boxchart(g, y)
    ylabel(nametitle)
end

function plot_corr_subplots(CorrStruct, centralityNames, figTitle)

    classFields = fieldnames(CorrStruct);
    nClass = numel(classFields);
    nMeasures = numel(centralityNames);

    nCols = ceil(sqrt(nClass));
    nRows = ceil(nClass / nCols);

    figure('Units','normalized','Position',[0.05 0.05 0.85 0.8]);
    tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

    for k = 1:nClass
        cname = classFields{k};
        R = CorrStruct.(cname);

        nexttile;
        imagesc(R);
        axis square;
        caxis([0.5 1]);
        colormap(flipud(autumn));

        title(strrep(cname,'_',' '), 'Interpreter','none');

        set(gca, ...
            'XTick', 1:nMeasures, ...
            'YTick', 1:nMeasures, ...
            'XTickLabel', centralityNames, ...
            'YTickLabel', centralityNames, ...
            'XTickLabelRotation', 45, ...
            'FontSize', 9);

        if mod(k-1, nCols) ~= 0
            set(gca, 'YTickLabel', []);
        end
        if k <= (nRows-1)*nCols
            set(gca, 'XTickLabel', []);
        end
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Correlation coefficient';

    sgtitle(figTitle, 'FontWeight','bold');
end

function plot_corr_allplots(CorrStruct, centralityNames, figTitle)
% Plots exactly two correlation matrices (e.g. High vs Low dominance)

    setNames = fieldnames(CorrStruct);
    nSets = numel(setNames);
    assert(nSets == 3, 'CorrStruct must contain exactly three fields.');

    nMeasures = numel(centralityNames);

    figure('Units','normalized','Position',[0.15 0.25 0.7 0.45]);
    tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

    % Build colormap (blue -> white -> red)
    n = 256;
    n1 = floor(n/2);
    blue = [59, 76,192]/255;   % tweak colors if desired
    red  = [180,  4, 38]/255;
    white = [1 1 1];

    cmap = [
        interp1([0 1],[blue;white], linspace(0,1,n1), 'linear');
        interp1([0 1],[white;red],  linspace(0,1,n-n1), 'linear')
    ];

    for k = 1:3
        setName = setNames{k};
        R = CorrStruct.(setName);

        nexttile;
        h = imagesc(R);
        axis square;
        caxis([-1 1]);              % preserve your chosen range
        colormap(cmap)

        % ---- mask NaNs ----
        h.AlphaData = ~isnan(R);      % NaNs become transparent
        set(gca, 'Color', [1 1 1]);   % background color for NaNs

        title(strrep(setName,'_','-'), 'Interpreter','none');

        set(gca, ...
            'XTick', 1:nMeasures, ...
            'YTick', 1:nMeasures, ...
            'XTickLabel', centralityNames, ...
            'YTickLabel', centralityNames, ...
            'XTickLabelRotation', 45, ...
            'FontSize', 10);
    end

    % Shared colorbar
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'Correlation coefficient';

    sgtitle(figTitle, 'FontWeight','bold');
end


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

function [nodeIDs_out, positions_out, TIdx] = gexf_create_positions(filename, includeIDs)
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

    G = digraph(A);
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
    % ms = ms ./ sum(ms) * scaled; %ms / max(ms) * scaled;

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
        plot(digraph(A), 'Layout','force', ...
            'MarkerSize', ms, ... % scaling by centrality
            'NodeColor', nodeColors, ...                              % color by class
            'EdgeAlpha', 0.01, 'EdgeColor', [0 0 0], 'HandleVisibility','off', ...
            'NodeLabel', {});

    else
        plot(digraph(A), 'XData', X(:,1), 'YData', X(:,2), ...
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
