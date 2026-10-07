%% main_analysis.m
% Reproducible analysis script for "A local eigenvector centrality"
%
% Usage: run this script from its folder. It will load the csv file,
% compute centralities, and produce plots corresponding with the paper.

clear; clc; close all;

%% --- Add path to helper functions ---
addpath('../');

% CSV file 
csv_file = 'thiers_2011.csv';
metafile = "metadata_HS_2011.txt";

%% --- Load adjacency matrices ---
[A_sparse, A_dense, nodeIDs_nonzero] = csv_to_adjacency(csv_file);
[A_classes, nodeIDs_classes] = csv_to_class_adjacency(csv_file);

%% --- Create clustered positions ---
[nodeIDs_pos, positions_tbl] = csv_create_positions(csv_file, nodeIDs_nonzero);
Pos = positions_tbl{:, {'X','Y','Z'}};

%% --- Adjust positions for teachers (manual tweaks) ---
stringClass = string(positions_tbl.Class);
TMask = stringClass == "teacher";
TIdx = find(TMask);
if length(TIdx)>0
    Pos(TIdx(2),1:2) = [-1,-3];
    Pos(TIdx(7),1:2) = [0.25,-2.5];
end

%% --- Compute local eigenvector centrality ---
global_centrality = local_eigenvector_centrality(A_dense, Pos, true, 1);
LEC = local_eigenvector_centrality(A_dense, Pos, false, 5);
plot_centrality_colourvary(A_dense, LEC, Pos, nodeIDs_nonzero, 15, metafile);
title('Local (i=5)')
axis equal;

%% Compute class-level centralities and assemble ordered vector
fprintf('\n--- Computing class-level centralities (per-class local centrality) ---\n');
classList = string(fieldnames(nodeIDs_classes));
allCentralityData = table('Size', [0 2], 'VariableTypes', {'string','double'}, 'VariableNames', {'NodeID','Centrality'});

numClass = numel(classList);
for k = 1:numClass
    classStructname = classList(k);
    ids_1class = nodeIDs_classes.(classStructname);
    % Find nodeIDs_nonzero index matching ids_1class
    idx_nonzero = ismember(nodeIDs_nonzero, ids_1class);
    A_class = A_classes.(classStructname);

    if k < numClass || length(TIdx)==0
        % Compute class-level local centrality (Imax = 1 -> principal eigenvector of class subgraph)
        [class_centrality, ~] = local_eigenvector_centrality(A_class, Pos(idx_nonzero,:), false, 2);
    else % Teachers nodes given LEC values
        class_centrality = LEC(ismember(nodeIDs_nonzero,nodeIDs_classes.classteacher));
    end

    % Store results (NodeID order corresponds to ids_1class)
    T = table(string(ids_1class(:)), class_centrality(:), 'VariableNames', {'NodeID','Centrality'});
    allCentralityData = [allCentralityData; T];
end

% Ensure Teacher nodes present: use LEC values for missing teacher indices
% Map nodeIDs_nonzero to their centrality entries (fallback to LEC)
[commonIdx, ia, ib] = intersect(string(nodeIDs_nonzero), allCentralityData.NodeID, 'stable');
orderedCentrality = zeros(size(nodeIDs_nonzero));

% Fill known entries
if ~isempty(commonIdx)
    orderedCentrality(ia) = allCentralityData.Centrality(ib);
end

% Final plot of class-based centrality
plot_centrality_colourvary(A_dense, orderedCentrality, Pos, nodeIDs_nonzero, 15, metafile);
title('Class (per-class centrality)'); axis equal;

%% --- Compute PageRank centrality ---
G = graph(A_dense);
PRank_centrality = centrality(G,'pagerank','FollowProbability',0.85,'Importance',G.Edges.Weight);

%% --- Plot centrality scaled nodes ---
plot_centrality_colourvary(A_dense, PRank_centrality, Pos, nodeIDs_nonzero, 15, metafile);
title('PageRank')
axis equal;

%% --- Optimize power for LEC_adjust centrality ---
pow_list = 0.05:0.05:1;
wsd_power_values = zeros(size(pow_list));

for i = 1:length(pow_list)
    pow = pow_list(i);
    LEC_adjust = LEC.^pow;
    % LEC_adjust = LEC_adjust / sum(LEC_adjust);
    LEC_adjust = LEC_adjust ./ norm(LEC_adjust);
    difference_vector = PRank_centrality - LEC_adjust;
    wsd_power_values(i) = norm(difference_vector, 2);
end

[~, idx_min] = min(wsd_power_values);
pow_opt = pow_list(idx_min);

% LEC_adjust = LEC.^pow_opt ./ sum(LEC.^pow_opt);
LEC_adjust = LEC.^pow_opt ./ norm(LEC.^pow_opt);
plot_centrality_colourvary(A_dense, LEC_adjust, Pos, nodeIDs_nonzero, 15, metafile);
% Add title including pow_opt
title(sprintf('Local eigenvector centrality (p=%.2f)', pow_opt));
axis equal;

%% --- Normalization & difference vectors ---
Local_norm = (LEC - median(LEC)) / mad(LEC, 1);
PR_norm = (PRank_centrality - median(PRank_centrality)) / mad(PRank_centrality, 1);
LEC_adjust_norm = (LEC_adjust - median(LEC_adjust)) / mad(LEC_adjust, 1);

% normalize_vec = @(x) (x(:)-median(x(:)))./mad(x(:),1);

validIdx = setdiff(1:length(Local_norm), TIdx);

Local_valid = Local_norm(validIdx);
PR_valid = PR_norm(validIdx);
LEC_adjust_valid = LEC_adjust_norm(validIdx);

diff_PRank = Local_valid - PR_valid;
diff_LEC_adjust = LEC_adjust_valid - PR_valid;

[~, sortIdx] = sort(Local_valid);
Local_sorted = Local_valid(sortIdx);
PR_sorted = PR_valid(sortIdx);
Class_sorted = orderedCentrality(sortIdx);
LEC_adjust_sorted = LEC_adjust_valid(sortIdx);
x = (1:length(Local_sorted))';

%% --- Euclidean norms ---
EN_local = norm(PR_sorted - Local_sorted, 2);
EN_power = norm(PR_sorted - LEC_adjust_sorted, 2);
fprintf('Euclidean norm (LEC - PR): %.4f\n', EN_local);
fprintf('Euclidean norm (LEC_adjust - PR): %.4f\n', EN_power);

%% Area plots
% area_plot(Local_sorted,Class_sorted,LEC_adjust_sorted,PR_sorted,"Local","Class",sprintf('Local (p=%1.2f)', pow_opt),"PageRank")
area_plot(Local_sorted,PR_sorted,LEC_adjust_sorted,PR_sorted,"Local","PageRank",sprintf('LEC (p=%1.2f)', pow_opt),"PR")

% --- Boxplots ---
nexttile(7,[1 3]); hold on;
all_sets = [diff_PRank; diff_LEC_adjust];
grps = [repmat({'LEC − PR'}, numel(diff_PRank), 1);
        repmat({sprintf('LEC (p=%1.2f) − PR', pow_opt)}, numel(diff_LEC_adjust), 1)];

uniqueGroups = unique(grps, 'stable'); numGroups = numel(uniqueGroups);
positions = 1:numGroups;

colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980];
alpha = 0.3; lighterColors = colors + (1 - colors)*alpha;

for i = 1:numGroups
    idx = strcmp(grps, uniqueGroups{i});
    bc = boxchart(positions(i)*ones(sum(idx),1), all_sets(idx), 'BoxWidth',0.7);
    bc.BoxFaceColor = lighterColors(i,:);
    bc.LineWidth = 1.5; bc.MarkerColor = lighterColors(i,:);
end

[~,~,ic] = unique(grps,'stable');
for i = 1:numGroups
    Ind = find(ic==i);
    scatter(positions(i)+randn(size(Ind))*0.1, all_sets(Ind), 25, colors(i,:), 'filled', ...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.3);
end

ylabel('Normalised Centrality Difference');
set(gca,'XTick',positions,'XTickLabel',uniqueGroups); grid off;

%% --- End of main_analysis.m ---

%% ----------------------
% Local helper functions
%% ----------------------

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
    
    % --- Area between LEC_adjust and PR ---
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
    fprintf('Comparison between %s: t-test p-value = %.6f\n', nametitle, p);
    
    figure;
    boxchart(g, y)
    ylabel(nametitle)
end


function [A_sparse, A_dense, nodeIDs_nonzero] = csv_to_adjacency(filename)
    contactData = readtable(filename,'FileType','text','Delimiter','\t','ReadVariableNames',false);
    contactData.Properties.VariableNames = {'t','i','j','Ci','Cj'};
    
    nodeIDs = unique([contactData.i; contactData.j]);
    numNodes = numel(nodeIDs);
    
    [~, sources] = ismember(contactData.i, nodeIDs);
    [~, targets] = ismember(contactData.j, nodeIDs);
    weights = ones(size(sources));
    
    A_sparse = sparse(sources, targets, weights, numNodes, numNodes);
    A_sparse = A_sparse + A_sparse';
    A_sparse(A_sparse>0 & A_sparse<1)=1;
    
    deg = full(sum(A_sparse,2)+sum(A_sparse,1)');
    nonzero_idx = find(deg>0);
    nodeIDs_nonzero = nodeIDs(nonzero_idx);
    
    A_dense = full(A_sparse(nonzero_idx, nonzero_idx));
end

function [A_class, nodeIDs_classes] = csv_to_class_adjacency(filename)
    contactData = readtable(filename,'FileType','text','Delimiter','\t','ReadVariableNames',false);
    contactData.Properties.VariableNames = {'t','i','j','Ci','Cj'};
    nodeIDs = unique([contactData.i; contactData.j]);
    numNodes = numel(nodeIDs);
    [~, sources] = ismember(contactData.i, nodeIDs);
    [~, targets] = ismember(contactData.j, nodeIDs);
    weights = ones(size(sources));
    
    A_full = sparse(sources, targets, weights, numNodes, numNodes);
    A_full = A_full + A_full';
    A_full(A_full>0 & A_full<1)=1;
    
    classMap = containers.Map('KeyType','double','ValueType','char');
    for r=1:height(contactData)
        i=contactData.i(r); j=contactData.j(r);
        Ci=string(contactData.Ci(r)); Cj=string(contactData.Cj(r));
        if ~isKey(classMap,i), classMap(i)=Ci; end
        if ~isKey(classMap,j), classMap(j)=Cj; end
    end
    
    nodeClasses = strings(numNodes,1);
    for n=1:numNodes
        id=nodeIDs(n);
        if isKey(classMap,id), nodeClasses(n)=classMap(id); else nodeClasses(n)="Unknown"; end
    end
    
    uniqueClasses = unique(nodeClasses);
    A_class = struct(); nodeIDs_classes = struct();
    for c=1:numel(uniqueClasses)
        className = uniqueClasses(c);
        if className=="Unknown", continue; end
        classMask = nodeClasses==className; classIdx=find(classMask);
        if isempty(classIdx), continue; end
        A_sub = A_full(classIdx,classIdx);
        deg = full(sum(A_sub,2)+sum(A_sub,1)');
        nonzeroIdx = find(deg>0);
        if isempty(nonzeroIdx), continue; end
        A_sub = A_sub(nonzeroIdx,nonzeroIdx);
        nodeIDs_sub = nodeIDs(classIdx(nonzeroIdx));
        classStructname = matlab.lang.makeValidName("class"+className);
        A_class.(classStructname)=A_sub;
        nodeIDs_classes.(classStructname)=nodeIDs_sub;
    end
end

function [nodeIDs_out, positions_out] = csv_create_positions(filename, includeIDs)
    if iscell(includeIDs), includeIDs=string(includeIDs); end
    contactData = readtable(filename,'FileType','text','Delimiter','\t','ReadVariableNames',false);
    contactData.Properties.VariableNames={'t','i','j','Ci','Cj'};
    nodeList = unique([contactData.i; contactData.j]);
    classMap = containers.Map('KeyType','double','ValueType','char');
    for r=1:height(contactData)
        i=contactData.i(r); j=contactData.j(r);
        Ci=string(contactData.Ci(r)); Cj=string(contactData.Cj(r));
        if ~isKey(classMap,i), classMap(i)=Ci; end
        if ~isKey(classMap,j), classMap(j)=Cj; end
    end
    uniqueClasses = unique(string(values(classMap)))';
    xSpacing=9; ySpacing=6;
    centers = containers.Map;
    for c=1:numel(uniqueClasses)
        if c~=4, xPos=mod(c-1,3)*xSpacing - xSpacing; yPos=mod(c-1,2)*ySpacing;
        else xPos=0; yPos=-ySpacing; end
        centers(uniqueClasses(c))=[xPos, yPos];
    end
    rng(1);
    nodeIDs_out={}; positions_out=[]; classes_out={};
    for n=1:numel(nodeList)
        nodeID=nodeList(n);
        if ~ismember(nodeID, includeIDs), continue; end
        if isKey(classMap,double(nodeList(n))), nodeClass=string(classMap(double(nodeList(n))));
        else nodeClass="Unknown"; end
        if isKey(centers,nodeClass), base=centers(nodeClass); else base=[0,0]; end
        jitter=(rand(1,2)-0.5)*5; pos2D=base + [0.9*jitter(1),0.9*jitter(2)]; z=0;
        nodeIDs_out{end+1,1}=nodeID;
        positions_out(end+1,:)=[pos2D,z];
        classes_out{end+1,1}=nodeClass;
    end
    positions_out=table(nodeIDs_out,positions_out(:,1),positions_out(:,2),positions_out(:,3),classes_out, ...
        'VariableNames',{'NodeID','X','Y','Z','Class'});
end

function plot_centrality_colourvary(A, centrality, X, nodeIDs_out, scaled, metafile)
% Plot graph with node colours based on class (from metadata) and marker
% sizes scaled from centrality. 'scaled' controls maximum marker multiplication.
    if nargin == 6
        %% Load metadata
        meta = readtable(metafile, 'Delimiter','\t','ReadVariableNames',false);
        meta.Properties.VariableNames = {'ID','Class','Gender'};
    end

    G = graph(A);

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

    classOrder = unique(nodeClasses, 'stable');
    cmap = lines(numel(classOrder));
    nodeColors = zeros(numNodes,3);
    for i = 1:numel(classOrder)
        nodeColors(nodeClasses==classOrder(i), :) = repmat(cmap(i,:), sum(nodeClasses==classOrder(i)), 1);
    end

    % Marker sizing: consistent scaling across networks
    ms = scale_markers(centrality, 6, 35);
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

function msizes = scale_markers(values, minSize, maxSize)
% SCALE_MARKERS Linearly scales values to a marker size range [minSize,maxSize].
    if nargin < 2, minSize = 6; end
    if nargin < 3, maxSize = 20; end
    v = abs(values(:));
    v = v - min(v);
    if max(v) > 0
        v = v ./ max(v);
    end
    msizes = minSize + v * (maxSize - minSize);
end