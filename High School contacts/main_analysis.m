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
local_centrality_5 = local_eigenvector_centrality(A_dense, Pos, false, 5);
plot_centrality_colourvary(A_dense, local_centrality_5, Pos, nodeIDs_nonzero, 15, metafile);
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
    A_class = A_classes.(classStructname);

    if k < numClass || length(TIdx)==0
        % Compute class-level local centrality (Imax = 1 -> principal eigenvector of class subgraph)
        [class_centrality, ~] = local_eigenvector_centrality(A_class, [], false, 1);
    else % Teachers nodes given local_centrality_10 values
        class_centrality = local_centrality_5(ismember(nodeIDs_nonzero,nodeIDs_classes.classteacher));
    end

    % Store results (NodeID order corresponds to ids_1class)
    T = table(string(ids_1class(:)), class_centrality(:), 'VariableNames', {'NodeID','Centrality'});
    allCentralityData = [allCentralityData; T];
end

% Ensure Teacher nodes present: use local_centrality_10 values for missing teacher indices
% Map nodeIDs_nonzero to their centrality entries (fallback to local_centrality_10)
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

%% --- Optimize power for warped centrality ---
pow_list = 0.05:0.05:1;
wsd_power_values = zeros(size(pow_list));

for i = 1:length(pow_list)
    pow = pow_list(i);
    warped = local_centrality_5.^pow;
    warped = warped / sum(warped);
    difference_vector = PRank_centrality - warped;
    wsd_power_values(i) = norm(difference_vector, 2);
end

[~, idx_min] = min(wsd_power_values);
pow_opt = pow_list(idx_min);

warped = local_centrality_5.^pow_opt ./ sum(local_centrality_5.^pow_opt);
plot_centrality_colourvary(A_dense, warped, Pos, nodeIDs_nonzero, 15, metafile);
% Add title including pow_opt
title(sprintf('Local eigenvector centrality (p=%.2f)', pow_opt));
axis equal;

%% --- Normalization & difference vectors ---
Local_norm = (local_centrality_5 - median(local_centrality_5)) / mad(local_centrality_5, 1);
PR_norm = (PRank_centrality - median(PRank_centrality)) / mad(PRank_centrality, 1);
warped_norm = (warped - median(warped)) / mad(warped, 1);

validIdx = setdiff(1:length(Local_norm), TIdx);

Local_valid = Local_norm(validIdx);
PR_valid = PR_norm(validIdx);
warped_valid = warped_norm(validIdx);

diff_PRank = Local_valid - PR_valid;
diff_warped = warped_valid - PR_valid;

[~, sortIdx] = sort(Local_valid);
Local_sorted = Local_valid(sortIdx);
PR_sorted = PR_valid(sortIdx);
warped_sorted = warped_valid(sortIdx);
x = (1:length(Local_sorted))';

%% --- Euclidean norms ---
EN_local = norm(PR_sorted - Local_sorted, 2);
EN_power = norm(PR_sorted - warped_sorted, 2);
fprintf('Euclidean norm (Local - PageRank): %.4f\n', EN_local);
fprintf('Euclidean norm (Warped - PageRank): %.4f\n', EN_power);

%% --- Area plots and boxplots ---
figure;

% --- Area between Local and PR ---
t = tiledlayout(1,9, 'TileSpacing', 'compact', 'Padding', 'compact');

% Difference area plots on left
nexttile(1,[1 6]); hold on;
X_fill = [x; flipud(x)];
Y_fill = [Local_sorted; flipud(PR_sorted)];
fill(X_fill, Y_fill, [0 0.4470 0.7410], 'FaceAlpha',0.3,'EdgeColor','none');

% --- Area between Warped and PR ---
Y_fill = [warped_sorted; flipud(PR_sorted)];
fill(X_fill, Y_fill, [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,'EdgeColor','none');

legend('Local − PageRank area', sprintf('Local (p=%1.2f) − PageRank area', pow_opt),'Location','northwest');
ylabel('Normalised Centrality'); xlabel('Sorted Index'); grid off; xlim([0 120]);

% --- Boxplots ---
nexttile(7,[1 3]); hold on;
all_sets = [diff_PRank; diff_warped];
grps = [repmat({'Local − PageRank'}, numel(diff_PRank), 1);
        repmat({sprintf('Local (p=%1.2f) − PageRank', pow_opt)}, numel(diff_warped), 1)];

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

%% =======================
%% Helper functions below
%% =======================

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

function [A_class, nodeIDs_class] = csv_to_class_adjacency(filename)
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
    A_class = struct(); nodeIDs_class = struct();
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
        nodeIDs_class.(classStructname)=nodeIDs_sub;
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
    xSpacing=7; ySpacing=4;
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