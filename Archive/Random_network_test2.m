%% ================= CENTRALITY COMPARISONS =================
clear;
clc;
close all;

%% Parameters
nodesPerHub = 20;
numHubs     = 6;
pInternal   = 0.15;

totalNodes = nodesPerHub * numHubs;

%% STEP 1: Create connected base hub adjacency
baseAdj = zeros(nodesPerHub);

perm = randperm(nodesPerHub);

% Spanning tree ensures connectivity
for i = 2:nodesPerHub
    
    n1 = perm(i);
    n2 = perm(randi(i-1));
    
    baseAdj(n1,n2) = 1;
    baseAdj(n2,n1) = 1;
end

% Add extra random edges
extraEdges = rand(nodesPerHub) < pInternal;

extraEdges = triu(extraEdges,1);
extraEdges = extraEdges + extraEdges.';

baseAdj = baseAdj | extraEdges;

baseAdj(1:nodesPerHub+1:end) = 0;

%% STEP 2: Replicate adjacency matrix
A = zeros(totalNodes);

for h = 1:numHubs
    
    idx = (h-1)*nodesPerHub + (1:nodesPerHub);
    A(idx,idx) = baseAdj;
end

%% STEP 3: Add one random inter-hub edge per hub pair
for h1 = 1:numHubs-1
    
    idx1 = (h1-1)*nodesPerHub + (1:nodesPerHub);
    
    for h2 = h1+1:numHubs
        
        idx2 = (h2-1)*nodesPerHub + (1:nodesPerHub);
        
        n1 = idx1(randi(nodesPerHub));
        n2 = idx2(randi(nodesPerHub));
        
        A(n1,n2) = 1;
        A(n2,n1) = 1;
    end
end

%% STEP 4: Define ONE relative hub layout
% Increase internal spread of nodes

theta = linspace(0,2*pi,nodesPerHub+1)';
theta(end) = [];

% Larger hub radius
hubRadius = 2;

hubX = hubRadius*cos(theta);
hubY = hubRadius*sin(theta);

% Slightly larger perturbation
hubX = hubX + 0.35*randn(nodesPerHub,1);
hubY = hubY + 0.35*randn(nodesPerHub,1);

%% STEP 5: Place 6 hubs in larger arrangement
% Bring hubs closer together

bigTheta = linspace(0,2*pi,numHubs+1)';
bigTheta(end) = [];

% Smaller spacing between hubs
bigRadius = 6;

centerX = bigRadius*cos(bigTheta);
centerY = bigRadius*sin(bigTheta);

%% STEP 6: Build full node coordinates
X = zeros(totalNodes,1);
Y = zeros(totalNodes,1);
Z = zeros(totalNodes,1);

for h = 1:numHubs
    
    idx = (h-1)*nodesPerHub + (1:nodesPerHub);
    
    % Reuse same hub geometry
    X(idx) = hubX + centerX(h);
    Y(idx) = hubY + centerY(h);
end

%% STEP 7: Create graph
G = graph(A);

%% STEP 8: Plot using fixed coordinates
figure('Color','w');

p = plot(G,...
    'XData',X,...
    'YData',Y,...
    'MarkerSize',5,...
    'LineWidth',1);

title('Repeated Hub Geometry Network');

%% Color nodes by hub
hubColors = lines(numHubs);

nodeColors = zeros(totalNodes,3);

for h = 1:numHubs
    
    idx = (h-1)*nodesPerHub + (1:nodesPerHub);
    
    nodeColors(idx,:) = repmat(hubColors(h,:),nodesPerHub,1);
end

p.NodeColor = nodeColors;

%% STEP 9: Export positions
positions_tbl = table(X,Y,Z,...
    'VariableNames', {'X','Y','Z'});

Pos = positions_tbl{:, {'X','Y','Z'}};

% Existing measures
global_centrality = local_eigenvector_centrality(A, Pos, false, 1);
local_centrality_6 = local_eigenvector_centrality(A, Pos, true, 6);

%% Graph object
G = graph(A);

%% Standard centralities

% Degree centrality
degree_c = centrality(G,'degree');

% Eigenvector centrality
eigen_c = centrality(G,'eigenvector');

% Betweenness centrality
between_c = centrality(G,'betweenness');

% Closeness centrality
close_c = centrality(G,'closeness');

%% Katz Centrality
% ---------------------------------------------------------
% Katz = (I - alpha*A)^-1 * 1
% alpha must be < 1/lambda_max

lambda_max = max(abs(eigs(sparse(A),1)));

alpha = 0.85 / lambda_max;

I = speye(size(A));

katz_c = (I - alpha*A) \ ones(totalNodes,1);

% Normalize
katz_c = katz_c ./ norm(katz_c);

%% PageRank
% ---------------------------------------------------------
% Convert to directed graph for pagerank

DG = digraph(A);

pagerank_c = centrality(DG,'pagerank');

%% Class Centrality
% ---------------------------------------------------------
% Compute eigenvector centrality INSIDE each 20-node hub only
% Ignores all inter-hub connections

class_centrality = zeros(totalNodes,1);

for h = 1:numHubs
    
    % Nodes belonging to this hub
    idx = (h-1)*nodesPerHub + (1:nodesPerHub);
    
    % Extract isolated hub adjacency matrix
    Ahub = A(idx,idx);
    
    % Graph for isolated hub
    Ghub = graph(Ahub);
    
    % Eigenvector centrality within hub only
    cHub = centrality(Ghub,'eigenvector');
    
    % Normalize within class
    cHub = cHub ./ max(cHub);
    
    % Store back into global vector
    class_centrality(idx) = cHub;
end

%% Normalize all measures for comparison
normalize_vec = @(x) x(:) ./ max(x);

degree_c      = normalize_vec(degree_c);
eigen_c       = normalize_vec(eigen_c);
between_c     = normalize_vec(between_c);
close_c       = normalize_vec(close_c);
katz_c        = normalize_vec(katz_c);
pagerank_c    = normalize_vec(pagerank_c);
global_centrality = normalize_vec(global_centrality);
local_centrality_6 = normalize_vec(local_centrality_6);
class_centrality = normalize_vec(class_centrality);

%% Create comparison table
centrality_tbl = table(...
    (1:totalNodes)',...
    degree_c,...
    eigen_c,...
    between_c,...
    close_c,...
    katz_c,...
    pagerank_c,...
    class_centrality,...
    global_centrality,...
    local_centrality_6,...
    'VariableNames',...
    {'Node',...
     'Degree',...
     'Eigenvector',...
     'Betweenness',...
     'Closeness',...
     'Katz',...
     'PageRank',...
     'ClassEigen',...
     'GlobalLocalEigen',...
     'LocalEigen_i6'});

disp(centrality_tbl(1:10,:));

%% ================= VISUALISATIONS =================

%% Global Eigenvector
plot_centrality_colourvary(A, global_centrality, Pos, 15);
title('Global Eigenvector')
axis equal;

%% Local Eigenvector
plot_centrality_colourvary(A, local_centrality_6, Pos, 15);
title('Local Eigenvector (i=6)')
axis equal;

%% Katz
plot_centrality_colourvary(A, katz_c, Pos, 15);
title('Katz Centrality')
axis equal;

%% PageRank
plot_centrality_colourvary(A, pagerank_c, Pos, 15);
title('PageRank')
axis equal;

%% Betweenness
plot_centrality_colourvary(A, between_c, Pos, 15);
title('Betweenness Centrality')
axis equal;

%% Degree
plot_centrality_colourvary(A, degree_c, Pos, 15);
title('Degree Centrality')
axis equal;

%% Class Eigenvector
plot_centrality_colourvary(A, class_centrality, Pos, 15);
title('Class Eigenvector Centrality')
axis equal;

%% ================= CORRELATION ANALYSIS =================

allCentralities = [...
    degree_c,...
    eigen_c,...
    between_c,...
    close_c,...
    katz_c,...
    pagerank_c,...
    class_centrality,...
    global_centrality,...
    local_centrality_6];

labels = {...
    'Degree',...
    'Eigenvector',...
    'Betweenness',...
    'Closeness',...
    'Katz',...
    'PageRank',...
    'ClassEigen',...
    'GlobalLocalEigen',...
    'LocalEigen_i6'};

corrMat = corrcoef(allCentralities);

figure;
imagesc(corrMat);

% Autumn colormap
colormap(flip(autumn));

% Optional: fix colour scaling to correlation range
% caxis([0 1]);
clim([0.9 1])

colorbar;
axis square;

xticks(1:length(labels));
yticks(1:length(labels));

xticklabels(labels);
yticklabels(labels);

xtickangle(45);

title('Centrality Correlation Matrix');

%% ================= TOP NODES =================

N = 10;

fprintf('\nTop %d Katz Centrality Nodes:\n',N);
[~,idx] = sort(katz_c,'descend');
disp(idx(1:N));

fprintf('\nTop %d PageRank Nodes:\n',N);
[~,idx] = sort(pagerank_c,'descend');
disp(idx(1:N));

fprintf('\nTop %d Betweenness Nodes:\n',N);
[~,idx] = sort(between_c,'descend');
disp(idx(1:N));

function plot_centrality_colourvary(A, centrality, X, scaled)
% Plot graph with node colours based on class (from metadata) and marker
% sizes scaled from centrality. 'scaled' controls maximum marker multiplication.

    G = graph(A);

    % Marker sizing: consistent scaling across networks
    ms = scale_markers(centrality, 6, 35);
    ms = ms / max(ms) * scaled;

    figure;
    if isempty(X)
        p = plot(G, 'Layout', 'force', 'MarkerSize', ms, 'EdgeAlpha', 0.01, 'EdgeColor', [0,0,0],'HandleVisibility','off');
    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), 'MarkerSize', ms, 'EdgeAlpha', 0.01, 'EdgeColor', [0,0,0],'HandleVisibility','off');
    end
    axis off; box on; grid on;
    % Legend
    hold on;
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