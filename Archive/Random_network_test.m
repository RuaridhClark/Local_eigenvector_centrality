% Random modular network with repeated hub geometry
% -------------------------------------------------
% - One connected 20-node hub structure
% - Same hub repeated 6 times
% - SAME relative node positions used in every hub
% - One random inter-hub edge between each hub pair

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
% Create reusable geometry for 20 nodes

theta = linspace(0,2*pi,nodesPerHub+1)';
theta(end) = [];

radius = 1;

hubX = radius*cos(theta);
hubY = radius*sin(theta);

% Small random perturbation for organic appearance
hubX = hubX + 0.15*randn(nodesPerHub,1);
hubY = hubY + 0.15*randn(nodesPerHub,1);

%% STEP 5: Place 6 hubs in larger arrangement
% Arrange hubs on large circle

bigTheta = linspace(0,2*pi,numHubs+1)';
bigTheta(end) = [];

bigRadius = 8;

centerX = bigRadius*cos(bigTheta);
centerY = bigRadius*sin(bigTheta);

%% STEP 6: Build full node coordinates
X = zeros(totalNodes,1);
Y = zeros(totalNodes,1);
Z = zeros(totalNodes,1);

for h = 1:numHubs
    
    idx = (h-1)*nodesPerHub + (1:nodesPerHub);
    
    % Same relative geometry reused for every hub
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

% %% Optional saves
% writetable(positions_tbl,'node_positions.csv');
% writematrix(A,'adjacency_matrix.csv');

%% Display
disp(Pos(1:10,:));

fprintf('Node positions exported as Pos\n');
fprintf('Positions saved to node_positions.csv\n');

global_centrality = local_eigenvector_centrality(A, Pos, false, 1);
local_centrality_6 = local_eigenvector_centrality(A, Pos, true, 6);
plot_centrality_colourvary(A, global_centrality, Pos, 15);
title('Global')
axis equal;
plot_centrality_colourvary(A, local_centrality_6, Pos, 15);
title('Local (i=6)')
axis equal;

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