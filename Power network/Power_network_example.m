clc; clear all
addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Local_eigenvector_centrality'

% --- Read data ---
edges = readtable('edge_df_anonymized.csv');
vertices = readtable('vertice_price_df_anonymized.csv');

% Normalize variable names
edges.Properties.VariableNames = {'Source','Target'};
vertices.Properties.VariableNames = {'Node','Price','Country'};

% Convert node columns to cellstr (works for string or cell arrays)
if isstring(edges.Source) || ischar(edges.Source)
    edges.Source = cellstr(edges.Source);
end
if isstring(edges.Target) || ischar(edges.Target)
    edges.Target = cellstr(edges.Target);
end
if isstring(vertices.Node) || ischar(vertices.Node)
    vertices.Node = cellstr(vertices.Node);
end
if isstring(vertices.Country) || ischar(vertices.Country)
    vertices.Country = cellstr(vertices.Country);
end

% All unique nodes (preserve as cell array)
allNodes = unique([edges.Source; edges.Target]);

% Create fast name -> index map
n = numel(allNodes);
nodeIndexMap = containers.Map(allNodes, 1:n);

% Sparse adjacency
A = sparse(n, n);

% Price lookup (ensure vertex node names are cellstr)
priceMap = containers.Map(vertices.Node, vertices.Price);

% --- Fill adjacency matrix (weight = price of source) ---
% This uses the same loop structure you wrote but with robust extraction
for i = 1:height(edges)
    srcCell = edges{i, 'Source'};    % may be 1x1 cell
    dstCell = edges{i, 'Target'};    % use i for the correct row
    % Extract strings if necessary
    if iscell(srcCell), src = srcCell{1}; else src = srcCell; end
    if iscell(dstCell), dst = dstCell{1}; else dst = dstCell; end

    % Only add edge if both nodes are known and price available
    if isKey(nodeIndexMap, src) && isKey(nodeIndexMap, dst) && isKey(priceMap, src)
        srcIdx = nodeIndexMap(src);
        dstIdx = nodeIndexMap(dst);
        A(srcIdx, dstIdx) = 1; %priceMap(src);
    end
end

% % --- Identify and remove zero-degree (isolated) nodes ---
% % degree = out-degree + in-degree
% outDeg = sum(A ~= 0, 2);    % column vector (n x 1)
% inDeg  = sum(A ~= 0, 1)';   % column vector (n x 1)
% deg = outDeg + inDeg;
% 
% keepMask = deg > 0;         % logical mask of nodes to keep
% keptNodes = allNodes(keepMask);      % cell of node names kept
% removedNodes = allNodes(~keepMask);  % cell of removed node names
% 
% keptOriginalIdx = find(keepMask);    % original numeric indices of kept nodes
% removedOriginalIdx = find(~keepMask);
% 
% % Reduced adjacency
% A_reduced = A(keepMask, keepMask);

% % Save mapping (optional)
% keptTable = table(keptNodes(:), keptOriginalIdx(:), 'VariableNames', {'Node','OriginalIndex'});
% removedTable = table(removedNodes(:), removedOriginalIdx(:), 'VariableNames', {'Node','OriginalIndex'});
% writetable(keptTable, 'kept_nodes.csv');
% writetable(removedTable, 'removed_nodes.csv');

% --- Create graph (directed) using kept nodes ---
G = digraph(A, string(allNodes));

% --- Compute centrality using your local function (note transpose) ---
% local_eigenvector_centrality expects adjacency of directed edges in a particular orientation
centrality = local_eigenvector_centrality(A, [], 0);  % returns row or column vector

% Ensure centrality is a column vector aligned with G.Nodes
centrality = centrality(:);

% --- Align vertex data (country/price) to the reduced node order ---
% vertices.Node contains the full list; use ismember to map
[tf, loc] = ismember(G.Nodes.Name, string(vertices.Node));
if ~all(tf)
    warning('Some graph nodes have no matching entry in vertices table.');
end
nodeCountry = vertices.Country(loc);
nodePrice   = vertices.Price(loc);

% --- Color and size mapping ---
[uniqueCountries, ~, countryIdx] = unique(nodeCountry);
colors = lines(numel(uniqueCountries));   % colormap for countries

% Node sizes scaled by centrality (tweak as needed)
% nodeSizes = .1 + 100 * normalize(centrality, 'range');  % ensures positive marker sizes
nodeSizes = .1 + 100 * (centrality.^.5)./sum(centrality.^.5);  % ensures positive marker sizes

figure;
plot(G, 'Layout', 'force','NodeCData', countryIdx, ...
         'MarkerSize', nodeSizes, ...
         'EdgeAlpha', 0.05, ...
         'NodeLabel', {});

colormap(colors);
colorbar('Ticks', 1:numel(uniqueCountries), 'TickLabels', uniqueCountries);
title('Network Graph (weights = source price)');


