%% Contact Network Construction from data.dat
% File format: t   i   j   Si   Sj
% Each line represents a contact active in the 20s interval [t-20, t]

% Load the data file (tab-separated)
filename = 'detailed_list_of_contacts_Hospital.dat';

% Read the tab-delimited file: numeric + string columns
opts = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
opts = setvartype(opts, {'Var4','Var5'}, 'string'); % ensure status vars are strings
T = readtable(filename, opts);

% Extract variables
t  = T{:,1};  % time in seconds
i  = T{:,2};  % person i
j  = T{:,3};  % person j
Si = T{:,4};  % status of i
Sj = T{:,5};  % status of j

%% Build weighted contact network
% Each contact lasts 20 seconds per line entry
contact_duration = 20; % seconds

%% Reindex node IDs to avoid empty nodes
% Extract all unique IDs appearing in i or j
all_nodes = unique([i; j]);
num_nodes = numel(all_nodes);

% Create mapping from original IDs to consecutive integers (1..num_nodes)
[~, i_new] = ismember(i, all_nodes);
[~, j_new] = ismember(j, all_nodes);

% Combine into pairs and sum contact times as before
pairs = [i_new, j_new];
pairs = sort(pairs, 2); % ensure undirected

contact_duration = 20; % seconds per contact line
PairTable = table(pairs(:,1), pairs(:,2), repmat(contact_duration, size(pairs,1), 1), ...
    'VariableNames', {'Node1','Node2','Time'});

% Aggregate total contact time per unique pair
Gsummary = groupsummary(PairTable, {'Node1','Node2'}, 'sum', 'Time');

%% Create graph using reindexed IDs
G = graph(Gsummary.Node1, Gsummary.Node2, Gsummary.sum_Time);

% Assign node labels as the original IDs
G.Nodes.OriginalID = all_nodes;

%% (Optional) Store node status info
% Build a map of node IDs -> status
all_nodes = unique([i; j]);
status_map = containers.Map('KeyType','double','ValueType','char');

for n = 1:numel(all_nodes)
    idx = find(i == all_nodes(n) | j == all_nodes(n), 1, 'first');
    if ~isempty(idx)
        if i(idx) == all_nodes(n)
            status_map(all_nodes(n)) = Si(idx);
        else
            status_map(all_nodes(n)) = Sj(idx);
        end
    end
end

% Example: display status of first 5 nodes
disp('Example node status mapping:');
for n = 1:min(5, numel(all_nodes))
    fprintf('Node %d -> %s\n', all_nodes(n), status_map(all_nodes(n)));
end

%% Create weighted adjacency matrix (only connected nodes)

% Compute the adjacency matrix (weighted)
A = adjacency(G, 'weighted');

% Identify nodes with at least one connection (nonzero degree)
deg = sum(A > 0, 2);      % degree vector
connected_nodes = find(deg > 0);

% Extract submatrix with only connected nodes
A_connected = A(connected_nodes, connected_nodes);

% %% Display results
% fprintf('Total nodes before filtering: %d\n', numnodes(G));
% fprintf('Nodes with at least one connection: %d\n', numel(connected_nodes));
% 
% % Show small part of adjacency matrix
% disp('Weighted adjacency matrix (subset):');
% disp(A_connected(1:min(5,end), 1:min(5,end)));
% 
% %% (Optional) Save to file
% writematrix(A_connected, 'adjacency_connected.csv');

[centrality] = local_eigenvector_centrality(A_connected,[],1);