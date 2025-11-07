clc; clear; close all;

% Parameters
N = 50;  % Number of nodes
k = 5;    % Number of nearest neighbors

% Generate random 2D points
X = rand(N,2);  % Random (x,y) coordinates for nodes

% Compute pairwise Euclidean distances
Dist = pdist2(X, X);

% Sort distances and get k-nearest neighbors for each node
[~, idx] = sort(Dist, 2, 'ascend'); % Sort distances row-wise

% Build adjacency matrix
A = zeros(N, N);
for i = 1:N
    A(i, idx(i, 2:k+1)) = 1;  % Connect to k nearest (excluding itself)
end

[centrality] = local_eigenvector_centrality(A,X,1);