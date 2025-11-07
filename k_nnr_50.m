load('adj_50_1.mat')

X = [x_pos, y_pos];
A = adj;

[centrality] = local_eigenvector_centrality(A,X,1);

% Standardize rows (optional but common)
A_std = zscore(A, 0, 2); % Z-score along each row (node-wise)

% Apply PCA
[coeff, score, latent, ~, explained] = pca(A_std);

n = length(A);

figure;
gscatter(score(:,1), score(:,2), 1:n)
text(score(:,1)+0.02, score(:,2), string(1:n))
xlabel('PC 1')
ylabel('PC 2')
title('PCA of Graph Adjacency Matrix with Node Labels')
grid on

D = diag(sum(A,2));      % Degree matrix
L = D - A;               % Unnormalized Laplacian
L_std = zscore(L, 0, 2); % Standardize

[coeff2, score2, ~, ~, explained2] = pca(L_std);

figure;
gscatter(score2(:,1), score2(:,2), 1:n)
text(score2(:,1)+0.02, score2(:,2), string(1:n))
xlabel('PC 1')
ylabel('PC 2')
title('PCA of Graph Adjacency Matrix with Node Labels')
grid on