function [centrality, details] = local_eigenvector_centrality(A, X, plotall, Imax)
% LOCAL_EIGENVECTOR_CENTRALITY
% Computes local eigenvector centrality of an adjacency matrix A.
%
% Usage:
%   centrality = local_eigenvector_centrality(A)
%   centrality = local_eigenvector_centrality(A, X, plotall)
%   centrality = local_eigenvector_centrality(A, X, plotall, Imax)
%
% Inputs:
%   A       - Adjacency matrix (sparse or full)
%   X       - (optional) Node coordinates [n x 2] for plotting
%   plotall - (optional) If true, generate plots via helper function
%   Imax    - (optional) Maximum number of eigenvectors to include
%
% Outputs:
%   centrality - Local eigenvector centrality values
%   details    - Struct containing eigenvalues, eigenvectors, eigengap,
%                and global eigenvector centrality for analysis/plotting
%
% Reference:
%   [Author(s)], "A Local Eigenvector Centrality", [Journal, Year].

    if nargin < 2, X = []; end
    if nargin < 3, plotall = false; end

    % === Eigen decomposition ===
    if issparse(A)
        [V, D] = eigs(A, 50, 'largestreal');
    else
        [V, D] = eig(A);
    end

    % Sort eigenvalues and eigenvectors
    [Dsort, idx] = sort(real(diag(D)), 'descend');
    V = V(:, idx);

    % Remove zero eigenvalue components
    zero_eig = (Dsort == 0);
    V(:, zero_eig) = 0;

    % Compute eigengap
    eigengap = Dsort(1:end-1) - Dsort(2:end);

    % === Determine Imax ===
    if nargin < 4 || isempty(Imax)
        [~, Imax] = max(eigengap);
        fprintf('Imax not provided — using largest eigengap at index %d.\n', Imax);
    else
        Imax = min(Imax, length(Dsort));
        fprintf('Using user-provided Imax = %d.\n', Imax);
    end

    % === Compute Local Centrality ===
    secondIMaxPeak = 0;
    centrality = vecnorm(V(:, secondIMaxPeak + 1:Imax), 2, 2);
    centrality = real(centrality);
    centrality(centrality <= 0) = eps;

    % === Global Centrality for Comparison ===
    global_centrality = abs(real(V(:, 1)));
    global_centrality(global_centrality <= 0) = eps;

    % === Package Results ===
    details = struct( ...
        'V', V, ...
        'D', Dsort, ...
        'eigengap', eigengap, ...
        'Imax', Imax, ...
        'global_centrality', global_centrality, ...
        'X', X, ...
        'A', A);

    % === Optional Visualization ===
    if plotall
        plot_local_eigenvector_centrality(details, centrality);
    end
end

function plot_local_eigenvector_centrality(details, local_centrality)
% PLOT_LOCAL_EIGENVECTOR_CENTRALITY
% Generates diagnostic plots for local vs global eigenvector centrality.
%
% Inputs:
%   details         - Struct from local_eigenvector_centrality()
%   local_centrality - Vector of local eigenvector centrality values

    A = details.A;
    X = details.X;
    V = details.V;
    D = details.D;
    eigengap = details.eigengap;
    Imax = details.Imax;
    global_centrality = details.global_centrality;

    G = digraph(A);

    % === Eigenvalue Spectrum ===
    figure; plot(D, '-o');
    ylabel('Eigenvalue'); xlabel('Index');
    title('Sorted Eigenvalues'); box off;

    figure; plot(eigengap, 'k-o');
    xlabel('i'); ylabel('\lambda_i - \lambda_{i+1}');
    title('Eigengap Spectrum'); box off;

    % === Plot eigenvectors up to Imax ===
    cmax = max(max(V(:, 1:Imax)));
    cmin = min(min(V(:, 1:Imax)));
    cabs = max([cmax, abs(cmin)]);

    for i = 1:Imax
        figure;
        ev = real(V(:, i));
        ev(ev == 0) = eps;

        if isempty(X)
            p = plot(G, 'Layout', 'force', ...
                'MarkerSize', (abs(ev) ./ sum(abs(ev))) * 100, ...
                'NodeCData', ev, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]);
        else
            p = plot(G, 'XData', X(:, 1), 'YData', X(:, 2), ...
                'MarkerSize', (abs(ev) ./ sum(abs(ev))) * 100, ...
                'NodeCData', ev, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]);
        end

        colormap(flip(summer));
        colorbar; caxis([-cabs cabs]);
        title(sprintf('Eigenvector %d', i)); axis off;
    end

    % === Compare Global vs Local Centrality ===
    values = {global_centrality, local_centrality};
    labels = {'Global Eigenvector Centrality', 'Local Eigenvector Centrality'};

    for i = 1:numel(values)
        v = values{i} ./ norm(values{i});


        figure;
        if isempty(X)
            plot(G, 'Layout', 'force', ...
                'MarkerSize', (abs(v) ./ sum(abs(v))) * 100, ...
                'NodeCData', v, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]);
        else
            plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
                'MarkerSize', (abs(v) ./ sum(abs(v))) * 100, ...
                'NodeCData', v, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]);
        end
        colormap(autumn); colorbar; title(labels{i});
        axis off;
    end

    % === Difference Map ===
    v1 = global_centrality ./ sum(global_centrality);
    vc = local_centrality ./ sum(local_centrality);
    pos = max(vc - v1, 1e-6);
    neg = max(v1 - vc, 1e-6);

    figure;
    if isempty(X)
        plot(G, 'Layout', 'force', 'MarkerSize', (abs(pos) ./ sum(abs(pos))) * 100, ...
             'NodeColor', [1 0 0], 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]); hold on;
        plot(G, 'Layout', 'force', 'MarkerSize', (abs(neg) ./ sum(abs(neg))) * 100, ...
             'NodeColor', [0 0 1], 'EdgeAlpha', 0, 'EdgeColor', [0,0,0]);
    else
        plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
             'MarkerSize', (abs(pos) ./ sum(abs(pos))) * 100, 'NodeColor', [1 0 0], 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0]); hold on;
        plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
             'MarkerSize', (abs(neg) ./ sum(abs(neg))) * 100, 'NodeColor', [0 0 1], 'EdgeAlpha', 0, 'EdgeColor', [0,0,0]);
    end

    legend('Local > Global', 'Global > Local');
    title('Difference Between Global and Local Centrality');
    axis off;
end
