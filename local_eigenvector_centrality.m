function [centrality] = local_eigenvector_centrality(A, X, plotall, Imax)
% local_eigenvector_centrality
% Computes and plots local eigenvector centrality of adjacency matrix A.
%
% Usage:
%   centrality = local_eigenvector_centrality(A)
%   centrality = local_eigenvector_centrality(A, X, plotall)
%   centrality = local_eigenvector_centrality(A, X, plotall, Imax)
%
% If Imax is not provided, it is determined automatically from the eigengap.

    % Check if A is sparse
    if issparse(A)
        [V, D] = eigs(A, 50, 'largestreal');  % Compute up to 34 eigenvalues/eigenvectors
    else
        [V, D] = eig(A);
    end
    
    % Extract and sort eigenvalues in descending order
    [Dsort, idx] = sort(real(diag(D)), 'descend'); 
    V = V(:, idx);

    % Handle possible complex conjugate eigenvectors
    Diff_Dsort = diff(real(Dsort));
    CC_I = find(Diff_Dsort == 0);
    for i = 1:length(CC_I)
        if imag(D(idx(CC_I(i)), idx(CC_I(i)))) > 0
            V(:, CC_I(i))   = real(V(:, CC_I(i)));
            V(:, CC_I(i)+1) = imag(V(:, CC_I(i)+1));
        end
    end

    % Remove zero eigenvalue components
    zero_eig = find(Dsort == 0);
    for i = 1:length(zero_eig)
        V(:, zero_eig(i)) = zeros(1, length(V(:,1)));
    end

    % Compute eigengap
    eigengap = Dsort(1:end-1) - Dsort(2:end);
    
    % === Determine Imax ===
    if nargin < 4 || isempty(Imax)
        [~, Imax] = max(eigengap);
        fprintf('Imax not provided — using largest eigengap at index %d.\n', Imax);
    else
        Imax = min(Imax, length(Dsort)); % safety check
        fprintf('Using user-provided Imax = %d.\n', Imax);
    end

    if plotall == 1
        % === Plot eigenvalues ===
        figure;
        plot(Dsort(:), '-o');
        ylabel('Eigenvalues');
        xlabel('Index');
        title('Sorted Eigenvalues');
    
        % === Plot eigengap ===
        figure;
        plot(eigengap, 'k-o');
        xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 14);
        ylabel('${\lambda}_i - {\lambda}_{i+1}$', 'Interpreter', 'latex', 'FontSize', 14);
        title('Eigengap Spectrum');
        box off;
    end

    % === Compute centrality ===
    secondIMaxPeak = 0;
    centrality = vecnorm((V(:, secondIMaxPeak + 1:Imax)), 2, 2);

    % === Create graph ===
    G = digraph(A);

    % For color scaling
    cmax = max(max(V(:, secondIMaxPeak + 1:Imax)));
    cmin = min(min(V(:, secondIMaxPeak + 1:Imax)));
    cabs = max([cmax, abs(cmin)]);

    % === Plot each eigenvector if requested ===
    if plotall == 1
        for i = 1:Imax
            figure;
            an_eigenvector = real(V(:,i));
            an_eigenvector(an_eigenvector == 0) = 0.00001;
        
            % Layout handling
            if nargin < 2 || isempty(X)
                p = plot(G, 'Layout', 'force', ...
                    'MarkerSize', (abs(an_eigenvector) ./ sum(abs(an_eigenvector))) * 50, ...
                    'NodeCData', an_eigenvector, 'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            else
                p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
                    'MarkerSize', (abs(an_eigenvector) ./ sum(abs(an_eigenvector))) * 50, ...
                    'NodeCData', an_eigenvector, 'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            end

            colormap(flip(summer));
            cb = colorbar;
            caxis([min(-cabs) max(cabs)]);
            ylabel(cb, sprintf('v_{%d}', i), 'FontSize', 12, 'Rotation', 270);
            axis off;
        end
    end

    % === Global vs Local Eigenvector Centrality ===
    global_centrality = real(V(:,1));
    global_centrality(global_centrality <= 0) = 0.00001;
    centrality(centrality == 0) = 0.00001;
    labels_toplot = {'Eigenvector centrality','Local eigenvector centrality'};
    values = {global_centrality, centrality};

    if plotall == 1
        for i = 1:length(values)
            values{i} = values{i} ./ norm(values{i});
            figure;
    
            if nargin < 2 || isempty(X)
                p = plot(G, 'Layout', 'force', ...
                    'MarkerSize', (values{i} ./ sum(values{i})) * 250, ...
                    'NodeCData', values{i}, 'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            else
                p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
                    'MarkerSize', (values{i} ./ sum(values{i})) * 250, ...
                    'NodeCData', values{i}, 'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            end
            
            colormap(autumn);
            cb = colorbar;
            ylabel(cb, labels_toplot{i}, 'FontSize', 12, 'Rotation', 270);
            axis off;
        end

    % === Compare local vs global ===
        v1 = (global_centrality) ./ sum(global_centrality);
        vc = (centrality) ./ sum(centrality);

        neg = v1 - vc;
        pos = vc - v1;
        neg(neg <= 0) = 0.00001;
        pos(pos <= 0) = 0.00001;

        figure;
        if nargin < 2 || isempty(X)
            p = plot(G, 'Layout', 'force', ...
                'MarkerSize', pos * 500, 'NodeColor', [1 0 0], ...
                'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            hold on
            plot(G, 'Layout', 'force', ...
                'MarkerSize', neg * 500, 'NodeColor', [0 0 1], ...
                'EdgeAlpha', 0, 'EdgeColor', [0 0 0]);
        else
            p = plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
                'MarkerSize', pos * 500, 'NodeColor', [1 0 0], ...
                'EdgeAlpha', 0.1, 'EdgeColor', [0 0 0]);
            hold on
            plot(G, 'XData', X(:,1), 'YData', X(:,2), ...
                'MarkerSize', neg * 500, 'NodeColor', [0 0 1], ...
                'EdgeAlpha', 0, 'EdgeColor', [0 0 0]);
        end
        
        legend('positive','negative');
        title('Difference Between Global and Local Centrality');
        axis off;
    end
end
