function results = evaluate_centrality_model()

clc;
close all;

%% ================= PARAMETERS =================
nodesPerHub = 20;
numHubs     = 6;

pInternal   = 0.15;
nInterList  = linspace(1,5,5);   % sweep inter-hub density

nP = length(nInterList);

%% ================= STORAGE =================
corr_local   = zeros(nP,1);
corr_global  = zeros(nP,1);
corr_katz    = zeros(nP,1);
corr_pr      = zeros(nP,1);

spearman_local  = zeros(nP,1);
rmse_local      = zeros(nP,1);
rankerr_local   = zeros(nP,1);

symmetry_local  = zeros(nP,1);
symmetry_global = zeros(nP,1);
symmetry_katz   = zeros(nP,1);
symmetry_pr     = zeros(nP,1);

%% ================= MAIN SWEEP =================
for t = 1:nP

    nInter = nInterList(t);

    %% ---- Build network ----
    seed = 1;
    [A, Pos, hubIdx] = build_modular_network(nodesPerHub, numHubs, pInternal, nInter, seed);

    G = graph(A);

    totalNodes = size(A,1);

    %% Normalize all measures for comparison
    normalize_vec = @(x) x(:) ./ max(x);

    %% ---- KNOWN-COMMUNITY EIGENVECTOR CENTRALITY (ground truth) ----
    KCEC = zeros(totalNodes,1);

    for h = 1:numHubs
        idx = hubIdx{h};
        Ahub = A(idx,idx);
        KCEC(idx) = centrality(graph(Ahub),'eigenvector');
    end
    KCEC = normalize_vec(KCEC);

    %% ---- LOCAL EIGENVECTOR ----
    LEC = local_eigenvector_centrality(A, Pos, false, 6);
    LEC = normalize_vec(LEC);

    %% ---- GLOBAL METRICS ----
    eigen_c = normalize_vec(centrality(G,'eigenvector'));

    DG = digraph(A);
    pagerank_c = normalize_vec(centrality(DG,'pagerank','FollowProbability',0.85));

    lambda_max = max(abs(eigs(sparse(A),1)));
    alpha = 0.85 / lambda_max;
    katz_c = normalize_vec((speye(totalNodes) - alpha*A) \ ones(totalNodes,1));

    [pcc] = principal_component_centrality(A);
    pcc = normalize_vec(pcc);

    if t==nP
        % polar_method_comparison(methods, methodNames)
        polar_method_comparison(KCEC, "KCEC", 'b')
        polar_method_comparison(eigen_c, "EC", 'g')
        polar_method_comparison(LEC, "LEC", 'g')
        polar_method_comparison(pagerank_c, "PageRank", 'g')
        polar_method_comparison(pcc, "PCC", 'g')
        polar_method_comparison(katz_c, "Katz", 'g')
    end

    % if t==nP
    %     plot_centrality_overlap(A, LEC, 30, 'b');
    %     title('LEC (i=6)')
    %     axis equal;
    % 
    %     plot_centrality_overlap(A, KCEC, 30, 'g');
    %     title('KCEC')
    %     axis equal;
    % 
    %     plot_centrality_overlap(A, eigen_c, 30, 'b');
    %     title('EC')
    %     axis equal;
    % 
    %     plot_centrality_overlap(A, pagerank_c, 30, 'b');
    %     title('PageRank')
    %     axis equal;
    % 
    %     plot_centrality_overlap(A, katz_c, 30, 'b');
    %     title('Katz')
    %     axis equal;
    % 
    %     plot_centrality_overlap(A, pcc, 30, 'b');
    %     title('PCC')
    %     axis equal;
    % end

    % if t ==1
    %     plot_centrality_colourvary(A, LEC, Pos, 15);
    %     title('Local Eigenvector (i=6)')
    %     axis equal;
    %     plot_centrality_colourvary(A, pcc, Pos, 15);
    %     title('PCC')
    %     axis equal;
    % else
    % 
    if t==nP
        plot_centrality_colourvary(A, LEC, Pos, 15);
        title('Local Eigenvector (i=6)')
        axis equal;
        plot_centrality_colourvary(A, pcc, Pos, 15);
        title('PCC')
        axis equal;
    end

    %% ================= EVALUATION =================

    % Pearson correlations using corrcoef
    C = corrcoef(KCEC, LEC);
    corr_local(t) = C(1,2);
    
    C = corrcoef(KCEC, eigen_c);
    corr_global(t) = C(1,2);
    
    C = corrcoef(KCEC, katz_c);
    corr_katz(t) = C(1,2);
    
    C = corrcoef(KCEC, pagerank_c);
    corr_pr(t) = C(1,2);

    C = corrcoef(KCEC, pcc);
    corr_pcc(t) = C(1,2);
    
    %% Spearman correlation
    % corrcoef does not support Spearman directly,
    % so convert to ranks first
    
    [~,r1] = sort(KCEC);
    [~,r1] = sort(r1);
    
    [~,r2] = sort(LEC);
    [~,r2] = sort(r2);
    
    C = corrcoef(r1, r2);
    spearman_local(t) = C(1,2);
    
    %% RMSE
    rmse_local(t) = sqrt(mean((KCEC - LEC).^2));
    
    %% Rank error
    [~,trueR] = sort(KCEC,'descend');
    [~,predR] = sort(LEC,'descend');
    
    rankerr_local(t) = mean(abs(trueR - predR));
    
    %% ================= SYMMETRY CONSISTENCY =================
    % Compare equivalent nodes across hubs
    % Lower = better preservation of modular symmetry
    
    nodeStd_local  = zeros(nodesPerHub,1);
    nodeStd_global = zeros(nodesPerHub,1);
    nodeStd_katz   = zeros(nodesPerHub,1);
    nodeStd_pr     = zeros(nodesPerHub,1);
    nodeStd_pcc    = zeros(nodesPerHub,1);
    
    for k = 1:nodesPerHub
    
        vals_local  = zeros(numHubs,1);
        vals_global = zeros(numHubs,1);
        vals_katz   = zeros(numHubs,1);
        vals_pr     = zeros(numHubs,1);
        vals_pcc     = zeros(numHubs,1);
    
        for h = 1:numHubs
    
            idx = (h-1)*nodesPerHub + k;
    
            vals_local(h)  = LEC(idx);
            vals_global(h) = eigen_c(idx);
            vals_katz(h)   = katz_c(idx);
            vals_pr(h)     = pagerank_c(idx);
            vals_pcc(h)     = pcc(idx);
    
        end
    
        % Variation across equivalent nodes
        nodeStd_local(k)  = std(vals_local);
        nodeStd_global(k) = std(vals_global);
        nodeStd_katz(k)   = std(vals_katz);
        nodeStd_pr(k)     = std(vals_pr);
        nodeStd_pcc(k)     = std(vals_pcc);
    
    end
    
    % Mean symmetry breaking
    symmetry_local(t)  = mean(nodeStd_local);
    symmetry_global(t) = mean(nodeStd_global);
    symmetry_katz(t)   = mean(nodeStd_katz);
    symmetry_pr(t)     = mean(nodeStd_pr);
    symmetry_pcc(t)     = mean(nodeStd_pcc);
end

%% ================= SUMMARY TABLE =================

results = table(nInterList', ...
    corr_local, corr_global, corr_katz, corr_pr, corr_pcc', ...
    spearman_local, rmse_local, rankerr_local, ...
    symmetry_local, symmetry_global, ...
    symmetry_katz, symmetry_pr, symmetry_pcc', ...
    'VariableNames', ...
    {'nInter', ...
     'Corr_Local', 'Corr_Global', 'Corr_Katz', 'Corr_PageRank', 'Corr_PCC', ...
     'Spearman_Local', 'RMSE_Local', 'RankError_Local', ...
     'Symmetry_Local', 'Symmetry_Global', ...
     'Symmetry_Katz', 'Symmetry_PageRank', 'Symmetry_PCC'});

disp(results);

%% ================= PLOTS =================

figure;
plot(nInterList, corr_local,'-o','LineWidth',2); hold on;
plot(nInterList, corr_global,'-o','LineWidth',2);
plot(nInterList, corr_katz,'-o','LineWidth',2);
plot(nInterList, corr_pr,'-o','LineWidth',2);
plot(nInterList, corr_pcc,'-o','LineWidth',2);

xlabel('No. of inter-hub connections','FontSize',12);
ylabel('Correlation with KCEC','FontSize',12);
legend('Local Eigenvector','Global Eigenvector','Katz','PageRank','PCC');
% title('Similarity to KCEC');
grid on;

% %% Spearman
% figure;
% plot(nInterList, spearman_local,'-o','LineWidth',2);
% xlabel('Inter-hub probability');
% ylabel('Spearman correlation');
% title('Rank Agreement: LEC vs KCEC');
% grid on;
% 
% %% RMSE
% figure;
% plot(nInterList, rmse_local,'-o','LineWidth',2);
% xlabel('Inter-hub probability');
% ylabel('RMSE');
% title('Error: LEC vs KCEC');
% grid on;

%% ================= SYMMETRY PRESERVATION =================

figure;

plot(nInterList, symmetry_local,'-o','LineWidth',2); hold on;
plot(nInterList, symmetry_global,'-o','LineWidth',2);
plot(nInterList, symmetry_katz,'-o','LineWidth',2);
plot(nInterList, symmetry_pr,'-o','LineWidth',2);
plot(nInterList, symmetry_pcc,'-o','LineWidth',2);

xlabel('No. of inter-hub connections','FontSize',12);
ylabel('Node centrality variation','FontSize',12);

legend('Local Eigenvector', ...
       'Global Eigenvector', ...
       'Katz', ...
       'PageRank', ...
       'PCC');

% title('Symmetry Preservation Across Identical Hubs');

grid on;

end

function [A, Pos, hubIdx] = build_modular_network( ...
    nodesPerHub, numHubs, pInternal, nInter, seed)

    % Optional argument
    if nargin >= 5 && ~isempty(seed)
        rng(seed);
    end

    totalNodes = nodesPerHub * numHubs;

    %% --- base connected hub ---
    baseAdj = zeros(nodesPerHub);

    perm = randperm(nodesPerHub);
    for i = 2:nodesPerHub
        n1 = perm(i);
        n2 = perm(randi(i-1));
        baseAdj(n1,n2)=1;
        baseAdj(n2,n1)=1;
    end

    extra = rand(nodesPerHub) < pInternal;
    extra = triu(extra,1);
    extra = extra + extra.';
    baseAdj = baseAdj | extra;
    baseAdj(1:nodesPerHub+1:end)=0;

    %% --- replicate ---
    A = zeros(totalNodes);
    hubIdx = cell(numHubs,1);

    for h = 1:numHubs
        idx = (h-1)*nodesPerHub + (1:nodesPerHub);
        A(idx,idx)=baseAdj;
        hubIdx{h}=idx;
    end

    %% --- inter hub edges ---
    for k = 1:nInter
        for h1=1:numHubs-1
            for h2=h1+1:numHubs
                i1 = hubIdx{h1}(randi(nodesPerHub));
                i2 = hubIdx{h2}(randi(nodesPerHub));

                A(i1,i2)=1;
                A(i2,i1)=1;
            end
        end
    end

    % %% --- geometry ---
    % theta = linspace(0,2*pi,nodesPerHub+1)';
    % theta(end)=[];
    % 
    % hubR = 2;
    % 
    % hubX = hubR*cos(theta);
    % hubY = hubR*sin(theta);
    % 
    % bigTheta = linspace(0,2*pi,numHubs+1)';
    % bigTheta(end)=[];
    % 
    % bigR = 6;
    % 
    % cx = bigR*cos(bigTheta);
    % cy = bigR*sin(bigTheta);
    % 
    % Pos = zeros(totalNodes,2);
    % 
    % for h=1:numHubs
    %     idx = hubIdx{h};
    %     Pos(idx,1)=hubX+cx(h);
    %     Pos(idx,2)=hubY+cy(h);
    % end

    %% --- geometry ---
    theta = linspace(0,2*pi,nodesPerHub+1)';
    theta(end) = [];
    
    hubR = 2;
    hubX = hubR*cos(theta);
    hubY = hubR*sin(theta);
    
    % Hub centre positions: 2 rows × 3 columns
    dx = 5;   % horizontal spacing between hub centres
    dy = 2.5;   % vertical spacing between hub centres
    
    cx = [-dx  0  dx  -dx  0  dx]';
    cy = [ dy  dy dy  -dy -dy -dy]';
    
    Pos = zeros(totalNodes,2);
    
    for h = 1:numHubs
        idx = hubIdx{h};
        Pos(idx,1) = hubX + cx(h);
        Pos(idx,2) = hubY + cy(h);
    end

end

function polar_method_comparison(method, methodName, colour)

    nNodes = 20;
    nComm  = 6;

    theta = linspace(0,2*pi,nNodes+1);
    theta(end) = [];

    figure
    hold on

    % Draw node spokes
    for k = 1:nNodes

        plot([0 cos(theta(k))], ...
             [0 sin(theta(k))], ...
             ':', ...
             'Color',[0.8 0.8 0.8]);

        text(1.1*cos(theta(k)), ...
             1.1*sin(theta(k)), ...
             num2str(k), ...
             'HorizontalAlignment','center');
    end

    % Draw concentric circles
    t = linspace(0,2*pi,200);

    for r = [0.25 0.5 0.75 1.0]
        plot(r*cos(t),r*sin(t),...
            ':','Color',[0.8 0.8 0.8]);
    end

    % colours = lines(nComm);

    for c = 1:nComm

        idx = (c-1)*nNodes + (1:nNodes);

        r = method(idx);
        r = r/max(r);

        r(end+1) = r(1);

        theta2 = [theta theta(1)];

        x = r'.*cos(theta2);
        y = r'.*sin(theta2);

        plot(x,y,...
             'LineWidth',2,...
             'Color',colour);
    end

    axis equal
    axis off

    title(methodName,'Units','normalized')
    h = get(gca,'Title');
    h.Position(2) = 1.1; % default is roughly 1.0

    % legend(compose('Community %d',1:nComm),...
    %        'Location','eastoutside');
end

function plot_centrality_colourvary(A, centrality, X, scaled)
% Plot graph with marker sizes scaled from centrality. 'scaled' controls maximum marker multiplication.

    G = graph(A);

    % Marker sizing: consistent scaling across networks
    ms = scale_markers(centrality, 6, 35);
    ms = ms / max(ms) * scaled;

    figure;
    if isempty(X)
        p = plot(G, 'Layout', 'force', 'MarkerSize', ms, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0],'HandleVisibility','off');
    else
        p = plot(G, 'XData', X(:,1), 'YData', X(:,2), 'MarkerSize', ms, 'EdgeAlpha', 0.1, 'EdgeColor', [0,0,0],'HandleVisibility','off');
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

function [pcc] = principal_component_centrality(A)
    
    %% ---------------------------
    % 2. Eigen decomposition
    % ---------------------------
    [X, D] = eig(A);
    lambda = diag(D);
    
    % Sort eigenvalues (descending magnitude)
    [~, idx] = sort(abs(lambda), 'descend');
    
    lambda = lambda(idx);
    X = X(:, idx);

    N = length(lambda);
    
    %% ---------------------------
    % 2a. Compute Eigenvector Centrality (baseline)
    % ---------------------------
    v = X(:,1);                 % principal eigenvector
    C_evc = abs(v) / norm(v);
    
    %% ---------------------------
    % 2b. Phase-angle stability curve
    % ---------------------------
    P_max = min(N, 100);        % limit for efficiency (adjust as needed)
    phi = zeros(P_max,1);
    
    for p = 1:P_max
        
        % Select top p eigenvectors
        Xp_temp = X(:,1:p);
        lambda_temp = lambda(1:p);
        
        % Compute PCC for this p
        Z = Xp_temp .* lambda_temp';
        C_p = sqrt(sum(Z.^2, 2));
        
        % Compute phase angle with EVC
        phi(p) = acos( dot(C_p, C_evc) / (norm(C_p)*norm(C_evc)) );
    end
    
    %% ---------------------------
    % 2c. Select P based on plateau detection
    % ---------------------------
    % Compute first derivative (change in phase angle)
    dphi = abs(diff(phi));
    
    % Threshold for "stability"
    tol = 1e-3;
    
    % Find first point where changes become small
    P_phase = find(dphi < tol, 1, 'first') + 1;
    
    % Safety fallback
    if isempty(P_phase)
        P_phase = round(P_max/2);
    end
    
    fprintf('Phase-based P = %d\n', P_phase);
    
    %% ---------------------------
    % 2e. Choose final P
    % ---------------------------
    % Prefer phase-based selection (as in paper)
    P = P_phase;
    
    fprintf('Final selected P = %d\n', P);
    
    %% ---------------------------
    % 2b. Extract top P components
    % ---------------------------
    Xp = X(:, 1:P);
    lambda_p = lambda(1:P);

    
    %% ---------------------------
    % 3. Compute PCC
    % ---------------------------
    % PCC definition:
    % C(i) = sqrt( sum_{k=1}^{P} (lambda_k * x_k(i))^2 )
    % Equivalent to L2 norm in eigenspace

    Z = Xp .* lambda_p';     % scale eigenvectors
    pcc = sqrt(sum(Z.^2, 2));

end