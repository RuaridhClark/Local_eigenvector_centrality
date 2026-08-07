%% ================= PARAMETERS =================
clc;
close all;

N = 500;
pInternal   = 0.15;

seedList    = 1:20;
nSeeds      = length(seedList);

n_percomm = 100;   % nodes per community
n_comms = 2;     % communities
gamma = 1.5;   % power law exponent
d_min = 1;    % min degree


%% ===== LOCALISATION TEST =====

eliteWeight = 100;    % very strong weight

%% ================= STORAGE =================

corr_LEC  = zeros(1,nSeeds);
corr_mitigate = zeros(1,nSeeds);
corr_pr = zeros(1,nSeeds);

%% ================= MAIN SWEEP =================

for s = 1:nSeeds

    seed = seedList(s);

    %% Erdos-Renyi topology + heavy-tailed weights

    rng(seed)

    N = 200;
    p = 0.03;
    
    AT = rand(N) < p;
    AT = triu(AT,1);
    AT = AT + AT';
    edges = find(triu(AT,1));
    
    % Fixed normal samples
    Z = randn(length(edges),1);

    sigma = 0.01;
    w = exp(sigma*Z - sigma^2/2);
    A0 = zeros(N);
    A0(edges) = w;
    A0 = A0 + A0';
    
    % N = 120;       % Number of nodes
    % gamma = 15;    % Power-law exponent
    
    % [G0,A0] = generate_uncorrelated_scale_free(N, gamma);

    % [A0,edges]=build_WRR_base(...
    %     k,...
    %     n1,...
    %     n2,...
    %     ncandidates,...
    %     1);

    % [A0,Pos,C,k]=build_scale_free_communities(...
    %     n_percomm, n_comms, gamma, d_min, 1, 1);

    totalNodes0 = size(A0,1);
    
    % normalize_vec = @(x) x(:)./max(abs(x));
    normalize_vec = @(x) x(:)./norm(x);

    G0 = graph(A0);


    %% ================= ORIGINAL NETWORK =================

    %% -------- LEC --------

    LEC0 = local_eigenvector_centrality(A0,[],false);
    LEC0 = normalize_vec(LEC0);


    %% ================= LOCALISATION ATTACK =================

    % [A,Pos,C,k]=build_scale_free_communities(...
    %     n_percomm, n_comms, gamma, d_min, 30, 1);
   
    % [G,A] = reduce_localization(G0);

    sigma = 1;
    w = exp(sigma*Z - sigma^2/2);
    A = zeros(N);
    A(edges) = w;
    A = A + A';
   
    totalNodes = size(A,1);

    % interDeg = count_additional_edges(A,nodesPerHub,numHubs,hubIdx);

    G = graph(A);

    %% ================= AUGMENTED NETWORK =================

    %% -------- LEC --------
    M = ones(length(A),length(A));

    LEC = local_eigenvector_centrality(A,[],false);
    LEC = normalize_vec(LEC);

    %% -------- OTHER --------

    DG = digraph(A);
    pr = normalize_vec( ...
        centrality(DG,'pagerank', ...
        'FollowProbability',0.85));

    %% --- Optimize power for LEC_adjust centrality ---
    pow_list = 0.05:0.05:1;
    wsd_power_values = zeros(size(pow_list));

    LECtemp = local_eigenvector_centrality(A+M,[],false);
    
    for i = 1:length(pow_list)
        pow = pow_list(i);
        LEC_adjust = LECtemp.^pow;
        LEC_adjust = normalize_vec(LEC_adjust);
        difference_vector = pr - LEC_adjust;
        wsd_power_values(i) = norm(difference_vector, 2);
    end
    
    [~, idx_min] = min(wsd_power_values);
    pow_opt = pow_list(idx_min);
    
    LEC_adjust = LECtemp.^pow_opt ./ sum(LECtemp.^pow_opt);

    LEC_adjust = normalize_vec(LEC_adjust);        

    %% -------- CORRELATION --------

    C = corrcoef(LEC0, LEC);
    corr_LEC(1,s) = C(1,2);
    C = corrcoef(LEC0, LEC_adjust);
    corr_mitigate(1,s) = C(1,2);
    C = corrcoef(LEC0, pr);
    corr_pr(1,s) = C(1,2);

    % Plot LEC0, LEC, LEC_adjust, pr
    figure;
    hold on;
    plot(LEC0, 'DisplayName', 'Original LEC', 'LineWidth', 1.5);
    plot(LEC, 'DisplayName', 'LEC', 'LineWidth', 1.5);
    plot(LEC_adjust, 'DisplayName', 'LEC Adjusted', 'LineWidth', 1.5);
    plot(pr, 'DisplayName', 'PageRank', 'LineWidth', 1.5);
    legend('show');
    title(sprintf('LEC Comparison for Seed = %d', seed));
    xlabel('Node Index');
    ylabel('LEC Value');
    grid on;

    %%%%%%%%%
    c_vals = [LEC0, LEC, LEC_adjust, pr];
    for i = 1:4
        c = c_vals(:,i);
        ms = 5 + 20*(c-min(c))/(max(c)-min(c));

        figure
        p = plot(G,'Layout','force');

        p.MarkerSize = ms;
        p.NodeCData = c;

        ew = G.Edges.Weight;
        p.LineWidth = 0.5 + 2*rescale(log(ew));

        colormap(autumn)
        cb = colorbar;
        cb.Label.String = 'Eigenvector centrality';

        title('Weighted Erdos-Rényi network')
        axis off
    end


end

%% ================= SUMMARY =================

corr_LEC_mu  = mean(corr_LEC,2);
corr_mitigate_mu = mean(corr_mitigate,2);

corr_LEC_sd  = std(corr_LEC,0,2);
corr_mitigate_sd = std(corr_mitigate,0,2);

% Boxplot of corr_LEC vs corr_mitigate
% Prepare data for boxchart comparing corr_LEC and corr_mitigate across nInterList

% Aggregated values across all nInter 
y_all_LEC = corr_LEC(:);
y_all_mitigate = corr_mitigate(:);
y_all_pr = corr_pr(:);

% Group labels
g = [repmat("LEC", numel(y_all_LEC), 1); repmat("LEC (mitigate)", numel(y_all_mitigate), 1); repmat("PageRank", numel(y_all_pr), 1)];
y = [y_all_LEC; y_all_mitigate; y_all_pr];
g = categorical(g);

% Statistical test (two-sample t-test) on pooled distributions
[~, p] = ttest2(y_all_LEC, y_all_mitigate);
fprintf('Comparison LEC vs LEC (mitigate): t-test p-value = %.6f\n', p);

% Boxchart
figure;
boxchart(g, y);
ylabel('Correlation with original LEC');
title('LEC vs LEC\_mitigate (pooled across seeds)');
grid on;


%% ================= HELPER =================
function [G,A] = generate_uncorrelated_scale_free(N, gamma, k_min)
    % 1. Set default minimum degree if not provided
    if nargin < 3
        k_min = 3;
    end
    
    % Structural cutoff to prevent degree correlations: k_max <= sqrt(N)
    k_max = floor(sqrt(N)); 
    
    % 2. Generate valid power-law degree sequence
    degrees = zeros(N, 1);
    count = 0;
    while count < N
        % Generate candidates using Inverse Transform Sampling for Power-Law
        u = rand(N - count, 1);
        k = floor(k_min * (1 - u).^(-1 / (gamma - 1)));
        
        % Filter candidates within the allowed UCM range
        valid_idx = (k >= k_min) & (k <= k_max);
        valid_k = k(valid_idx);
        
        % Store valid degrees
        num_valid = min(length(valid_k), N - count);
        degrees(count + 1 : count + num_valid) = valid_k(1:num_valid);
        count = count + num_valid;
    end
    
    % Enforce an even sum of degrees for the stub matching procedure
    if mod(sum(degrees), 2) ~= 0
        [~, max_idx] = max(degrees);
        degrees(max_idx) = degrees(max_idx) + 1;
    end
    
    % 3. Build network using the Configuration Model (Stub Matching)
    stubs = repeating_indices(degrees);
    stubs = stubs(randperm(length(stubs))); % Randomly shuffle the stubs
    
    % Pair adjacent stubs into edges
    source_nodes = stubs(1:2:end);
    target_nodes = stubs(2:2:end);
    
    % Create graph and clean up structural defects (self-loops & multiple edges)
    G = graph(source_nodes, target_nodes);
    G = simplify(G); 

    % Create adjacency matrix from G
    A = adjacency(G);

end

function stubs = repeating_indices(degrees)
    % Helper function to replicate node IDs based on their degree
    stubs = zeros(sum(degrees), 1);
    idx = 1;
    for node = 1:length(degrees)
        deg = degrees(node);
        stubs(idx : idx + deg - 1) = node;
        idx = idx + deg;
    end
end

function [G_delocalized,A] = reduce_localization(G, num_swaps)
    % 1. Extract the existing edge list from your graph
    edges = G.Edges.EndNodes;
    num_edges = size(edges, 1);
    
    % Default swaps: typically 10 to 20 times the number of edges for full mixing
    if nargin < 2
        num_swaps = 10 * num_edges;
    end
    
    successful_swaps = 0;
    
    % 2. Perform Markov Chain edge swapping
    for i = 1:num_swaps
        % Select two random distinct edges
        idx1 = randi(num_edges);
        idx2 = randi(num_edges);
        while idx1 == idx2
            idx2 = randi(num_edges);
        end
        
        A = edges(idx1, 1); B = edges(idx1, 2);
        C = edges(idx2, 1); D = edges(idx2, 2);
        
        % Avoid creating self-loops or multi-edges during the swap
        % Check if new connections (A-D) or (C-B) already exist or overlap
        if (A ~= D) && (C ~= B) && ...
           ~any((edges(:,1) == A & edges(:,2) == D) | (edges(:,1) == D & edges(:,2) == A)) && ...
           ~any((edges(:,1) == C & edges(:,2) == B) | (edges(:,1) == B & edges(:,2) == C))
            
            % Execute the swap
            edges(idx1, 2) = D;
            edges(idx2, 2) = B;
            successful_swaps = successful_swaps + 1;
        end
    end
    
    % 3. Rebuild the modified graph structure
    G_delocalized = graph(edges(:,1), edges(:,2));

    % Create adjacency matrix from G
    A = adjacency(G_delocalized);
    
    fprintf('Completed %d successful degree-preserving swaps.\n', successful_swaps);
end


function [A,Pos,community,degrees] = ...
    build_scale_free_communities(...
    nodesPerCommunity,...
    numCommunities,...
    gamma,...
    kmin,...
    interEdges,...
    seed)

% BUILD_SCALE_FREE_COMMUNITIES
%
% Creates a modular network where each community has
% scale-free degree statistics.
%
% Inputs:
%
% nodesPerCommunity : nodes in each community
% numCommunities    : number of communities
% gamma             : power-law exponent (2-3 typical)
% kmin              : minimum degree
% interEdges        : number of edges between communities
% seed              : random seed
%
% Outputs:
%
% A          adjacency matrix
% Pos        node positions
% community  community labels
% degrees    node degrees


rng(seed);


%% -------------------------------------------------------
% Initialise
% --------------------------------------------------------

N = nodesPerCommunity*numCommunities;

A=zeros(N);

community=zeros(N,1);


%% -------------------------------------------------------
% Generate communities
% --------------------------------------------------------

for c=1:numCommunities


    idx=(c-1)*nodesPerCommunity+(1:nodesPerCommunity);

    community(idx)=c;


    % create community graph
    m0 = 5;
    m  = 2;

    Ac=generate_scale_free_graph(...
        nodesPerCommunity,...
        m0,...
        m,...
        seed+c);


    A(idx,idx)=Ac;


end



%% -------------------------------------------------------
% Add inter-community connections
% --------------------------------------------------------

for e=1:interEdges


    c1=randi(numCommunities);

    c2=randi(numCommunities);


    while c2==c1
        c2=randi(numCommunities);
    end


    n1=find(community==c1);
    n2=find(community==c2);


    i=n1(randi(length(n1)));
    j=n2(randi(length(n2)));


    A(i,j)=1;
    A(j,i)=1;


end



%% -------------------------------------------------------
% Degrees
% --------------------------------------------------------

degrees=sum(A,2);



%% -------------------------------------------------------
% Simple geometry
% --------------------------------------------------------

Pos=zeros(N,2);


for c=1:numCommunities

    idx=find(community==c);


    theta=2*pi*rand;

    centre=[5*cos(theta),5*sin(theta)];


    Pos(idx,:)=centre+0.5*randn(length(idx),2);

end


end

function A = generate_scale_free_graph(N,m0,m,seed)

rng(seed);


if m >= m0
    error('m must be smaller than m0')
end


A=zeros(N);


%% Initial fully connected core

for i=1:m0

    for j=i+1:m0

        A(i,j)=1;
        A(j,i)=1;

    end

end


degree=sum(A,2);



%% Preferential attachment

for new=m0+1:N


    targets=[];


    while length(targets)<m

        p=degree(1:new-1)/sum(degree(1:new-1));

        candidate=randsample(new-1,1,true,p);


        if ~ismember(candidate,targets)

            targets=[targets candidate];

        end

    end



    for t=targets

        A(new,t)=1;
        A(t,new)=1;

    end


    degree=sum(A,2);

end


end


% function A = apply_localisation_attack(A,nEliteNodes,w)
% 
%     N = size(A,1);
% 
%     elite = randperm(N,nEliteNodes);
% 
%     for i = 1:nEliteNodes-1
% 
%         for j = i+1:nEliteNodes
% 
%             A(elite(i),elite(j)) = w;
%             A(elite(j),elite(i)) = w;
% 
%         end
% 
%     end
% 
% end

function [A,eliteNodes] = apply_localisation_attack(A,hubIdx,LEC0,w)

    numHubs = numel(hubIdx);
    
    eliteNodes = zeros(numHubs,1);
    
    %% Select highest-LEC node in each hub
    for h = 1:numHubs
    
        idx = hubIdx{h};
    
        [~,m] = max(LEC0(idx));
    
        eliteNodes(h) = idx(m);
    
    end
    
    % %% Create weighted clique between elite nodes
    % for i = 1:numHubs
    % 
    %     nbrs = find(A(eliteNodes(i),:));
    % 
    %     A(eliteNodes(i),nbrs) = w;
    %     A(nbrs,eliteNodes(i)) = w;
    % 
    % end

    for i = 1:numHubs-1
    
        for j = i+1:numHubs
    
            A(eliteNodes(i),eliteNodes(j)) = ...
                w;

            A(eliteNodes(j),eliteNodes(i)) = ...
                w;
    
        end
    
    end

end

function shaded_line(x,mu,sigma,col)

fill([x fliplr(x)], ...
     [(mu-sigma)' fliplr((mu+sigma)')], ...
     col,...
     'FaceAlpha',0.2,...
     'EdgeColor','none');

hold on

plot(x,mu,...
    'Color',col,...
    'LineWidth',2);

end

function A = generate_random_regular(N,k)

success=false;


while ~success

    stubs=repmat(1:N,1,k);

    stubs=stubs(randperm(length(stubs)));


    A=zeros(N);


    success=true;


    for i=1:2:length(stubs)

        a=stubs(i);
        b=stubs(i+1);


        if a==b || A(a,b)==1

            success=false;
            break

        end


        A(a,b)=1;
        A(b,a)=1;

    end


end

end

function [A0,candidateEdges] = ...
    build_WRR_base(k,n1,n2,nCandidates,seed)

    rng(seed)
    
    % create wheel
    W=zeros(n1);
    
    for i=2:n1
        W(1,i)=1;
        W(i,1)=1;
    end
    
    for i=2:n1-1
        W(i,i+1)=1;
        W(i+1,i)=1;
    end
    
    W(2,n1)=1;
    W(n1,2)=1;
    
    
    % random regular
    R=generate_random_regular(n2,k);
    
    
    N=n1+n2;
    
    A0=zeros(N);
    
    A0(1:n1,1:n1)=W;
    
    A0(n1+1:end,n1+1:end)=R;
    
    
    %% possible extra links
    
    candidateEdges=zeros(nCandidates,2);
    
    
    for e=1:nCandidates
        
        i=randi(n1);
        j=n1+randi(n2);
    
        candidateEdges(e,:)=[i j];
    
    end


end

function [A,Pos,pev,IPR,lambda1,lambda2,deg] = ...
    build_Wheel_RandomRegular(k,n1,n2,seed)

    % BUILD_WHEEL_RANDOMREGULAR
    %
    % Creates the wheel-random-regular graph used for eigenvector localisation
    % studies.
    %
    % Inputs:
    %
    % k     average degree of random regular graph
    % n1    wheel graph size
    % n2    random regular graph size
    % seed  random seed
    %
    % Outputs:
    %
    % A       adjacency matrix
    % Pos     node positions
    % pev     principal eigenvector
    % IPR     inverse participation ratio
    % lambda1 largest eigenvalue
    % lambda2 second eigenvalue
    % deg     degree sequence
    
    
    rng(seed);
    
    
    %% =========================================================
    % 1. Create wheel graph
    % ==========================================================
    
    W = zeros(n1);
    
    
    % hub node = 1
    
    for i=2:n1
        
        W(1,i)=1;
        W(i,1)=1;
        
    end
    
    
    % outer ring
    
    for i=2:n1-1
        
        W(i,i+1)=1;
        W(i+1,i)=1;
        
    end
    
    
    W(2,n1)=1;
    W(n1,2)=1;
    
    
    
    %% =========================================================
    % 2. Create random regular graph
    % ==========================================================
    
    % k*n2 must be even
    
    if mod(k*n2,2)~=0
        error('k*n2 must be even')
    end
    
    
    R = generate_random_regular(n2,k);
    
    
    
    %% =========================================================
    % 3. Combine networks
    % ==========================================================
    
    N = n1+n2+1;
    
    A=zeros(N);
    
    
    % wheel
    
    A(1:n1,1:n1)=W;
    
    
    % bridge node
    
    bridge=n1+1;
    
    A(n1,bridge)=1;
    A(bridge,n1)=1;

    % Random nodes to study mitigation of localisation
    for e=1:m
    
        i=randi(n1);
        j=n1+1+randi(n2);
        
        A(i,j)=1;
        A(j,i)=1;
        
    end
    
    
    % random regular block
    
    idx=(n1+2):N;
    
    A(idx,idx)=R;
    
    
    % bridge into random regular graph
    
    A(bridge,n1+2)=1;
    A(n1+2,bridge)=1;
    
    
    
    %% =========================================================
    % 4. Eigenvector localisation
    % ==========================================================
    
    [V,D]=eig(A);
    
    
    [eigenvalues,order]=sort(diag(D));
    
    
    lambda1=eigenvalues(end);
    lambda2=eigenvalues(end-1);
    
    
    pev=abs(V(:,order(end)));
    
    
    IPR=sum(pev.^4);
    
    
    
    %% =========================================================
    % 5. Degree sequence
    % =========================================================
    
    deg=sum(A,2);
    
    
    
    %% =========================================================
    % 6. Geometry
    % ==========================================================
    
    Pos=zeros(N,2);
    
    
    % wheel coordinates
    
    theta=linspace(0,2*pi,n1-1);
    
    Pos(1,:)=[0 0];
    
    Pos(2:n1,1)=cos(theta);
    Pos(2:n1,2)=sin(theta);
    
    
    
    % random regular graph placement
    
    theta2=2*pi*rand(n2,1);
    
    r=5+rand(n2,1);
    
    Pos(n1+2:end,1)=r.*cos(theta2);
    Pos(n1+2:end,2)=r.*sin(theta2);
    
    
    
    % bridge
    
    Pos(bridge,:)=[2 0];


end

function [A,Pos,core,degree] = build_richclub_localisation(...
    N,gamma,kmin,coreFrac,rho,alpha,seed)

    % BUILD_RICHCLUB_LOCALISATION
    %
    % Synthetic network with tunable eigenvector localisation.
    %
    % Inputs:
    %
    % N          number of nodes
    % gamma      power-law exponent
    % kmin       minimum degree
    % coreFrac   fraction of nodes in rich club
    % rho        rich-club density (0-1)
    % alpha      rewiring fraction (0-1)
    % seed       random seed
    %
    % Outputs:
    %
    % A          adjacency matrix
    % Pos        node coordinates
    % core       indices of rich-club nodes
    % degree     degree sequence
    
    
    rng(seed);
    
    
    %% =========================================================
    % 1. Power-law degree sequence
    % ==========================================================
    
    u = rand(N,1);
    
    degree = floor(kmin*(1-u).^(-1/(gamma-1)));
    
    degree = max(degree,2);
    
    % avoid unrealistic degrees
    degree(degree>N/4)=N/4;
    
    
    %% =========================================================
    % 2. Chung-Lu background network
    % ==========================================================
    
    A = zeros(N);
    
    S = sum(degree);
    
    
    for i=1:N-1
    
        for j=i+1:N
    
            p = degree(i)*degree(j)/S;
    
            p = min(p,0.95);
    
            if rand<p
                A(i,j)=1;
                A(j,i)=1;
            end
    
        end
    
    end
    
    
    
    %% =========================================================
    % 3. Create rich club
    % ==========================================================
    
    [~,order] = sort(sum(A,2),'descend');
    
    
    ncore = round(coreFrac*N);
    
    core = order(1:ncore);
    
    
    % Add edges inside core
    
    for i=1:ncore-1
    
        for j=i+1:ncore
    
            if rand < rho
    
                A(core(i),core(j))=1;
                A(core(j),core(i))=1;
    
            end
    
        end
    
    end
    
    
    
    %% =========================================================
    % 4. Degree-preserving rewiring outside rich club only
    % ==========================================================
    
    % Identify allowed edges:
    % at least one endpoint must be outside the core
    
    coreMask = false(N,1);
    coreMask(core)=true;
    
    
    [i,j] = find(triu(A));
    
    
    allowed = ~(coreMask(i) & coreMask(j));
    
    
    i = i(allowed);
    j = j(allowed);
    
    
    numSwap = round(alpha*length(i));
    
    
    for s=1:numSwap
    
    
        % refresh candidate edge list
        [i,j] = find(triu(A));
    
    
        allowed = ~(coreMask(i) & coreMask(j));
    
        i=i(allowed);
        j=j(allowed);
    
    
        if length(i)<2
            break
        end
    
    
        e = randperm(length(i),2);
    
    
        a=i(e(1));
        b=j(e(1));
    
        c=i(e(2));
        d=j(e(2));
    
    
        % prevent invalid swaps
    
        if length(unique([a b c d]))<4
            continue
        end
    
    
        if A(a,d) || A(c,b)
            continue
        end
    
    
        % perform swap
    
        A(a,b)=0;
        A(b,a)=0;
    
        A(c,d)=0;
        A(d,c)=0;
    
    
        A(a,d)=1;
        A(d,a)=1;
    
        A(c,b)=1;
        A(b,c)=1;
    
    end
    
    
    %% =========================================================
    % 5. Geometry
    % ==========================================================
    
    Pos=zeros(N,2);
    
    
    % place core centrally
    
    Pos(core,:) = 0.5*randn(length(core),2);
    
    
    % place remaining nodes around core
    
    outside=setdiff(1:N,core);
    
    
    theta=2*pi*rand(length(outside),1);
    
    r=5+5*rand(length(outside),1);
    
    
    Pos(outside,1)=r.*cos(theta);
    Pos(outside,2)=r.*sin(theta);

end

% function [A, Pos, hubIdx] = build_modular_network( ...
%     nodesPerHub, numHubs, pInternal, nInter, seed)
% 
%     % Optional argument
%     if nargin >= 5 && ~isempty(seed)
%         rng(seed);
%     end
% 
%     totalNodes = nodesPerHub * numHubs;
% 
%     %% --- base connected hub ---
%     baseAdj = zeros(nodesPerHub);
% 
%     perm = randperm(nodesPerHub);
%     for i = 2:nodesPerHub
%         n1 = perm(i);
%         n2 = perm(randi(i-1));
%         baseAdj(n1,n2)=1;
%         baseAdj(n2,n1)=1;
%     end
% 
%     extra = rand(nodesPerHub) < pInternal;
%     extra = triu(extra,1);
%     extra = extra + extra.';
%     baseAdj = baseAdj | extra;
%     baseAdj(1:nodesPerHub+1:end)=0;
% 
%     %% --- replicate ---
%     A = zeros(totalNodes);
%     hubIdx = cell(numHubs,1);
% 
%     for h = 1:numHubs
%         idx = (h-1)*nodesPerHub + (1:nodesPerHub);
%         A(idx,idx)=baseAdj;
%         hubIdx{h}=idx;
%     end
% 
%     %% --- inter hub edges ---
%     for k = 1:nInter
%         for h1=1:numHubs-1
%             for h2=h1+1:numHubs
%                 i1 = hubIdx{h1}(randi(nodesPerHub));
%                 i2 = hubIdx{h2}(randi(nodesPerHub));
% 
%                 A(i1,i2)=1;
%                 A(i2,i1)=1;
%             end
%         end
%     end
% 
%     %% --- geometry ---
%     theta = linspace(0,2*pi,nodesPerHub+1)';
%     theta(end)=[];
% 
%     hubR = 2;
% 
%     hubX = hubR*cos(theta);
%     hubY = hubR*sin(theta);
% 
%     bigTheta = linspace(0,2*pi,numHubs+1)';
%     bigTheta(end)=[];
% 
%     bigR = 6;
% 
%     cx = bigR*cos(bigTheta);
%     cy = bigR*sin(bigTheta);
% 
%     Pos = zeros(totalNodes,2);
% 
%     for h=1:numHubs
%         idx = hubIdx{h};
%         Pos(idx,1)=hubX+cx(h);
%         Pos(idx,2)=hubY+cy(h);
%     end
% end

function [interDeg] = count_additional_edges(A,nodesPerHub,numHubs,hubIdx)
    totalNodes = nodesPerHub*numHubs;
    interDeg = zeros(totalNodes,1);
    
    for h = 1:numHubs
        idx = hubIdx{h};
    
        % degree inside hub only
        degInternal = sum(A(idx,idx),2);
    
        % total degree
        degTotal = sum(A(idx,:),2);
    
        % extra inter-hub degree
        interDeg(idx) = degTotal - degInternal;
    end
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