%% ================= PARAMETERS =================
clc;
% close all;

nodesPerComm = 20;
numComms     = 6;
pInternal   = 0.15;

nInterList  = 1:5;
nP          = length(nInterList);

seedList    = 1:10;
nSeeds      = length(seedList);

%% ================= STORAGE =================

corr_local  = zeros(nP,nSeeds);
corr_global = zeros(nP,nSeeds);
corr_katz   = zeros(nP,nSeeds);
corr_pr     = zeros(nP,nSeeds);
corr_pcc    = zeros(nP,nSeeds);

sym_local   = zeros(nP,nSeeds);
sym_global  = zeros(nP,nSeeds);
sym_katz    = zeros(nP,nSeeds);
sym_pr      = zeros(nP,nSeeds);
sym_pcc     = zeros(nP,nSeeds);

r_EC   = zeros(nP,nSeeds);
r_PR   = zeros(nP,nSeeds);

%% ================= MAIN SWEEP =================

for t = 1:nP

    nInter = nInterList(t);

    for s = 1:nSeeds

        seed = seedList(s);

        [A, Pos, CommIdx] = build_modular_network( ...
            nodesPerComm, ...
            numComms, ...
            pInternal, ...
            nInter, ...
            seed);

        interDeg = count_additional_edges(A,nodesPerComm,numComms,CommIdx);

        G = graph(A);

        totalNodes = size(A,1);

        normalize_vec = @(x) abs(x(:))./sum(abs(x));

        %% -------- KCEC --------

        KCEC = zeros(totalNodes,1);

        for h = 1:numComms

            idx = CommIdx{h};

            AComm = A(idx,idx);

            KCEC(idx) = centrality(graph(AComm),'eigenvector');

        end

        KCEC = normalize_vec(KCEC);

        %% -------- LEC --------

        LEC = local_eigenvector_centrality(A,Pos,false,6);
        LEC = normalize_vec(LEC);

        %% -------- GLOBALS --------

        ec = normalize_vec(centrality(G,'eigenvector'));

        DG = digraph(A);
        pr = normalize_vec( ...
            centrality(DG,'pagerank', ...
            'FollowProbability',0.85));

        lambda_max = max(abs(eigs(sparse(A),1)));

        alpha = 0.85/lambda_max;

        katz = normalize_vec( ...
            (speye(totalNodes)-alpha*A')\ones(totalNodes,1));

        pcc = normalize_vec( ...
            principal_component_centrality(A));

        %% -------- CORRELATION --------

        C = corrcoef(KCEC, LEC);
        corr_local(t,s) = C(1,2);
        C = corrcoef(KCEC,ec);
        corr_global(t,s) = C(1,2);
        C = corrcoef(KCEC,katz);
        corr_katz(t,s) = C(1,2);
        C = corrcoef(KCEC,pr);
        corr_pr(t,s) = C(1,2);
        C = corrcoef(KCEC,pcc);
        corr_pcc(t,s) = C(1,2);

        %% ------- GLOBAL NODE CORR -----

        deltaLEC = LEC - KCEC;
        deltaEC = ec - KCEC;
        deltaPR = pr - KCEC;

        R = corrcoef(deltaLEC,ec);
        r_EC(t,s) = corr(deltaLEC,deltaEC,'Type','Spearman');
        r_PR(t,s) = corr(deltaLEC,deltaPR,'Type','Spearman');

        % %% Plot detlaLEC against ec
        % figure;
        % scatter(deltaEC, deltaLEC, 'filled');
        % xlabel('Eigenvector Centrality (EC)');
        % ylabel('Delta Local Eigenvector Centrality (deltaLEC)');
        % title(sprintf('Delta LEC vs EC for nInter = %d, seed = %d', nInter, seed));
        % grid on;


        %% -------- SYMMETRY --------

        std_local  = zeros(nodesPerComm,1);
        std_global = zeros(nodesPerComm,1);
        std_katz   = zeros(nodesPerComm,1);
        std_pr     = zeros(nodesPerComm,1);
        std_pcc    = zeros(nodesPerComm,1);

        for k = 1:nodesPerComm

            L = zeros(numComms,1);
            Gv = zeros(numComms,1);
            K = zeros(numComms,1);
            P = zeros(numComms,1);
            PC = zeros(numComms,1);

            for h = 1:numComms

                idx = (h-1)*nodesPerComm + k;

                L(h)  = LEC(idx);
                Gv(h) = ec(idx);
                K(h)  = katz(idx);
                P(h)  = pr(idx);
                PC(h) = pcc(idx);

            end

            std_local(k)  = std(L);
            std_global(k) = std(Gv);
            std_katz(k)   = std(K);
            std_pr(k)     = std(P);
            std_pcc(k)    = std(PC);

        end

        sym_local(t,s)  = mean(std_local);
        sym_global(t,s) = mean(std_global);
        sym_katz(t,s)   = mean(std_katz);
        sym_pr(t,s)     = mean(std_pr);
        sym_pcc(t,s)    = mean(std_pcc);

    end

end

%% ================= SUMMARY =================

corr_local_mu  = mean(corr_local,2);
corr_global_mu = mean(corr_global,2);
corr_katz_mu   = mean(corr_katz,2);
corr_pr_mu     = mean(corr_pr,2);
corr_pcc_mu    = mean(corr_pcc,2);

corr_local_sd  = std(corr_local,0,2);
corr_global_sd = std(corr_global,0,2);
corr_katz_sd   = std(corr_katz,0,2);
corr_pr_sd     = std(corr_pr,0,2);
corr_pcc_sd    = std(corr_pcc,0,2);

sym_local_mu  = mean(sym_local,2);
sym_global_mu = mean(sym_global,2);
sym_katz_mu   = mean(sym_katz,2);
sym_pr_mu     = mean(sym_pr,2);
sym_pcc_mu    = mean(sym_pcc,2);

sym_local_sd  = std(sym_local,0,2);
sym_global_sd = std(sym_global,0,2);
sym_katz_sd   = std(sym_katz,0,2);
sym_pr_sd     = std(sym_pr,0,2);
sym_pcc_sd    = std(sym_pcc,0,2);

muEC   = mean(r_EC  ,2);
muPR   = mean(r_PR  ,2);

sdEC   = std(r_EC  ,0,2);
sdPR   = std(r_PR  ,0,2);

%% ================= CORRELATION PLOT =================

figure;
hold on

cols = lines(5);

shaded_line(nInterList,corr_local_mu ,corr_local_sd ,cols(1,:));
shaded_line(nInterList,corr_global_mu,corr_global_sd,cols(2,:));
shaded_line(nInterList,corr_katz_mu  ,corr_katz_sd  ,cols(3,:));
shaded_line(nInterList,corr_pr_mu    ,corr_pr_sd    ,cols(4,:));
shaded_line(nInterList,corr_pcc_mu   ,corr_pcc_sd   ,cols(5,:));

xlabel('No. of inter-community connections','FontSize',14);
ylabel('Correlation with KCEC','FontSize',14);

% legend({'LEC','','EC','','Katz','','PageRank','','PCC',''},"Location","southoutside","Orientation","horizontal");
grid on

%% ================= SYMMETRY PLOT =================

figure;
hold on

shaded_line(nInterList,sym_local_mu ,sym_local_sd ,cols(1,:));
shaded_line(nInterList,sym_global_mu,sym_global_sd,cols(2,:));
shaded_line(nInterList,sym_katz_mu  ,sym_katz_sd  ,cols(3,:));
shaded_line(nInterList,sym_pr_mu    ,sym_pr_sd    ,cols(4,:));
shaded_line(nInterList,sym_pcc_mu   ,sym_pcc_sd   ,cols(5,:));

xlabel('No. of inter-community connections','FontSize',14);
ylabel('Community node variation','FontSize',14);

% legend({'','LEC','','EC','','Katz','','PageRank','','PCC'},"Location","southoutside","Orientation","horizontal");
grid on
xticks(1:5)

figure

methods = {'LEC','EC','Katz','PR','PCC'};

for k = 1:nP

    subplot(1,nP,k)
    hold on

    X = [ ...
        corr_local(k,:)' ...
        corr_global(k,:)' ...
        corr_katz(k,:)' ...
        corr_pr(k,:)' ...
        corr_pcc(k,:)' ];

    for m = 1:5

        vals = X(:,m);

        % horizontal jitter
        x = m + 0.08*randn(size(vals));

        % individual seeds
        scatter(x, vals, 25, ...
            'filled', ...
            'MarkerFaceAlpha',0.6);

        % mean ± std
        mu = mean(vals);
        sd = std(vals);

        errorbar(m, mu, sd, ...
            'k', ...
            'LineWidth',2, ...
            'CapSize',10);

        % mean marker
        plot(m, mu, 'ks', ...
            'MarkerFaceColor','y', ...
            'MarkerSize',6);
    end

    xlim([0.5 5.5])
    ylim([0 1])

    xticks(1:5)
    xticklabels(methods)

    title(sprintf('nInter = %g', nInterList(k)))

    if k == 1
        ylabel('Correlation with KCEC')
    end

    grid on

end

%% ============= DISAGREE PLOT ==============

figure
hold on

shaded_line(nInterList,muEC  ,sdEC  ,cols(2,:));
shaded_line(nInterList,muPR  ,sdPR  ,cols(4,:));

xlabel('No. of inter-community connections','FontSize',14)
ylabel('KCEC difference correlation','FontSize',14)
% legend({'LEC','','EC','','Katz','','PR','','PCC',''})
grid on
xticks(1:5)

% %% ================= BOXPLOTS =================
% 
% figure
% 
% for k = 1:nP
% 
%     subplot(1,nP,k)
% 
%     X = [ ...
%         corr_local(k,:)' ...
%         corr_global(k,:)' ...
%         corr_katz(k,:)' ...
%         corr_pr(k,:)' ...
%         corr_pcc(k,:)' ];
% 
%     boxplot(X,...
%         'Labels',{'LEC','EC','Katz','PR','PCC'});
% 
%     title(sprintf('nInter=%d',nInterList(k)))
% 
%     ylim([0 1])
% 
% end

%% ================= HELPER =================

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

function [A, Pos, CommIdx] = build_modular_network( ...
    nodesPerComm, numComms, pInternal, nInter, seed)

    % Optional argument
    if nargin >= 5 && ~isempty(seed)
        rng(seed);
    end

    totalNodes = nodesPerComm * numComms;

    %% --- base connected community ---
    baseAdj = zeros(nodesPerComm);

    perm = randperm(nodesPerComm);
    for i = 2:nodesPerComm
        n1 = perm(i);
        n2 = perm(randi(i-1));
        baseAdj(n1,n2)=1;
        baseAdj(n2,n1)=1;
    end

    extra = rand(nodesPerComm) < pInternal;
    extra = triu(extra,1);
    extra = extra + extra.';
    baseAdj = baseAdj | extra;
    baseAdj(1:nodesPerComm+1:end)=0;

    %% --- replicate ---
    A = zeros(totalNodes);
    CommIdx = cell(numComms,1);

    for h = 1:numComms
        idx = (h-1)*nodesPerComm + (1:nodesPerComm);
        A(idx,idx)=baseAdj;
        CommIdx{h}=idx;
    end

    %% --- inter community edges ---
    for k = 1:nInter
        for h1=1:numComms-1
            for h2=h1+1:numComms
                i1 = CommIdx{h1}(randi(nodesPerComm));
                i2 = CommIdx{h2}(randi(nodesPerComm));

                A(i1,i2)=1;
                A(i2,i1)=1;
            end
        end
    end

    %% --- geometry ---
    theta = linspace(0,2*pi,nodesPerComm+1)';
    theta(end)=[];

    ComR = 2;

    ComX = ComR*cos(theta);
    ComY = ComR*sin(theta);

    bigTheta = linspace(0,2*pi,numComms+1)';
    bigTheta(end)=[];

    bigR = 6;

    cx = bigR*cos(bigTheta);
    cy = bigR*sin(bigTheta);

    Pos = zeros(totalNodes,2);

    for h=1:numComms
        idx = CommIdx{h};
        Pos(idx,1)=ComX+cx(h);
        Pos(idx,2)=ComY+cy(h);
    end
end

function [interDeg] = count_additional_edges(A,nodesPerComm,numComms,CommIdx)
    totalNodes = nodesPerComm*numComms;
    interDeg = zeros(totalNodes,1);
    
    for h = 1:numComms
        idx = CommIdx{h};
    
        % degree inside community only
        degInternal = sum(A(idx,idx),2);
    
        % total degree
        degTotal = sum(A(idx,:),2);
    
        % extra inter-community degree
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