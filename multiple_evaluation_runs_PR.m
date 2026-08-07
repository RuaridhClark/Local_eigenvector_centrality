%% ================= PARAMETERS =================
clc;
close all;

nodesPerHub = 20;
numHubs     = 6;
pInternal   = 0.15;

nInterList  = 1:5;
nP          = length(nInterList);

seedList    = 1:10;
nSeeds      = length(seedList);

%% ================= STORAGE =================

fpList = [0.01 0.05:0.05:1.0];
nFP = length(fpList);

corr_local = zeros(nP,nSeeds);
corr_pr    = zeros(nP,nSeeds,nFP);

sym_local = zeros(nP,nSeeds);
sym_pr = zeros(nP,nSeeds,nFP);

std_local = zeros(nodesPerHub,1);
std_pr = zeros(nodesPerHub,nFP);

%% ================= MAIN SWEEP =================

for t = 1:nP

    nInter = nInterList(t);

    for s = 1:nSeeds

        seed = seedList(s);

        [A, Pos, hubIdx] = build_modular_network( ...
            nodesPerHub, ...
            numHubs, ...
            pInternal, ...
            nInter, ...
            seed);

        interDeg = count_additional_edges(A,nodesPerHub,numHubs,hubIdx);

        rowSum = sum(A, 2);
        A = A ./ rowSum;
        A(rowSum == 0, :) = 0;  

        G = digraph(A);

        totalNodes = size(A,1);

        normalize_vec = @(x) x(:)./max(x);

        %% -------- KCEC --------

        KCEC = zeros(totalNodes,1);

        for h = 1:numHubs

            idx = hubIdx{h};

            Ahub = A(idx,idx);

            % KCEC(idx) = centrality(digraph(Ahub),'eigenvector');
            [KCEC(idx),  ~] = local_eigenvector_centrality(Ahub', Pos, false, 1);

        end

        KCEC = normalize_vec(KCEC);

        %% -------- LEC --------

        LEC = local_eigenvector_centrality(A',Pos,false,6);
        LEC = normalize_vec(LEC);

        %% -------- GLOBALS --------

        DG = digraph(A);

        prAll = zeros(totalNodes,nFP);
        
        for fpIdx = 1:nFP
            fp = fpList(fpIdx);
            prAll(:,fpIdx) = normalize_vec( ...
                centrality(DG,'pagerank', ...
                'FollowProbability',fp));
        end

        %% -------- CORRELATION --------

        C = corrcoef(KCEC,LEC);
        corr_local(t,s) = C(1,2);
        
        for fpIdx = 1:nFP
        
            C = corrcoef(KCEC,prAll(:,fpIdx));
            corr_pr(t,s,fpIdx) = C(1,2);
        
        end

        %% Symmetry

        for k = 1:nodesPerHub

            L = zeros(numHubs,1);
        
            P = zeros(numHubs,nFP);
        
            for h = 1:numHubs
        
                idx = (h-1)*nodesPerHub + k;
        
                L(h) = LEC(idx);
        
                for fpIdx = 1:nFP
                    P(h,fpIdx) = prAll(idx,fpIdx);
                end
        
            end
        
            std_local(k) = std(L);
        
            for fpIdx = 1:nFP
                std_pr(k,fpIdx) = std(P(:,fpIdx));
            end
        
        end
        
        sym_local(t,s) = mean(std_local);
        
        for fpIdx = 1:nFP
            sym_pr(t,s,fpIdx) = mean(std_pr(:,fpIdx));
        end


    end

end

%% ================= SUMMARY =================

corr_local_mu = mean(corr_local,2);
corr_local_sd = std(corr_local,0,2);

corr_pr_mu = squeeze(mean(corr_pr,2));
corr_pr_sd = squeeze(std(corr_pr,0,2));

muEC   = mean(r_EC  ,2);
muPR   = mean(r_PR  ,2);

sdEC   = std(r_EC  ,0,2);
sdPR   = std(r_PR  ,0,2);

sym_local_mu = mean(sym_local,2);
sym_local_sd = std(sym_local,0,2);

sym_pr_mu = squeeze(mean(sym_pr,2));
sym_pr_sd = squeeze(std(sym_pr,0,2));

%% ================= CORRELATION PLOT =================

cols = lines(5);

figure;
hold on

%% --- LEC shaded line ---

shaded_line( ...
    nInterList,...
    corr_local_mu,...
    corr_local_sd,...
    cols(1,:), ...
    0.2);

%% --- PR shaded line ---
shaded_line( ...
    nInterList,...
    corr_pr_mu(:,1:nFP),...
    corr_pr_sd(:,1:nFP),...
    cols(4,:), ...
    0.05);

%% --- Individual PR curves ---

for fpIdx = 1:nFP

    plot( ...
        nInterList,...
        corr_pr_mu(:,fpIdx),...
        'Color',cols(4,:),...
        'LineWidth',1);

end

%% --- Text labels on curves ---

labelPos = length(nInterList);  % right-most point

for fpIdx = 1:nFP

    x = nInterList(labelPos);

    y = corr_pr_mu(labelPos,fpIdx);

    % text( ...
    %     x + 0.08,...
    %     y,...
    %     sprintf('%.2f',fpList(fpIdx)),...
    %     'Color',cols(fpIdx,:),...
    %     'FontSize',8,...
    %     'HorizontalAlignment','left');

end

labelIdx = [1 6 11 16 21];

for fpIdx = labelIdx

    text( ...
        nInterList(end)+0.08,...
        corr_pr_mu(end,fpIdx),...
        sprintf('%.2f',fpList(fpIdx)),...
        'Color',cols(4,:),...
        'FontWeight','bold');
end

xlabel('No. of inter-community connections','FontSize',12)
ylabel('Correlation with KCEC','FontSize',12)

xlim([min(nInterList) max(nInterList)+0.8])

xticks(1:5)

grid on
box on

axis tight

%% Sym plot

%% ================= SYMMETRY PLOT =================

cols = lines(5);

figure;
hold on



% %% --- PR envelope (all FollowProbability values) ---
% 
% prMin = min(sym_pr_mu,[],2);
% prMax = max(sym_pr_mu,[],2);
% 
% fill([nInterList fliplr(nInterList)], ...
%      [prMin' fliplr(prMax')], ...
%      cols(4,:), ...
%      'FaceAlpha',0.05, ...
%      'EdgeColor','none');

%% --- PR shaded line ---
shaded_line( ...
    nInterList,...
    sym_pr_mu(:,1:nFP),...
    sym_pr_sd(:,1:nFP),...
    cols(4,:), ...
    0.1);

%% --- Individual PR curves ---

for fpIdx = 1:nFP

    plot( ...
        nInterList,...
        sym_pr_mu(:,fpIdx),...
        'Color',cols(4,:),...
        'LineWidth',1);

end

%% --- LEC shaded line ---

shaded_line( ...
    nInterList,...
    sym_local_mu,...
    sym_local_sd,...
    cols(1,:), ...
    0.2);

%% --- Labels for selected PR curves ---

labelIdx = [1 6 11 16 21];

for fpIdx = labelIdx

    text( ...
        nInterList(end)+0.08,...
        sym_pr_mu(end,fpIdx),...
        sprintf('%.2f',fpList(fpIdx)),...
        'Color',cols(4,:),...
        'FontWeight','bold');

end

xlabel('No. of inter-hub connections','FontSize',12)
ylabel('Node centrality variation','FontSize',12)

xlim([min(nInterList) max(nInterList)+0.8])

xticks(1:5)

grid on
box on

axis tight
ylim([0 0.15])
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

function shaded_line(x,mu,sigma,col,alpha)

    fill([x fliplr(x)], ...
         [(mu-sigma)' fliplr((mu+sigma)')], ...
         col,...
         'FaceAlpha',alpha,...
         'EdgeColor','none');

    hold on

    plot(x,mu,...
        'Color',col,...
        'LineWidth',2);

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

    %% --- geometry ---
    theta = linspace(0,2*pi,nodesPerHub+1)';
    theta(end)=[];

    hubR = 2;

    hubX = hubR*cos(theta);
    hubY = hubR*sin(theta);

    bigTheta = linspace(0,2*pi,numHubs+1)';
    bigTheta(end)=[];

    bigR = 6;

    cx = bigR*cos(bigTheta);
    cy = bigR*sin(bigTheta);

    Pos = zeros(totalNodes,2);

    for h=1:numHubs
        idx = hubIdx{h};
        Pos(idx,1)=hubX+cx(h);
        Pos(idx,2)=hubY+cy(h);
    end
end

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