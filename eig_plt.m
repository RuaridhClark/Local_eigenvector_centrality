%%
load('eigenvalues_Chicago.mat')

% Extract and sort eigenvalues in descending order
[Dsort, idx] = sort(eigenvalues, 'descend'); 

% Compute eigengap
eigengap = Dsort(1:end-2) - Dsort(2:end-1);

% Find peaks
[eigenpeaks, locs] = findpeaks(eigengap, 1:length(eigengap));

% Plot eigengap
figure;
plot(Dsort(:));

% Plot eigengap
figure;
plot(eigengap,'k-o');
xlabel('$i$','interpreter','latex','FontSize',14)
ylabel('${\lambda}_i - {\lambda}_{i+1}$','interpreter','latex','FontSize',14)
box off

%%
load('eigenvalues_Glasgow_nL.mat')

% Extract and sort eigenvalues in descending order
[Dsort, idx] = sort(eigenvalues, 'ascend'); 

% Compute eigengap
eigengap = abs(Dsort(1:end-2) - Dsort(2:end-1));

% Find peaks
[eigenpeaks, locs] = findpeaks(eigengap, 1:length(eigengap));

% Plot eigengap
figure;
plot(Dsort(:));

% Plot eigengap
figure;
plot(eigengap,'k-o');
xlabel('$i$','interpreter','latex','FontSize',14)
ylabel('${\lambda}_i - {\lambda}_{i+1}$','interpreter','latex','FontSize',14)
box off