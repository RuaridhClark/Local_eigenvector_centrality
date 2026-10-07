% Data
Transportation = [1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.01];

warped   = [0.740071597 0.773477041 0.79096287 0.804692762 0.810352858 ...
            0.813315389 0.816690696 0.816420607 0.815171404 0.812637322 0.808250572];

global_  = [0.993678608 0.978317117 0.961886072 0.946287504 0.93371644 ...
            0.923817722 0.91194892 0.902608271 0.891722253 0.882315208 0.872117213];

local2   = [0.935378776 0.919735877 0.903471152 0.88825399 0.876173817 ...
            0.866826554 0.855616614 0.846910326 0.836736379 0.828034168 0.818589762];

local11  = [0.698274034 0.728040587 0.741700174 0.748298604 0.750268546 ...
            0.74997165 0.748696853 0.746294928 0.743122304 0.739308207 0.734834453];

classAcc = [0.740921788 0.754330279 0.753411765 0.746595638 0.737159461 ...
            0.726931469 0.71518268 0.704513687 0.692679519 0.681549351 0.669805272];

% Plot
figure; hold on; grid on;
% plot(Transportation, warped,   '-o', 'LineWidth', 2);
plot(Transportation, global_,  '-s', 'LineWidth', 2);
plot(Transportation, local2,   '-^', 'LineWidth', 2);
plot(Transportation, local11,  '-d', 'LineWidth', 2);
plot(Transportation, classAcc,'-x', 'LineWidth', 2);

% Formatting
% set(gca, 'XDir', 'reverse');  % common for transportation parameter
xlabel('$\alpha$','Interpreter','latex','FontSize',18);
ylabel('$c$','Interpreter','latex','FontSize',18);
% title('Performance vs Teleportation factor');
legend({'Eigenvector centrality','Local ($i=2$)','Local ($i=11$)','Class'},'Interpreter','latex', ...
       'Location', 'northwest');

ylim([0.65 1.0]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2);
box on;
