function generate_cf_visualizations(A_cf, B_cf, run_number)
%% Generate visualizations for UP vs FP counterfactual analysis

if nargin < 3
    run_number = 1;
end

% Filter entered bidders only for bid analysis
B_entered = B_cf(B_cf.entered == 1, :);

% Get scale factor from auction data (assuming it's the same for all)
scale_s = A_cf.scale_s(1);  % All auctions should have same scale

%% 1. Bid Distribution Comparison using KDE
figure('Position', [100, 100, 800, 400]);

% Subplot 1: KDE comparison
subplot(1,2,1);
% Apply rescaling by dividing by scale_s
UP_scores = B_entered.score(B_entered.mechanism == "UP") / scale_s;
FP_scores = B_entered.score(B_entered.mechanism == "FP") / scale_s;

% Create KDE plots
[f_up, xi_up] = ksdensity(UP_scores);
[f_fp, xi_fp] = ksdensity(FP_scores);

plot(xi_up, f_up, 'LineWidth', 2, 'DisplayName', 'UP');
hold on;
plot(xi_fp, f_fp, 'LineWidth', 2, 'DisplayName', 'FP');
xlabel('Score/Bid'); ylabel('Density'); 
title('Bid Distribution (KDE)');
legend; grid on;

% Subplot 2: Risk state comparison using KDE
subplot(1,2,2);
% Create lookup table for auction risk states
A_risk_lookup = containers.Map('KeyType', 'double', 'ValueType', 'char');
for i = 1:height(A_cf)
    A_risk_lookup(A_cf.auction_id(i)) = char(A_cf.risk(i));
end

% Get risk state for each bidder
B_UP_entered = B_entered(B_entered.mechanism == "UP", :);
risk_states_UP = cellfun(@(x) A_risk_lookup(x), num2cell(B_UP_entered.auction_id), 'UniformOutput', false);
UP_L = B_UP_entered.score(strcmp(risk_states_UP, 'L')) / scale_s;
UP_H = B_UP_entered.score(strcmp(risk_states_UP, 'H')) / scale_s;

B_FP_entered = B_entered(B_entered.mechanism == "FP", :);
risk_states_FP = cellfun(@(x) A_risk_lookup(x), num2cell(B_FP_entered.auction_id), 'UniformOutput', false);
FP_L = B_FP_entered.score(strcmp(risk_states_FP, 'L')) / scale_s;
FP_H = B_FP_entered.score(strcmp(risk_states_FP, 'H')) / scale_s;

if ~isempty(UP_L) && ~isempty(UP_H)
    [f_up_l, xi_up_l] = ksdensity(UP_L);
    [f_up_h, xi_up_h] = ksdensity(UP_H);
    plot(xi_up_l, f_up_l, '--', 'LineWidth', 1.5, 'DisplayName', 'UP-Low Risk');
    hold on;
    plot(xi_up_h, f_up_h, '-', 'LineWidth', 1.5, 'DisplayName', 'UP-High Risk');
end

if ~isempty(FP_L) && ~isempty(FP_H)
    [f_fp_l, xi_fp_l] = ksdensity(FP_L);
    [f_fp_h, xi_fp_h] = ksdensity(FP_H);
    plot(xi_fp_l, f_fp_l, ':', 'LineWidth', 1.5, 'DisplayName', 'FP-Low Risk');
    plot(xi_fp_h, f_fp_h, '-.', 'LineWidth', 1.5, 'DisplayName', 'FP-High Risk');
end

xlabel('Score/Bid'); ylabel('Density'); 
title('Bid Distribution by Risk State');
legend; grid on;

filename = sprintf('output/fig/%02d_bid_distributions.png', run_number);
saveas(gcf, filename);
close;

%% 2. Box Plot (separate figure)
figure('Position', [100, 100, 600, 400]);
% Apply rescaling to boxplot data
boxplot([B_entered.score(B_entered.mechanism == "UP") / scale_s; ...
         B_entered.score(B_entered.mechanism == "FP") / scale_s], ...
        [repmat({'UP'}, sum(B_entered.mechanism == "UP"), 1); ...
         repmat({'FP'}, sum(B_entered.mechanism == "FP"), 1)]);
ylabel('Score/Bid'); title('Bid Distribution Box Plot');
grid on;

filename = sprintf('output/fig/%02d_bid_boxplot.png', run_number);
saveas(gcf, filename);
close;

%% 3. Entry and Competition Effects
A_UP = A_cf(A_cf.mechanism == "UP", :);
A_FP = A_cf(A_cf.mechanism == "FP", :);

figure('Position', [100, 100, 800, 400]);

subplot(1,2,1);
% Use solid filled histograms for entry distribution
histogram(A_UP.n, 'BinMethod', 'integers', 'Normalization', 'probability', ...
         'FaceAlpha', 1.0, 'EdgeColor', 'black', 'DisplayName', 'UP');
hold on;
histogram(A_FP.n, 'BinMethod', 'integers', 'Normalization', 'probability', ...
         'FaceAlpha', 0.7, 'EdgeColor', 'black', 'DisplayName', 'FP');
xlabel('Number of Entrants'); ylabel('Probability');
title('Entry Distribution');
legend; grid on;

subplot(1,2,2);
% Apply rescaling to final payments
scatter(A_UP.n, A_UP.pay_final / scale_s, 30, 'filled', 'DisplayName', 'UP', 'MarkerFaceAlpha', 0.6);
hold on;
scatter(A_FP.n, A_FP.pay_final / scale_s, 30, 'filled', 'DisplayName', 'FP', 'MarkerFaceAlpha', 0.6);
xlabel('Number of Entrants'); ylabel('Final Payment');
title('Competition vs Payment');
legend; grid on;

filename = sprintf('output/fig/%02d_entry_competition.png', run_number);
saveas(gcf, filename);
close;

fprintf('Visualizations saved to output/fig/ with run number %02d\n', run_number);

end