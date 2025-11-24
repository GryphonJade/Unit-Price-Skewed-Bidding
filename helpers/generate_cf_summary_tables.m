function generate_cf_summary_tables(A_cf, B_cf, run_number)
%% Generate summary tables for UP vs FP counterfactual analysis

if nargin < 3
    run_number = 1;
end

% Filter entered bidders only
B_entered = B_cf(B_cf.entered == 1, :);

% Get scale factor for rescaling monetary values
scale_s = A_cf.scale_s(1);  % All auctions should have same scale

%% Table 1: Basic Bid Statistics Comparison
fprintf('\n=== Table 1: Basic Bid Statistics ===\n');

% Separate UP and FP data
B_UP = B_entered(B_entered.mechanism == "UP", :);
B_FP = B_entered(B_entered.mechanism == "FP", :);
A_UP = A_cf(A_cf.mechanism == "UP", :);
A_FP = A_cf(A_cf.mechanism == "FP", :);

% Basic statistics (apply rescaling)
stats_basic = table();
stats_basic.Mechanism = {'UP'; 'FP'};
stats_basic.Mean_Bid = [mean(B_UP.score) / scale_s; mean(B_FP.score) / scale_s];
stats_basic.Std_Bid = [std(B_UP.score) / scale_s; std(B_FP.score) / scale_s];
stats_basic.Mean_Winning_Bid = [mean(A_UP.win_score) / scale_s; mean(A_FP.win_score) / scale_s];
stats_basic.Std_Winning_Bid = [std(A_UP.win_score) / scale_s; std(A_FP.win_score) / scale_s];
stats_basic.Mean_Final_Payment = [mean(A_UP.pay_final) / scale_s; mean(A_FP.pay_final) / scale_s];
stats_basic.Std_Final_Payment = [std(A_UP.pay_final) / scale_s; std(A_FP.pay_final) / scale_s];

disp(stats_basic);
filename = sprintf('output/tab/%02d_basic_bid_statistics.csv', run_number);
writetable(stats_basic, filename);

%% Table 2: UP Component Bid Statistics  
fprintf('\n=== Table 2: UP Component Bid Statistics ===\n');

% By risk state
B_UP_L = B_UP(ismember(B_UP.auction_id, A_UP.auction_id(A_UP.risk == "L")), :);
B_UP_H = B_UP(ismember(B_UP.auction_id, A_UP.auction_id(A_UP.risk == "H")), :);

stats_up_components = table();
stats_up_components.Risk_State = {'Low'; 'High'; 'Pooled'};
stats_up_components.Mean_s = [mean(B_UP_L.score) / scale_s; mean(B_UP_H.score) / scale_s; mean(B_UP.score) / scale_s];
stats_up_components.Std_s = [std(B_UP_L.score) / scale_s; std(B_UP_H.score) / scale_s; std(B_UP.score) / scale_s];
stats_up_components.Mean_b0 = [mean(B_UP_L.b0) / scale_s; mean(B_UP_H.b0) / scale_s; mean(B_UP.b0) / scale_s];
stats_up_components.Std_b0 = [std(B_UP_L.b0) / scale_s; std(B_UP_H.b0) / scale_s; std(B_UP.b0) / scale_s];
stats_up_components.Mean_b2 = [mean(B_UP_L.b2) / scale_s; mean(B_UP_H.b2) / scale_s; mean(B_UP.b2) / scale_s];
stats_up_components.Std_b2 = [std(B_UP_L.b2) / scale_s; std(B_UP_H.b2) / scale_s; std(B_UP.b2) / scale_s];
stats_up_components.Corr_b0_b2 = [corr(B_UP_L.b0, B_UP_L.b2); corr(B_UP_H.b0, B_UP_H.b2); corr(B_UP.b0, B_UP.b2)];

disp(stats_up_components);
filename = sprintf('output/tab/%02d_up_component_statistics.csv', run_number);
writetable(stats_up_components, filename);

%% Table 3: UP to FP Conversion Effects (Changes in %)
fprintf('\n=== Table 3: UP to FP Conversion Effects ===\n');

% Match auctions by auction_id (assuming same auction simulated under both mechanisms)
n_auctions = height(A_UP);
conversion_effects = table();
conversion_effects.Auction_ID = A_UP.auction_id;
conversion_effects.Risk_State = A_UP.risk;

% Calculate percentage changes (percentages don't need rescaling as they are ratios)
conversion_effects.Winning_Bid_Change_Pct = 100 * (A_FP.win_score - A_UP.win_score) ./ A_UP.win_score;
conversion_effects.Final_Payment_Change_Pct = 100 * (A_FP.pay_final - A_UP.pay_final) ./ A_UP.pay_final;
conversion_effects.Entry_Change = A_FP.n - A_UP.n;
conversion_effects.Overrun_Change = A_FP.overrun - A_UP.overrun;

% Summary statistics of changes
summary_changes = table();
summary_changes.Metric = {'Mean_Winning_Bid_Change_Pct'; 'Mean_Final_Payment_Change_Pct'; 'Mean_Entry_Change'; 'Mean_Overrun_Change'};
summary_changes.All = [mean(conversion_effects.Winning_Bid_Change_Pct); 
                      mean(conversion_effects.Final_Payment_Change_Pct);
                      mean(conversion_effects.Entry_Change);
                      mean(conversion_effects.Overrun_Change)];
summary_changes.Low_Risk = [mean(conversion_effects.Winning_Bid_Change_Pct(conversion_effects.Risk_State == "L"));
                           mean(conversion_effects.Final_Payment_Change_Pct(conversion_effects.Risk_State == "L"));
                           mean(conversion_effects.Entry_Change(conversion_effects.Risk_State == "L"));
                           mean(conversion_effects.Overrun_Change(conversion_effects.Risk_State == "L"))];
summary_changes.High_Risk = [mean(conversion_effects.Winning_Bid_Change_Pct(conversion_effects.Risk_State == "H"));
                            mean(conversion_effects.Final_Payment_Change_Pct(conversion_effects.Risk_State == "H"));
                            mean(conversion_effects.Entry_Change(conversion_effects.Risk_State == "H"));
                            mean(conversion_effects.Overrun_Change(conversion_effects.Risk_State == "H"))];

disp(summary_changes);
filename = sprintf('output/tab/%02d_conversion_effects_detailed.csv', run_number);
writetable(conversion_effects, filename);
filename = sprintf('output/tab/%02d_conversion_effects_summary.csv', run_number);
writetable(summary_changes, filename);

%% Table 4: Risk State Analysis
fprintf('\n=== Table 4: Risk State Effects ===\n');

risk_analysis = table();
risk_analysis.Mechanism = {'UP_Low'; 'UP_High'; 'FP_Low'; 'FP_High'};

A_UP_L = A_UP(A_UP.risk == "L", :);
A_UP_H = A_UP(A_UP.risk == "H", :);
A_FP_L = A_FP(A_FP.risk == "L", :);
A_FP_H = A_FP(A_FP.risk == "H", :);

risk_analysis.N_Auctions = [height(A_UP_L); height(A_UP_H); height(A_FP_L); height(A_FP_H)];
risk_analysis.Mean_Entry = [mean(A_UP_L.n); mean(A_UP_H.n); mean(A_FP_L.n); mean(A_FP_H.n)];
risk_analysis.Mean_Winning_Bid = [mean(A_UP_L.win_score) / scale_s; mean(A_UP_H.win_score) / scale_s; mean(A_FP_L.win_score) / scale_s; mean(A_FP_H.win_score) / scale_s];
risk_analysis.Mean_Final_Payment = [mean(A_UP_L.pay_final) / scale_s; mean(A_UP_H.pay_final) / scale_s; mean(A_FP_L.pay_final) / scale_s; mean(A_FP_H.pay_final) / scale_s];
risk_analysis.Std_Final_Payment = [std(A_UP_L.pay_final) / scale_s; std(A_UP_H.pay_final) / scale_s; std(A_FP_L.pay_final) / scale_s; std(A_FP_H.pay_final) / scale_s];
risk_analysis.Mean_Abs_Overrun = [mean(abs(A_UP_L.overrun)); mean(abs(A_UP_H.overrun)); mean(abs(A_FP_L.overrun)); mean(abs(A_FP_H.overrun))];

disp(risk_analysis);
filename = sprintf('output/tab/%02d_risk_state_analysis.csv', run_number);
writetable(risk_analysis, filename);

%% Table 5: Entry and Competition Analysis
fprintf('\n=== Table 5: Entry and Competition Analysis ===\n');

% Pre-create table with proper structure
entry_analysis = table('Size', [0, 6], ...
    'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'N_Entrants', 'Mechanism', 'N_Obs', 'Mean_Winning_Bid', 'Mean_Final_Payment', 'Std_Final_Payment'});

entry_levels = [2; 3; 4; 5; 6]; % Common entry levels

for i = 1:length(entry_levels)
    n_level = entry_levels(i);
    
    % UP statistics for this entry level
    A_UP_n = A_UP(A_UP.n == n_level, :);
    A_FP_n = A_FP(A_FP.n == n_level, :);
    
    if height(A_UP_n) > 0 && height(A_FP_n) > 0
        entry_analysis = [entry_analysis; {n_level, "UP", height(A_UP_n), mean(A_UP_n.win_score) / scale_s, mean(A_UP_n.pay_final) / scale_s, std(A_UP_n.pay_final) / scale_s}];
        entry_analysis = [entry_analysis; {n_level, "FP", height(A_FP_n), mean(A_FP_n.win_score) / scale_s, mean(A_FP_n.pay_final) / scale_s, std(A_FP_n.pay_final) / scale_s}];
    end
end
disp(entry_analysis);
filename = sprintf('output/tab/%02d_entry_competition_analysis.csv', run_number);
writetable(entry_analysis, filename);

%% Table 6: Efficiency and Welfare Metrics
fprintf('\n=== Table 6: Efficiency and Welfare Metrics ===\n');

welfare_metrics = table();
welfare_metrics.Mechanism = {'UP'; 'FP'};

% Payment efficiency (lower is better for procurer) - apply rescaling
welfare_metrics.Mean_Payment = [mean(A_UP.pay_final) / scale_s; mean(A_FP.pay_final) / scale_s];
welfare_metrics.Payment_CV = [std(A_UP.pay_final)/mean(A_UP.pay_final); std(A_FP.pay_final)/mean(A_FP.pay_final)];  % CV is ratio, no rescaling needed

% Entry rates (no rescaling needed - these are counts/rates)
welfare_metrics.Mean_Entry_Rate = [mean(A_UP.n)/mean(A_UP.N_potential); mean(A_FP.n)/mean(A_FP.N_potential)];

% Risk metrics (no rescaling needed for overrun)
welfare_metrics.Mean_Abs_Overrun = [mean(abs(A_UP.overrun)); mean(abs(A_FP.overrun))];
welfare_metrics.Overrun_Risk = [std(A_UP.overrun); std(A_FP.overrun)];

% Winner selection (lower winning score = more efficient in cost terms) - apply rescaling
welfare_metrics.Mean_Winner_Score = [mean(A_UP.win_score) / scale_s; mean(A_FP.win_score) / scale_s];

disp(welfare_metrics);
filename = sprintf('output/tab/%02d_welfare_efficiency_metrics.csv', run_number);
writetable(welfare_metrics, filename);

%% Print key insights
fprintf('\n=== KEY INSIGHTS ===\n');
fprintf('1. Payment Change (UP->FP): %.2f%%\n', mean(conversion_effects.Final_Payment_Change_Pct));
fprintf('2. UP attracts %.1f entrants on average vs %.1f for FP\n', mean(A_UP.n), mean(A_FP.n));
fprintf('3. Mean overrun: UP=%.4f, FP=%.4f\n', mean(A_UP.overrun), mean(A_FP.overrun));
fprintf('4. Payment volatility (CV): UP=%.3f, FP=%.3f\n', std(A_UP.pay_final)/mean(A_UP.pay_final), std(A_FP.pay_final)/mean(A_FP.pay_final));
fprintf('5. Mean final payment (rescaled): UP=%.4f, FP=%.4f\n', mean(A_UP.pay_final)/scale_s, mean(A_FP.pay_final)/scale_s);

fprintf('\nSummary tables saved to output/tab/ with run number %02d\n', run_number);

end