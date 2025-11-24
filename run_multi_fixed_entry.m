%% Multi Fixed-entry Comparison: n = 3,5,7,9
clear; clc;

addpath(genpath('config'));
addpath(genpath('utils'));
addpath(genpath('strategy'));
addpath(genpath('sim'));
addpath(genpath('stats'));
addpath(genpath('counterfactual'));

fprintf('Multi Fixed-entry Comparison Simulation...\n');

% Entry levels to test
n_levels = [3, 5, 7, 9];
all_results = [];

% Run simulations for each entry level
for i = 1:length(n_levels)
    n_target = n_levels(i);
    fprintf('Running simulations for n = %d...\n', n_target);
    
    % Call the existing function and capture results
    [results, bidder_stats] = cf_run_fixed_entry_silent(n_target);
    
    % Add n_level column to results
    results.n_level = repmat(n_target, height(results), 1);
    
    % Store bidder statistics for Table 5
    bidder_stats.n_level = n_target;
    if i == 1
        all_bidder_stats = bidder_stats;
    else
        all_bidder_stats = [all_bidder_stats; bidder_stats];
    end
    
    % Append to combined results
    all_results = [all_results; results];
end

% Create output directories
if ~exist('output','dir'); mkdir('output'); end
if ~exist('output/tab','dir'); mkdir('output/tab'); end

% Save complete results
writetable(all_results, 'output/multi_fixed_entry_complete.csv');

% Generate summary table pooled across risk states
fprintf('\n=== Summary Table: Entry Level Effects (Pooled Across Risk States) ===\n');

summary_table = table();
summary_table.N_Entrants = n_levels';
summary_table.UP_Mean_Winning_Score = zeros(length(n_levels), 1);
summary_table.UP_Mean_Final_Payment = zeros(length(n_levels), 1);
summary_table.FP_Mean_Winning_Score = zeros(length(n_levels), 1);
summary_table.FP_Mean_Final_Payment = zeros(length(n_levels), 1);
summary_table.Winning_Score_Change_Pct = zeros(length(n_levels), 1);
summary_table.Final_Payment_Change_Pct = zeros(length(n_levels), 1);

for i = 1:length(n_levels)
    n_level = n_levels(i);
    subset = all_results(all_results.n_level == n_level, :);
    
    % Calculate pooled means
    summary_table.UP_Mean_Winning_Score(i) = mean(subset.up_win_score);
    summary_table.UP_Mean_Final_Payment(i) = mean(subset.up_pay_final);
    summary_table.FP_Mean_Winning_Score(i) = mean(subset.fp_win_score);
    summary_table.FP_Mean_Final_Payment(i) = mean(subset.fp_pay_final);
    
    % Calculate change percentages (UP -> FP)
    summary_table.Winning_Score_Change_Pct(i) = 100 * (mean(subset.fp_win_score) - mean(subset.up_win_score)) / mean(subset.up_win_score);
    summary_table.Final_Payment_Change_Pct(i) = 100 * (mean(subset.fp_pay_final) - mean(subset.up_pay_final)) / mean(subset.up_pay_final);
end

disp(summary_table);
writetable(summary_table, 'output/tab/multi_fixed_entry_summary.csv');

% Generate detailed breakdown by risk state
fprintf('\n=== Detailed Table: By Risk State ===\n');

detailed_table = table();
row_count = 0;

for i = 1:length(n_levels)
    n_level = n_levels(i);
    for risk = {'L', 'H'}
        risk_char = risk{1};
        subset = all_results(all_results.n_level == n_level & all_results.risk == risk_char, :);
        
        if height(subset) > 0
            row_count = row_count + 1;
            
            detailed_table.N_Entrants(row_count) = n_level;
            detailed_table.Risk_State(row_count) = string(risk_char);
            detailed_table.N_Auctions(row_count) = height(subset);
            detailed_table.UP_Mean_Winning_Score(row_count) = mean(subset.up_win_score);
            detailed_table.UP_Mean_Final_Payment(row_count) = mean(subset.up_pay_final);
            detailed_table.FP_Mean_Winning_Score(row_count) = mean(subset.fp_win_score);
            detailed_table.FP_Mean_Final_Payment(row_count) = mean(subset.fp_pay_final);
            detailed_table.Winning_Score_Change_Pct(row_count) = 100 * (mean(subset.fp_win_score) - mean(subset.up_win_score)) / mean(subset.up_win_score);
            detailed_table.Final_Payment_Change_Pct(row_count) = 100 * (mean(subset.fp_pay_final) - mean(subset.up_pay_final)) / mean(subset.up_pay_final);
        end
    end
end

disp(detailed_table);
writetable(detailed_table, 'output/tab/multi_fixed_entry_detailed.csv');

% Generate Table 6 style summary (similar to paper's Table 6)
fprintf('\n=== Table 6: Distribution of Simulated Bids (UP mechanism) ===\n');

table6_style = table();
table6_style.N_Entrants = n_levels';
table6_style.Score_Mean = zeros(length(n_levels), 1);
table6_style.Score_Std = zeros(length(n_levels), 1);
table6_style.B1_Mean = zeros(length(n_levels), 1);
table6_style.B1_Std = zeros(length(n_levels), 1);

% Use the collected bidder statistics
for i = 1:length(n_levels)
    n_level = n_levels(i);
    % Find the matching statistics for this n_level
    stats_idx = find([all_bidder_stats.n_level] == n_level);
    
    if ~isempty(stats_idx)
        stats_subset = all_bidder_stats(stats_idx(1));  % Take first match
        table6_style.Score_Mean(i) = stats_subset.score_mean;
        table6_style.Score_Std(i) = stats_subset.score_std;
        table6_style.B1_Mean(i) = stats_subset.b1_mean;
        table6_style.B1_Std(i) = stats_subset.b1_std;
    end
end

% Display in format similar to Table 6
fprintf('\nNumber of bidders: n    UP (Score, B1) Statistics\n');
fprintf('                        (Mean, Std)    (Mean, Std)\n');
fprintf('---------------------------------------------------\n');
for i = 1:length(n_levels)
    fprintf('n = %d                 (%.3f, %.3f)  (%.3f, %.3f)\n', ...
        n_levels(i), table6_style.Score_Mean(i), table6_style.Score_Std(i), ...
        table6_style.B1_Mean(i), table6_style.B1_Std(i));
end

writetable(table6_style, 'output/tab/multi_fixed_entry_table6_style.csv');

% Print key insights
fprintf('\n=== KEY INSIGHTS: Competition Effects ===\n');
fprintf('1. As competition increases (n increases):\n');
for i = 1:length(n_levels)
    fprintf('   n=%d: UP payment=%.4f, FP payment=%.4f, change=%.2f%%\n', ...
        summary_table.N_Entrants(i), summary_table.UP_Mean_Final_Payment(i), ...
        summary_table.FP_Mean_Final_Payment(i), summary_table.Final_Payment_Change_Pct(i));
end

fprintf('\n2. Competition intensity (UP vs FP payment ratio):\n');
for i = 1:length(n_levels)
    ratio = summary_table.FP_Mean_Final_Payment(i) / summary_table.UP_Mean_Final_Payment(i);
    fprintf('   n=%d: FP/UP ratio = %.3f\n', summary_table.N_Entrants(i), ratio);
end

fprintf('\nCompleted. Results saved in output/ and output/tab/ directories\n');

%% Helper function to run cf_run_fixed_entry silently and return results
function [results, bidder_stats] = cf_run_fixed_entry_silent(n_target)
    % Modified version of cf_run_fixed_entry that runs silently and returns results
    
    cfg = config_up_1d();
    
    % derive effective parameters
    if isfield(cfg,'scale') && isfield(cfg.scale,'useProjectScale') && cfg.scale.useProjectScale
        s = cfg.scale.s;
    else
        s = 1.0;
    end
    cfg.eff.theta0 = cfg.model.theta0 * s;
    cfg.eff.theta1 = cfg.model.theta1 * s;
    cfg.eff.alpha  = cfg.model.alpha  / s;
    cfg.eff.sigma_L = cfg.model.sigma_L; 
    cfg.eff.sigma_H = cfg.model.sigma_H;
    
    % Pre-build distributions
    riskTags = {'L','H'};
    FF_UP = struct(); FF_FP = struct();
    for r = 1:numel(riskTags)
        tag   = riskTags{r};
        sigma = cfg.eff.(['sigma_' tag]);
        
        % UP pseudo-cost via sampling + KDE
        c_emp = sample_c_via_sobol(cfg.kernel.M_c, sigma, cfg.eff.alpha, cfg);
        [Fh_u, fh_u, supp_u] = build_Ff_splines(c_emp, cfg.kernel.bw_mult);
        FF_UP.(tag) = struct('Fh',Fh_u,'fh',fh_u,'supp',supp_u,'sigma',sigma);
    
        % FP exact normal
        [Fh_f, fh_f, supp_f] = build_Ff_fp_normal(cfg, sigma);
        FF_FP.(tag) = struct('Fh',Fh_f,'fh',fh_f,'supp',supp_f,'sigma',sigma);
    end
    
    % Strategy caches
    SCacheUP.L = containers.Map('KeyType','double','ValueType','any');
    SCacheUP.H = containers.Map('KeyType','double','ValueType','any');
    SCacheFP.L = containers.Map('KeyType','double','ValueType','any');
    SCacheFP.H = containers.Map('KeyType','double','ValueType','any');
    
    get_s_handle = @(tag,n) get_or_build_s_handle(cfg, FF_UP.(tag).Fh, FF_UP.(tag).fh, ...
                                                  FF_UP.(tag).supp, n, tag, SCacheUP);
    get_b_handle = @(tag,n) get_or_build_bfp_handle(cfg, FF_FP.(tag).Fh, FF_FP.(tag).fh, ...
                                                    FF_FP.(tag).supp, n, tag, SCacheFP);
    
    % Run simulations
    nAuctions = cfg.sim.nAuctions;
    A = table('Size',[0 9], 'VariableTypes', ...
              {'double','string','double','double','double','double','double','double','double'}, ...
              'VariableNames', {'auction_id','risk','sigma','n', ...
                                'up_win_score','up_pay_final','fp_win_score','fp_pay_final','scale_s'});
    
    % Collect all bidder data for statistics
    all_scores = [];
    all_b1 = [];
    
    auction_id = 0;
    for t = 1:nAuctions
        % Draw risk
        if cfg.sim.useProb
            tag = (rand < cfg.sim.P) * 'H' + (rand >= cfg.sim.P) * 'L';
            if ~ischar(tag), tag = ternary(rand<cfg.sim.P,'H','L'); end
        else
            sch = cfg.sim.riskSchedule; tag = sch{mod(t-1,numel(sch))+1};
        end
        sigma = cfg.eff.(['sigma_' tag]);
    
        auction_id = auction_id + 1;
        [out, bidder_data] = cf_fixed_entry_one_enhanced(cfg, tag, sigma, n_target, ...
                                 FF_UP.(tag).Fh, FF_UP.(tag).fh, FF_UP.(tag).supp, ...
                                 FF_FP.(tag).Fh, FF_FP.(tag).fh, FF_FP.(tag).supp, ...
                                 @(n)get_s_handle(tag,n), @(n)get_b_handle(tag,n));
    
        A = [A; {auction_id, string(tag), sigma, n_target, ...
                 out.UP.win_score, out.UP.pay_final, ...
                 out.FP.win_score, out.FP.pay_final, s}];
        
        % Collect bidder data for statistics (rescaled)
        all_scores = [all_scores; bidder_data.scores / s];
        all_b1 = [all_b1; bidder_data.b1 / s];
    end
    
    results = A;
    
    % Calculate bidder statistics
    bidder_stats = struct();
    bidder_stats.score_mean = mean(all_scores);
    bidder_stats.score_std = std(all_scores);
    bidder_stats.b1_mean = mean(all_b1);
    bidder_stats.b1_std = std(all_b1);

    % Helper functions (same as in original)
    function s_handle = get_or_build_s_handle(cfg, Fh, fh, supp, n, tag, SCache)
        m = SCache.(tag);
        key = n + 0.0;
        if isKey(m, key)
            s_handle = m(key);
        else
            s_handle = build_s_handle(cfg, Fh, fh, supp, n);
            m(key) = s_handle;
            SCache.(tag) = m;
        end
    end

    function y = ternary(cond,a,b), if cond, y=a; else, y=b; end; end

    function b_handle = get_or_build_bfp_handle(cfg, Fh, fh, supp, n, tag, SCacheFP)
        m = SCacheFP.(tag);
        key = n + 0.0;
        if isKey(m,key)
            b_handle = m(key);
        else
            b_handle = build_bfp_handle(cfg, Fh, fh, supp, n);
            m(key) = b_handle;
            SCacheFP.(tag) = m;
        end
    end

    % Enhanced version of cf_fixed_entry_one that returns all bidder data
    function [out, bidder_data] = cf_fixed_entry_one_enhanced(cfg, riskTag, sigma, n, ...
        Fh_UP, fh_UP, supp_UP, ...
        Fh_FP, fh_FP, supp_FP, ...
        get_s_handle, get_b_handle)
        
        alp = cfg.eff.alpha;
        th0 = cfg.eff.theta0;
        th1 = cfg.eff.theta1;

        % Draw the SAME entrants' private types (e0,e2)
        [e0, e2] = draw_types_e0e2(n, cfg.type.sd_e0, cfg.type.sd_e2, cfg.type.rho_e);

        % ===================== UP branch =====================
        % pseudo-cost (Eq. 5, 1D): c_u = th0*e0 + th1 - (e2-1)^2 / (2*alpha*sigma^2)
        c_u = th0.*e0 + th1 - ((e2 - 1).^2) ./ (2*alp*sigma^2);

        % s(c;n) (ODE (7)) from cache
        s_handle = get_s_handle(n);
        s_u = s_handle(c_u);

        % inner skew (Eq. 4, 1D): b2* = th1 + (e2-1)/(alpha*sigma^2), trunc to [0,s]
        b2_u = th1 + (e2 - 1) ./ (alp * sigma^2);
        b2_u = max(0, min(b2_u, s_u));
        b0_u = s_u - b2_u;

        % winner & EXPECTED final payment (E_eps=0): p = b0 + b2*e2
        [win_s_u, wU] = min(s_u);
        pay_u = b0_u(wU) + b2_u(wU) * e2(wU);

        % ===================== FP branch =====================
        % FP type: c_fp = th0*e0 + th1*e2 + (alpha/2)*th1^2*sigma^2
        RP  = 0.5 * alp * (th1^2) * (sigma^2);
        c_fp = th0.*e0 + th1.*e2 + RP;

        % b(c;n) from cache
        b_handle = get_b_handle(n);
        b_fp = b_handle(c_fp);

        % winner & EXPECTED final payment (no adjustment): p = b
        [win_b_fp, wF] = min(b_fp);
        pay_fp = win_b_fp;

        % ===================== pack results =====================
        out = struct();
        out.n       = n;
        out.riskTag = riskTag;

        out.UP  = struct('win_score',win_s_u, 'pay_final',pay_u, ...
                        'winner_e0',e0(wU), 'winner_e2',e2(wU));
        out.FP  = struct('win_score',win_b_fp, 'pay_final',pay_fp, ...
                        'winner_e0',e0(wF),  'winner_e2',e2(wF));

        % Return all bidder data
        bidder_data = struct();
        bidder_data.scores = s_u;  % All UP scores
        bidder_data.b1 = b2_u;     % All UP b1 (non-lumpsum) values
    end
end