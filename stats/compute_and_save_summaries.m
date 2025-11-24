function S = compute_and_save_summaries(A, B, outdir)
    % COMPUTE_AND_SAVE_SUMMARIES
    % Summaries at auction- and bidder-level:
    %  - Entry probability delta (by risk & pooled)
    %  - Auction-level: n, overrun, win_score, pay_final (by risk & pooled)
    %  - Bidder-level (entrants only): mean/median/std of score, b0, b2 (by risk & pooled)
    %  - Empirical entry rate vs model delta
    %  - NEW: Homogenized (divide by scale_s) summaries for bids & auction-level money vars
    %
    % Inputs:
    %   A : auction-level table. REQUIRED columns:
    %       {'auction_id','risk','N_potential','delta','kbar','skipped','scale_s',
    %        'n','overrun','win_score','pay_final','epsilon','win_b0','win_b2','win_e0','win_e2'}
    %   B : bidder-level long table (N rows per auction; non-entrants have NaNs for bids).
    %       REQUIRED columns:
    %       {'auction_id','bidder_id','entered','k','e0','e2','c','score','b0','b2','winner'}
    %   outdir : folder to write CSVs (e.g., 'output')
    %
    % Output:
    %   S : struct of tables with all summaries (original-units + homogenized)
    
    if ~exist(outdir,'dir'); mkdir(outdir); end
    
    % Join risk & scale_s onto bidder panel to group by risk and to create homogenized vars
    B2 = innerjoin(B, A(:,{'auction_id','risk','N_potential','delta','kbar','skipped','scale_s'}), ...
                   'Keys','auction_id');
    
    % Keep non-skipped auctions
    A_valid  = A(~A.skipped,:);
    B2_valid = B2(~B2.skipped,:);
    
    %% 1) Entry probability delta (auction-level)
    S.entry_by_risk = groupsummary(A_valid, 'risk', {'mean','median','std'}, 'delta');
    
    S.entry_all = table;
    S.entry_all.risk         = "ALL";
    S.entry_all.GroupCount   = height(A_valid);
    S.entry_all.mean_delta   = mean(A_valid.delta);
    S.entry_all.median_delta = median(A_valid.delta);
    S.entry_all.std_delta    = std(A_valid.delta);
    
    %% 2) Auction-level stats (original units)
    varsA = {'n','overrun','win_score','pay_final'};
    S.auc_by_risk = groupsummary(A_valid, 'risk', {'mean','median','std'}, varsA);
    
    S.auc_all = table;
    S.auc_all.risk           = "ALL";
    S.auc_all.GroupCount     = height(A_valid);
    S.auc_all.mean_n         = mean(A_valid.n);
    S.auc_all.median_n       = median(A_valid.n);
    S.auc_all.std_n          = std(A_valid.n);
    S.auc_all.mean_overrun   = mean(A_valid.overrun,'omitnan');
    S.auc_all.mean_win_score = mean(A_valid.win_score,'omitnan');
    S.auc_all.mean_pay_final = mean(A_valid.pay_final,'omitnan');
    
    %% 3) Entrants' bids (bidder-level; entrants only, original units)
    ENT = B2_valid.entered;
    varsBid = {'score','b0','b2'};
    S.bids_by_risk = groupsummary(B2_valid(ENT,:), 'risk', {'mean','median','std'}, varsBid);
    
    S.bids_all = table;
    S.bids_all.risk           = "ALL";
    S.bids_all.GroupCount     = sum(ENT);
    S.bids_all.mean_score     = mean(B2_valid.score(ENT),'omitnan');
    S.bids_all.median_score   = median(B2_valid.score(ENT),'omitnan');
    S.bids_all.std_score      = std(B2_valid.score(ENT),'omitnan');
    S.bids_all.mean_b0        = mean(B2_valid.b0(ENT),'omitnan');
    S.bids_all.mean_b2        = mean(B2_valid.b2(ENT),'omitnan');
    
    %% 3bis) NEW: Homogenized (divide by scale_s) entrants' bids & auction-level money vars
    % Create homogenized (÷s) columns
    B2_valid.score_h = B2_valid.score ./ B2_valid.scale_s;
    B2_valid.b0_h    = B2_valid.b0    ./ B2_valid.scale_s;
    B2_valid.b2_h    = B2_valid.b2    ./ B2_valid.scale_s;
    
    A_valid.win_score_h = A_valid.win_score ./ A_valid.scale_s;
    A_valid.pay_final_h = A_valid.pay_final ./ A_valid.scale_s;
    A_valid.overrun_h   = A_valid.overrun   ./ A_valid.scale_s;
    
    % Entrants' bids (homogenized) by risk & pooled
    S.bids_h_by_risk = groupsummary(B2_valid(ENT,:), 'risk', {'mean','median','std'}, {'score_h','b0_h','b2_h'});
    
    S.bids_h_all = table;
    S.bids_h_all.risk             = "ALL";
    S.bids_h_all.GroupCount       = sum(ENT);
    S.bids_h_all.mean_score_h     = mean(B2_valid.score_h(ENT),'omitnan');
    S.bids_h_all.median_score_h   = median(B2_valid.score_h(ENT),'omitnan');
    S.bids_h_all.std_score_h      = std(B2_valid.score_h(ENT),'omitnan');
    S.bids_h_all.mean_b0_h        = mean(B2_valid.b0_h(ENT),'omitnan');
    S.bids_h_all.mean_b2_h        = mean(B2_valid.b2_h(ENT),'omitnan');
    
    % Auction-level (homogenized) by risk & pooled (optional but useful)
    S.auc_h_by_risk = groupsummary(A_valid, 'risk', {'mean','median','std'}, {'win_score_h','pay_final_h','overrun_h'});
    
    S.auc_h_all = table;
    S.auc_h_all.risk               = "ALL";
    S.auc_h_all.GroupCount         = height(A_valid);
    S.auc_h_all.mean_win_score_h   = mean(A_valid.win_score_h,'omitnan');
    S.auc_h_all.median_win_score_h = median(A_valid.win_score_h,'omitnan');
    S.auc_h_all.std_win_score_h    = std(A_valid.win_score_h,'omitnan');
    S.auc_h_all.mean_pay_final_h   = mean(A_valid.pay_final_h,'omitnan');
    S.auc_h_all.mean_overrun_h     = mean(A_valid.overrun_h,'omitnan');
    
    % (Strongly recommended) Condition by realized n to compare with Table 6
    Tn = innerjoin(B2_valid(ENT,:), A_valid(:,{'auction_id','n'}), 'Keys','auction_id');
    S.bids_h_by_n = groupsummary(Tn, 'n', {'mean','median','std'}, {'score_h','b0_h','b2_h'});
    
    %% 4) Empirical entry rate vs model delta (sanity)
    emp_rate = groupsummary(B2_valid, 'auction_id', 'mean', 'entered');
    T_emp = innerjoin(emp_rate(:,{'auction_id','mean_entered'}), A_valid(:,{'auction_id','risk','delta'}), 'Keys','auction_id');
    T_emp.Properties.VariableNames{'mean_entered'} = 'emp_entry_rate';
    S.emp_by_risk = groupsummary(T_emp, 'risk', {'mean','std'}, {'emp_entry_rate','delta'});
    
    S.emp_all = table;
    S.emp_all.risk                  = "ALL";
    S.emp_all.GroupCount            = height(T_emp);
    S.emp_all.mean_emp_entry_rate   = mean(T_emp.emp_entry_rate);
    S.emp_all.mean_delta            = mean(T_emp.delta);
    
    %% 5) Save CSVs (original + homogenized)
    % Original units
    writetable(S.entry_by_risk, fullfile(outdir,'summary_entry_by_risk.csv'));
    writetable(S.auc_by_risk,   fullfile(outdir,'summary_auc_by_risk.csv'));
    writetable(S.bids_by_risk,  fullfile(outdir,'summary_bids_by_risk.csv'));
    
    T = S.entry_all; T.Properties.VariableNames = {'risk','GroupCount','mean_delta','median_delta','std_delta'};
    writetable(T, fullfile(outdir,'summary_entry_pooled.csv'));
    
    T = S.auc_all; T.Properties.VariableNames = {'risk','GroupCount','mean_n','median_n','std_n','mean_overrun','mean_win_score','mean_pay_final'};
    writetable(T, fullfile(outdir,'summary_auc_pooled.csv'));
    
    T = S.bids_all; T.Properties.VariableNames = {'risk','GroupCount','mean_score','median_score','std_score','mean_b0','mean_b2'};
    writetable(T, fullfile(outdir,'summary_bids_pooled.csv'));
    
    % Homogenized (÷s)
    writetable(S.bids_h_by_risk,  fullfile(outdir,'summary_bids_h_by_risk.csv'));
    T = S.bids_h_all; T.Properties.VariableNames = {'risk','GroupCount','mean_score_h','median_score_h','std_score_h','mean_b0_h','mean_b2_h'};
    writetable(T, fullfile(outdir,'summary_bids_h_pooled.csv'));
    
    writetable(S.auc_h_by_risk,   fullfile(outdir,'summary_auc_h_by_risk.csv'));
    T = S.auc_h_all; T.Properties.VariableNames = {'risk','GroupCount','mean_win_score_h','median_win_score_h','std_win_score_h','mean_pay_final_h','mean_overrun_h'};
    writetable(T, fullfile(outdir,'summary_auc_h_pooled.csv'));
    
    % By n (homogenized) to compare with Table 6
    writetable(S.bids_h_by_n,     fullfile(outdir,'summary_bids_h_by_n.csv'));
    
    % Also print a compact header to console
    disp('== Entry probability (by risk) ==');            disp(S.entry_by_risk);
    disp('== Entry probability (pooled) ==');             disp(S.entry_all);
    disp('== Auction-level (by risk, original) ==');      disp(S.auc_by_risk);
    disp('== Auction-level (pooled, original) ==');       disp(S.auc_all);
    disp('== Entrants'' bids (by risk, original) ==');    disp(S.bids_by_risk);
    disp('== Entrants'' bids (pooled, original) ==');     disp(S.bids_all);
    
    disp('== Entrants'' bids (homog., by risk) ==');      disp(S.bids_h_by_risk);
    disp('== Entrants'' bids (homog., pooled) ==');       disp(S.bids_h_all);
    disp('== Auction-level (homog., by risk) ==');        disp(S.auc_h_by_risk);
    disp('== Auction-level (homog., pooled) ==');         disp(S.auc_h_all);
    
    disp('== Empirical entry vs delta (by risk) ==');     disp(S.emp_by_risk);
    disp('== Empirical entry vs delta (pooled) ==');      disp(S.emp_all);
end
    