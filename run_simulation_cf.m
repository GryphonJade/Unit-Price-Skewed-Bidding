%% ======================= UP vs FP Counterfactual Runner ======================= %%
clear; clc; rng(1);

addpath(genpath('config')); addpath(genpath('utils'));
addpath(genpath('strategy')); addpath(genpath('sim'));
addpath(genpath('stats')); addpath(genpath('helpers'));

cfg = config_up_1d();

% ---- derive effective (scaled) parameters per Eq.(18) ----
if isfield(cfg,'scale') && isfield(cfg.scale,'useProjectScale') && cfg.scale.useProjectScale
    s = cfg.scale.s;
else
    s = 1.0;
end
cfg.eff.theta0  = cfg.model.theta0 * s;
cfg.eff.theta1  = cfg.model.theta1 * s;
cfg.eff.alpha   = cfg.model.alpha  / s;
cfg.eff.sigma_L = cfg.model.sigma_L; 
cfg.eff.sigma_H = cfg.model.sigma_H;

% =================== Pre-build F_u(c) for UP and F_fp(c) for FP ===================
riskTags = {'L','H'};
FF_UP  = struct();    % UP: pseudo-cost distribution via sample (non-Gaussian)
FF_FP  = struct();    % FP: exact normal for c^FP

for r = 1:numel(riskTags)
    tag   = riskTags{r};
    sigma = cfg.eff.(['sigma_' tag]);    % same as model; kept for clarity

    % UP (nonlinear pseudo-cost c_u): sample + KDE (you already have)
    c_emp = sample_c_via_sobol(cfg.kernel.M_c, sigma, cfg.eff.alpha, cfg);   % ensure uses eff.*
    [Fh_u, fh_u, supp_u] = build_Ff_splines(c_emp, cfg.kernel.bw_mult);
    FF_UP.(tag) = struct('Fh',Fh_u,'fh',fh_u,'supp',supp_u,'sigma',sigma);

    % FP (linear + normal risk premium): exact normal
    [Fh_f, fh_f, supp_f] = build_Ff_fp_normal(cfg, sigma);
    FF_FP.(tag) = struct('Fh',Fh_f,'fh',fh_f,'supp',supp_f,'sigma',sigma);
end

% =================== Precompute entry δ and cache strategies =====================
Entry_UP = struct(); Entry_FP = struct();
SCacheUP.L = containers.Map('KeyType','double','ValueType','any');
SCacheUP.H = containers.Map('KeyType','double','ValueType','any');
SCacheFP.L = containers.Map('KeyType','double','ValueType','any');
SCacheFP.H = containers.Map('KeyType','double','ValueType','any');

get_s_handle = @(tag,n) get_or_build_s_handle(cfg, FF_UP.(tag).Fh, FF_UP.(tag).fh, FF_UP.(tag).supp, n, tag, SCacheUP);
get_b_handle = @(tag,n) get_or_build_bfp_handle(cfg, FF_FP.(tag).Fh, FF_FP.(tag).fh, FF_FP.(tag).supp, n, tag, SCacheFP);

if cfg.cache.precomputeEntry
    for r = 1:numel(riskTags)
        tag = riskTags{r};
        % UP
        pack = FF_UP.(tag);
        [delta,kbar,CE,U,traceE] = compute_entry_probability_bisect(cfg, pack.sigma, pack.Fh, pack.fh, pack.supp);
        Entry_UP.(tag) = struct('delta',delta,'kbar',kbar,'CE',CE,'U',U,'trace',traceE);
        % FP
        pack = FF_FP.(tag);
        [delta,kbar,CE,U,traceE] = compute_entry_probability_bisect_fp(cfg, pack.sigma, pack.Fh, pack.fh, pack.supp);
        Entry_FP.(tag) = struct('delta',delta,'kbar',kbar,'CE',CE,'U',U,'trace',traceE);
    end
end

% =================== Output tables (add mechanism column) ========================
A_cols = {'auction_id','mechanism','risk','sigma','N_potential','delta','kbar','n',...
          'win_score','win_b0','win_b2','win_e0','win_e2','epsilon','pay_final','overrun','skipped','scale_s'};
B_cols = {'auction_id','mechanism','bidder_id','entered','k','e0','e2','c','score','b0','b2','winner'};

A_cf = table('Size',[0 numel(A_cols)], 'VariableTypes', ...
         {'double','string','string','double','double','double','double','double',...
          'double','double','double','double','double','double','double','double','logical','double'}, ...
          'VariableNames',A_cols);

B_cf = table('Size',[0 numel(B_cols)], 'VariableTypes', ...
         {'double','string','double','logical','double','double','double','double','double','double','double','logical'}, ...
          'VariableNames',B_cols);

% paths
if ~exist('output','dir'); mkdir('output'); end
save_mat_path = fullfile('output','dgp_results_cf.mat');

% =================== Loop over auctions: simulate UP & FP ========================
t0 = tic; auction_id = 0;
for t = 1:cfg.sim.nAuctions
    % risk state
    if cfg.sim.useProb
        tag = (rand < cfg.sim.P)*'H' + (rand >= cfg.sim.P)*'L';
        if ~ischar(tag); tag = rand < cfg.sim.P; tag = ternary(tag,'H','L'); end
    else
        sch = cfg.sim.riskSchedule; tag = sch{mod(t-1, numel(sch))+1};
    end
    auction_id = auction_id + 1;
    sigma = FF_UP.(tag).sigma;  s_meas = cfg.scale.s;

    % ---------- UP ----------
    if cfg.cache.precomputeEntry
        entry_up = Entry_UP.(tag); entry_up.precomputed = true;
    else
        entry_up = []; entry_up.precomputed=false;
    end
    outUP = simulate_one_auction_full(cfg, tag, sigma, FF_UP.(tag).Fh, FF_UP.(tag).fh, FF_UP.(tag).supp, ...
                                      @(n)get_s_handle(tag,n), entry_up);
    % write UP
    A_cf = [A_cf; {auction_id,"UP",string(tag),sigma,cfg.sim.N_potential, outUP.delta,outUP.kbar,outUP.n, ...
                   outUP.win_s,outUP.win_b0,outUP.win_b2,outUP.e0,outUP.e2, outUP.epsilon, outUP.pay_final,outUP.overrun, outUP.skipped, s_meas}];
    N = cfg.sim.N_potential; bidder_id = (1:N)';
    Bt = table( repmat(auction_id,N,1), repmat("UP",N,1), ...
                bidder_id, outUP.panel.entered(:), outUP.panel.k(:), ...
                outUP.panel.e0(:), outUP.panel.e2(:), outUP.panel.c(:), ...
                outUP.panel.s(:), outUP.panel.b0(:), outUP.panel.b2(:), outUP.panel.winner(:), ...
                'VariableNames', B_cols );
    B_cf = [B_cf; Bt];

    % ---------- FP ----------
    if cfg.cache.precomputeEntry
        entry_fp = Entry_FP.(tag); entry_fp.precomputed = true;
    else
        entry_fp = []; entry_fp.precomputed=false;
    end
    outFP = simulate_one_auction_full_fp(cfg, tag, sigma, FF_FP.(tag).Fh, FF_FP.(tag).fh, FF_FP.(tag).supp, ...
                                         @(n)get_b_handle(tag,n), entry_fp);
    A_cf = [A_cf; {auction_id,"FP",string(tag),sigma,cfg.sim.N_potential, outFP.delta,outFP.kbar,outFP.n, ...
                   outFP.win_s,outFP.win_b0,outFP.win_b2,outFP.e0,outFP.e2, outFP.epsilon, outFP.pay_final,outFP.overrun, outFP.skipped, s_meas}];
    Bt = table( repmat(auction_id,N,1), repmat("FP",N,1), ...
                bidder_id, outFP.panel.entered(:), outFP.panel.k(:), ...
                outFP.panel.e0(:), outFP.panel.e2(:), outFP.panel.c(:), ...
                outFP.panel.s(:), outFP.panel.b0(:), outFP.panel.b2(:), outFP.panel.winner(:), ...
                'VariableNames', B_cols );
    B_cf = [B_cf; Bt];
end

elapsed = toc(t0);
fprintf('CF simulation completed: %d auctions in %.2fs\n', cfg.sim.nAuctions, elapsed);

% Create output directories
if ~exist('output/fig','dir'); mkdir('output/fig'); end
if ~exist('output/tab','dir'); mkdir('output/tab'); end

% Save complete simulation data
save(save_mat_path, 'A_cf', 'B_cf', 'cfg', 'FF_UP', 'FF_FP', 'Entry_UP', 'Entry_FP', 'SCacheUP', 'SCacheFP', '-v7.3');
writetable(A_cf, fullfile('output','auctions_cf.csv'));
writetable(B_cf, fullfile('output','bids_cf.csv'));

% Generate run number for this session
run_number = get_next_run_number();
fprintf('Generating visualizations and summary tables (Run #%02d)...\n', run_number);

% Generate visualizations with run number
generate_cf_visualizations(A_cf, B_cf, run_number);

% Generate summary tables with run number
generate_cf_summary_tables(A_cf, B_cf, run_number);

fprintf('Analysis complete (Run #%02d). Check output/fig/ and output/tab/ folders.\n', run_number);

% -------- helper functions --------
function s_handle = get_or_build_s_handle(cfg, Fh, fh, supp, n, tag, SCache)
    % Return cached s(c;n) for (risk=tag, n); otherwise build & store.
    m = SCache.(tag);
    key = n + 0.0;
    if isKey(m, key)
        s_handle = m(key);
    else
        s_handle = build_s_handle(cfg, Fh, fh, supp, n);  % ODE (7) solver
        m(key)   = s_handle;
        SCache.(tag) = m;
    end
end

function y = ternary(cond,a,b), if cond, y=a; else, y=b; end; end

% strategy/get_or_build_bfp_handle.m
function b_handle = get_or_build_bfp_handle(cfg, Fh, fh, supp, n, tag, SCacheFP)
    m = SCacheFP.(tag);
    key = n + 0.0;
    if isKey(m,key)
        b_handle = m(key);
    else
        b_handle = build_bfp_handle(cfg, Fh, fh, supp, n);
        m(key)   = b_handle;
        SCacheFP.(tag) = m;
    end
    end

function run_number = get_next_run_number()
    % Get the next sequential run number for output files
    run_counter_file = 'output/.run_counter';
    
    if exist(run_counter_file, 'file')
        fid = fopen(run_counter_file, 'r');
        run_number = fscanf(fid, '%d');
        fclose(fid);
        if isempty(run_number)
            run_number = 1;
        else
            run_number = run_number + 1;
        end
    else
        run_number = 1;
    end
    
    % Save updated counter
    fid = fopen(run_counter_file, 'w');
    fprintf(fid, '%d', run_number);
    fclose(fid);
end
