%% ===================== FP (1D) DGP — Main Runner ===================== %%
clear; clc; rng(1);

% --- add paths
addpath(genpath('config')); addpath(genpath('utils'));
addpath(genpath('strategy')); addpath(genpath('sim'));
addpath(genpath('stats'));

% --- load config (hyper-parameters & switches)
cfg = config_up_1d();

% --- derive effective (scaled) parameters, per Eq.(18) ---
if isfield(cfg,'scale') && isfield(cfg.scale,'useProjectScale') && cfg.scale.useProjectScale
    scale_s = cfg.scale.s;
else
    scale_s = 1.0;
end
cfg.eff.theta0 = cfg.model.theta0 * scale_s;
cfg.eff.theta1 = cfg.model.theta1 * scale_s;
cfg.eff.alpha  = cfg.model.alpha  / scale_s;   % CARA rescaling per Eq.(18)

% =================== 0) Pre-build F_fp(c) per risk state ===================
% FP type: c^FP = theta0*e0 + theta1*e2 + (alpha/2)*theta1^2*sigma^2  (exact Normal)
riskTags = {'L','H'};
FFcache  = struct();
for r = 1:numel(riskTags)
    tag   = riskTags{r};
    sigma = cfg.model.(['sigma_' tag]);     % sigma_t(X) kept constant per Eq.(18)
    [Fh, fh, supp] = build_Ff_fp_normal(cfg, sigma);   % uses cfg.eff.{theta,alpha}
    FFcache.(tag) = struct('Fh',Fh,'fh',fh,'supp',supp,'sigma',sigma);
end

% =================== 1) Cache entry fixed point per risk (FP) ===================
EntryCache = struct();
if cfg.cache.precomputeEntry
    for r = 1:numel(riskTags)
        tag = riskTags{r};
        pack = FFcache.(tag);
        % IMPORTANT: this FP solver must internally use expected_utility_integral_weighted
        [delta,kbar,CE,U,traceE] = compute_entry_probability_bisect_fp( ...
            cfg, pack.sigma, pack.Fh, pack.fh, pack.supp);
        EntryCache.(tag) = struct('delta',delta,'kbar',kbar,'CE',CE,'U',U,'trace',traceE);
        if cfg.trace.verbose
            fprintf('[ENTRY-FP] risk=%s  delta=%.4f  kbar=%.5f\n', tag, delta, kbar);
        end
    end
end

% =================== 2) Lazy cache b(c;n) per (risk, n) ===================
SCache = struct();                 % maps riskTag -> containers.Map(n -> b_handle)
SCache.L = containers.Map('KeyType','double','ValueType','any');
SCache.H = containers.Map('KeyType','double','ValueType','any');

% Provide a closure returning cached b_handle, or build-and-store it.
get_b_handle = @(tag,n) get_or_build_bfp_handle(cfg, FFcache.(tag).Fh, FFcache.(tag).fh, ...
                                                FFcache.(tag).supp, n, tag, SCache);

% =================== 3) Output containers (auction- & bidder-level) =======
% Keep exactly the same columns as your UP runner (so estimation code can reuse)
A_cols = {'auction_id','risk','sigma','N_potential','delta','kbar','n', ...
          'win_score','win_b0','win_b2','win_e0','win_e2','epsilon','pay_final','overrun', ...
          'skipped','scale_s'};
A = table('Size',[0 numel(A_cols)], 'VariableTypes', ...
          {'double','string','double','double','double','double','double', ...
           'double','double','double','double','double','double','double','double', ...
           'logical','double'}, ...
          'VariableNames',A_cols);

B_cols = {'auction_id','bidder_id','entered','k','e0','e2','c','score','b0','b2','winner'};
B = table('Size',[0 numel(B_cols)], 'VariableTypes', ...
          {'double','double','logical','double','double','double','double','double','double','double','logical'}, ...
          'VariableNames',B_cols);

% Save paths
if ~exist('output','dir'); mkdir('output'); end
save_mat_path = fullfile('output','dgp_results_fp.mat');
save_csv_auct = fullfile('output','auctions_fp.csv');
save_csv_bids = fullfile('output','bids_fp.csv');

% =================== 4) Simulation loop over auctions =====================
t0 = tic;
auction_id = 0;
skipped_cnt = 0;

for t = 1:cfg.sim.nAuctions
    % ----- choose risk state per P or schedule -----
    if cfg.sim.useProb
        tag = (rand < cfg.sim.P) * 'H' + (rand >= cfg.sim.P) * 'L';
        if ~ischar(tag); tag = rand < cfg.sim.P; tag = ternary(tag,'H','L'); end
    else
        sch = cfg.sim.riskSchedule;
        tag = sch{mod(t-1, numel(sch))+1};
    end
    pack  = FFcache.(tag);
    sigma = pack.sigma;

    % ----- choose entry params (precomputed vs compute-on-the-fly) -----
    if cfg.cache.precomputeEntry
        entry.delta = EntryCache.(tag).delta;
        entry.kbar  = EntryCache.(tag).kbar;
        entry.CE    = EntryCache.(tag).CE;
        entry.U     = EntryCache.(tag).U;
        entry.precomputed = true;
    else
        entry = []; entry.precomputed = false; % simulate_one_auction_full_fp will compute
    end

    % ----- simulate ONE FP auction -----
    auction_id = auction_id + 1;
    out = simulate_one_auction_full_fp(cfg, tag, sigma, pack.Fh, pack.fh, pack.supp, ...
                                       @(n)get_b_handle(tag,n), entry);

    if out.skipped
        skipped_cnt = skipped_cnt + 1;
        % Record a minimal auction row for skipped case
        A = [A; {auction_id, string(tag), sigma, cfg.sim.N_potential, out.delta, out.kbar, ...
                 out.n, NaN,NaN,NaN,NaN,NaN, NaN, NaN, NaN, true, scale_s}];
        continue;
    end

    % ----- append auction-level row -----
    A = [A; {auction_id, string(tag), sigma, cfg.sim.N_potential, out.delta, out.kbar, ...
             out.n, out.win_s, out.win_b0, out.win_b2, out.e0, out.e2, ...
             out.epsilon, out.pay_final, out.overrun, false, scale_s}];

    % ----- append bidder-level rows (panel; N rows, non-entrants = NaN bids) -----
    N = cfg.sim.N_potential;
    bidder_id = (1:N)';

    Bt = table( repmat(auction_id,N,1), ...
                bidder_id, ...
                out.panel.entered(:), ...
                out.panel.k(:), ...
                out.panel.e0(:), out.panel.e2(:), out.panel.c(:), ...
                out.panel.s(:), out.panel.b0(:), out.panel.b2(:), ...
                out.panel.winner(:), ...
                'VariableNames', {'auction_id','bidder_id','entered','k','e0','e2','c','score','b0','b2','winner'} );
    B = [B; Bt];

    if cfg.trace.verbose && (mod(t, cfg.trace.printEvery)==0)
        fprintf('[FP %4d/%4d] risk=%s  n=%d  delta=%.3f\n', ...
            t, cfg.sim.nAuctions, tag, out.n, out.delta);
    end
end

elapsed = toc(t0);
fprintf('FP Simulation finished in %.2fs. Avg n = %.2f (skipped=%d)\n', ...
        elapsed, mean(A.n(~A.skipped)), skipped_cnt);

% =================== 5) Persist all data for estimation ====================
save(save_mat_path, 'A', 'B', 'cfg', 'FFcache', 'EntryCache', 'SCache', '-v7.3');
writetable(A, save_csv_auct);
writetable(B, save_csv_bids);

% ===== compute rich summaries + optional CSVs (by risk & pooled) =====
statsOut = compute_and_save_summaries(A, B, 'output');  %#ok<NASGU>

% ===== run consistency checks (assertions); set 'strict',true to error on fail =====
validate_dgp(A, B, cfg, 'strict', false);

% -------- tiny helpers --------
function y = ternary(cond, a, b), if cond, y=a; else, y=b; end
end

function b_handle = get_or_build_bfp_handle(cfg, Fh, fh, supp, n, tag, SCache)
% Return cached b(c;n) for (risk=tag, n); otherwise build & store.
m = SCache.(tag);
key = n + 0.0;
if isKey(m, key)
    b_handle = m(key);
else
    b_handle = build_bfp_handle(cfg, Fh, fh, supp, n);  % analytic FP bid curve
    m(key)   = b_handle;
    SCache.(tag) = m;
end
end
