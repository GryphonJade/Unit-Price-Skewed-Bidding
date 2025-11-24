%% ===================== UP (1D) DGP — Main Runner ===================== %%
clear; clc; rng(1);

% --- add paths
addpath(genpath('config')); addpath(genpath('utils'));
addpath(genpath('strategy')); addpath(genpath('sim'));
addpath(genpath('stats')); 
% --- load config (hyper-parameters & switches)
cfg = config_up_1d();


% derive effective (scaled) parameters, per Eq.(18)
if cfg.scale.useProjectScale
    scale_s = cfg.scale.s;
else
    scale_s = 1.0;
end
cfg.eff.theta0 = cfg.model.theta0 * scale_s;
cfg.eff.theta1 = cfg.model.theta1 * scale_s;
cfg.eff.alpha  = cfg.model.alpha  / scale_s;    % CARA rescaling



% =================== 0) Pre-build F_u(c) per risk state ===================
% DGP tie-in:
%  - Entry fixed point (Eq. (8)) needs E[u(s(c;n)-c)] under the UNCONDITIONAL type law F_u.
%  - s(c;n) ODE (Eq. (7)) also uses the same F, f.
riskTags = {'L','H'};              % we prepare both; sampling uses P or schedule
FFcache  = struct();               % holds Fh, fh, supp, and raw c draws per risk
for r = 1:numel(riskTags)
    tag   = riskTags{r};
    sigma = cfg.model.(['sigma_' tag]);
    c_emp = sample_c_via_sobol(cfg.kernel.M_c, sigma, cfg.eff.alpha, cfg);     % Eq.(5)
    [Fh, fh, supp] = build_Ff_splines(c_emp, cfg.kernel.bw_mult);                % smoothed F,f
    FFcache.(tag) = struct('Fh',Fh,'fh',fh,'supp',supp,'sigma',sigma,'cEmp',c_emp);
end

% =================== 1) Cache entry fixed point per risk ===================
% DGP tie-in:
%  - δ solves h(δ)=Φ((CE(δ)-μ_k)/σ_k) - δ = 0 where CE(δ)=u^{-1}(U(δ)), u(x)=-exp(-αx).
%  - U(δ)= Σ_n Pr(n|δ)*(E u(s(c;n)-c)); the EU depends on risk via F_u and on α, θ.
% Optimization:
%  - Compute once per risk state and re-use across auctions (same config).
EntryCache = struct();
if cfg.cache.precomputeEntry
    for r = 1:numel(riskTags)
        tag = riskTags{r};
        pack = FFcache.(tag);
        [delta,kbar,CE,U,traceE] = compute_entry_probability_bisect( ...
            cfg, pack.sigma, pack.Fh, pack.fh, pack.supp);
        EntryCache.(tag) = struct('delta',delta,'kbar',kbar,'CE',CE,'U',U,'trace',traceE);
        if cfg.trace.verbose
            fprintf('[ENTRY] risk=%s  delta=%.4f  kbar=%.5f\n', tag, delta, kbar);
        end
    end
end

% =================== 2) Lazy cache s(c;n) per (risk, n) ===================
% DGP tie-in:
%  - s(c;n) is determined by ODE (7) given F,f (risk) and n; it’s re-used many times.
SCache = struct();                 % maps riskTag -> containers.Map(n -> s_handle)
SCache.L = containers.Map('KeyType','double','ValueType','any');
SCache.H = containers.Map('KeyType','double','ValueType','any');

% Provide a closure that returns a cached s_handle, otherwise builds and stores it.
get_s_handle = @(tag,n) get_or_build_s_handle(cfg, FFcache.(tag).Fh, FFcache.(tag).fh, ...
                                              FFcache.(tag).supp, n, tag, SCache);

% =================== 3) Output containers (auction- & bidder-level) =======
% Design: (i) auction-level table; (ii) bidder-level long table (panel).
A_cols = {'auction_id','risk','sigma','N_potential','delta','kbar','n',...
          'win_score','win_b0','win_b2','win_e0','win_e2','epsilon','pay_final','overrun',...
          'skipped','scale_s'};                     
A = table('Size',[0 numel(A_cols)], 'VariableTypes', ...
          {'double','string','double','double','double','double','double',...
           'double','double','double','double','double','double','double','double',...
           'logical','double'}, ...                
          'VariableNames',A_cols);


B_cols = {'auction_id','bidder_id','entered','k','e0','e2','c','score','b0','b2','winner'};
B = table('Size',[0 numel(B_cols)], 'VariableTypes', ...
         {'double','double','logical','double','double','double','double','double','double','double','logical'}, ...
         'VariableNames',B_cols);
           


% Save paths
if ~exist('output','dir'); mkdir('output'); end
save_mat_path = fullfile('output','dgp_results.mat');
save_csv_auct = fullfile('output','auctions.csv');
save_csv_bids = fullfile('output','bids.csv');

% =================== 4) Simulation loop over auctions =====================
t0 = tic;
auction_id = 0;
skipped_cnt = 0;

for t = 1:cfg.sim.nAuctions
    % ----- choose risk state per P or schedule -----
    if cfg.sim.useProb
        tag = (rand < cfg.sim.P) * 'H' + (rand >= cfg.sim.P) * 'L';
        if ~ischar(tag); tag = rand < cfg.sim.P; tag = ternary(tag,'H','L'); end % robust cast
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
        entry = []; entry.precomputed = false; % simulate_one_auction will compute
    end

    % ----- simulate ONE auction -----
    % We pass a handle to obtain cached s(c;n): get_s_handle(tag,n)
    % and an (optional) precomputed entry package.
    auction_id = auction_id + 1;
    out = simulate_one_auction_full(cfg, tag, sigma, pack.Fh, pack.fh, pack.supp, ...
                                    @(n)get_s_handle(tag,n), entry);

    if out.skipped
        skipped_cnt = skipped_cnt + 1;
        % Record a minimal auction row so you can track skips in estimation if needed
        A = [A; {auction_id, string(tag), sigma, cfg.sim.N_potential, out.delta, out.kbar, ...
                 out.n, NaN,NaN,NaN,NaN,NaN, NaN, NaN, NaN, true, scale_s}];
        continue;
    end

    % ----- append auction-level row -----
    A = [A; {auction_id, string(tag), sigma, cfg.sim.N_potential, out.delta, out.kbar, ...
             out.n, out.win_s, out.win_b0, out.win_b2, out.e0, out.e2, ...
             out.epsilon, out.pay_final, out.overrun, false, scale_s}];

    % ----- append bidder-level rows (panel) -----
    % ----- append bidder-level rows (panel) -----
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
        fprintf('[%4d/%4d] risk=%s  n=%d  delta=%.3f  overrun=%.4f\n', ...
            t, cfg.sim.nAuctions, tag, out.n, out.delta, out.overrun);
    end
end

elapsed = toc(t0);
fprintf('Simulation finished in %.2fs. Avg n = %.2f (skipped=%d)\n', ...
        elapsed, mean(A.n(~A.skipped)), skipped_cnt);

% =================== 5) Persist all data for estimation ====================
% Save both MAT (lossless struct/tables) and CSV (interoperable) formats
save(save_mat_path, 'A', 'B', 'cfg', 'FFcache', 'EntryCache', 'SCache', '-v7.3');
writetable(A, save_csv_auct);
writetable(B, save_csv_bids);

% ===== sumaries =====
statsOut = compute_and_save_summaries(A, B, 'output');  % returns a struct of tables



    % -------- tiny helpers --------
    function y = ternary(cond, a, b), if cond, y=a; else, y=b; end
    end
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
    