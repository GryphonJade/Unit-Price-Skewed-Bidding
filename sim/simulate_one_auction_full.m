function out = simulate_one_auction_full(cfg, riskTag, sigma, Fh, fh, supp, get_s_handle, entry)
    % =============== Simulate one UP (1D) auction with endogenous entry =============== %
    % Inputs:
    %   cfg      : config struct
    %   riskTag  : 'L' or 'H' (only for labeling in trace)
    %   sigma    : project risk (std of demand shock)
    %   Fh, fh   : handles for CDF and PDF of pseudo-cost c (built on unconditional F_u)
    %   supp     : [L,R] support used for building s(c;n)
    %
    % Correct timeline per the paper:
    %   (1) Draw entry costs k only (types are NOT known before entry).
    %   (2) Solve endogenous entry δ via fixed point on Eq. (8) (we use rigorous bisection).
    %   (3) Realize entrants by the threshold rule: enter iff k <= k_bar (k_bar = CE(δ)).
    %   (4) Only for entrants, draw private types (e0,e2); compute pseudo-costs c (Eq. (5), 1D).
    %   (5) Build s(c;n) via ODE (7) and evaluate s_i = s(c_i;n).
    %   (6) Inner skew (Eq. (4), 1D): b2*, then b0* = s - b2*.
    %   (7) Winner, implementation shock ε ~ N(0,σ^2), final payment & overrun.
    
    % shorthand
    N   = cfg.sim.N_potential;
    alp = cfg.eff.alpha;
    th0 = cfg.eff.theta0;
    th1 = cfg.eff.theta1;
    
    % (1) Draw entry costs only (types are realized AFTER entry).
    k = draw_entry_costs_norm(N, cfg.entry.mu_k, cfg.entry.sd_k);
    
    % (2) Entry fixed point: precomputed vs compute now
    if isstruct(entry) && isfield(entry,'precomputed') && entry.precomputed
        delta = entry.delta; kbar = entry.kbar;
        CE_enter = entry.CE; U_enter = entry.U;
    else
        [delta, kbar, CE_enter, U_enter] = compute_entry_probability_bisect( ...
            cfg, sigma, Fh, fh, supp);
    end

    
    % (3) Threshold entry rule: enter iff k <= k_bar = CE(δ).
    entrants = (k <= kbar);
    n = sum(entrants);
    
    %  if there are < 2 entrants (i.e., n == 0), SKIP this auction.
    if n < 2
        % Return complete structure with default/NaN values for consistency
        out = struct();
        out.skipped   = true;
        out.n         = n;
        out.delta     = delta;
        out.kbar      = kbar;
        out.win_s     = NaN;       % No winner
        out.win_b0    = NaN;       
        out.win_b2    = NaN;       
        out.e0        = NaN;
        out.e2        = NaN;
        out.epsilon   = 0.0;       
        out.pay_final = NaN;
        out.overrun   = NaN;
        out.w_idx     = NaN;
        
        % Create empty panel with correct structure
        out.panel = struct('entered', entrants, ...  % use actual entrants (could be false(N,1))
                          'k', k(:), ...
                          'e0', nan(N,1), 'e2', nan(N,1), 'c', nan(N,1), ...
                          's', nan(N,1), 'b0', nan(N,1), 'b2', nan(N,1), ...
                          'winner', false(N,1));
        
        if cfg.trace.enabled
            out.trace = struct('riskTag', riskTag, 'U_enter', U_enter, 'CE_enter', CE_enter, ...
                               'delta', delta, 'kbar', kbar, 'n', n);
        else
            out.trace = [];
        end
        return;
    end
    
    % (4) Draw types only for entrants; compute c (Eq. 5)
    [e0_in, e2_in] = draw_types_e0e2(n, cfg.type.sd_e0, cfg.type.sd_e2, cfg.type.rho_e);
    %     Pseudo-cost (Eq. (5), 1D): c = θ0*e0 + θ1 - (e2-1)^2 / (2 α σ^2).
    c_in = th0.*e0_in + th1 - ((e2_in - 1).^2) ./ (2*alp*sigma^2);
    
    
    % (5) Scoring strategy s(c;n) from cache or build (Eq. 7)
    s_handle = get_s_handle(n);
    s_in = s_handle(c_in);
    
    % (6) Inner skew (Eq. (4), 1D): b2* = θ1 + (e2-1)/(α σ^2), truncated to [0, s]; b0* = s - b2*.
    b2_in = compute_optimal_skew_1d(s_in, e2_in, th1, alp, sigma);
    b0_in = s_in - b2_in;
    
    % (7) Winner (lowest score). Implementation shock ε ~ N(0,σ^2).
    %     Final payment (Eq. (2), 1D): p = b0 + b2 * (e2 + ε).
    %     Overrun = p - s = b2 * ( (e2 + ε) - 1 ).
    [win_s, w_local] = min(s_in);           % index within entrants
    epsi      = sigma * randn();
    pay_final = b0_in(w_local) + b2_in(w_local) * (e2_in(w_local) + epsi);
    overrun   = pay_final - win_s;
    
    % pack outputs
    out.skipped   = false;
    out.n         = n;
    out.delta     = delta;
    out.kbar      = kbar;
    out.win_s     = win_s;
    out.win_b0    = b0_in(w_local);
    out.win_b2    = b2_in(w_local);
    out.e0        = e0_in(w_local);
    out.e2        = e2_in(w_local);
    out.epsilon   = epsi;
    out.pay_final = pay_final;
    out.overrun   = overrun;
    out.w_idx     = w_local;                 % index among entrants
    
    % ---------------- build full-N panel (store ALL bidders) ----------------
    % Arrays of length N; fill entrants slots, NaN for non-entrants
    panel_e0 = nan(N,1); panel_e2 = nan(N,1); panel_c  = nan(N,1);
    panel_s  = nan(N,1); panel_b0 = nan(N,1); panel_b2 = nan(N,1);
    panel_entered = entrants;                % logical N×1
    panel_k  = k(:);                         % store entry costs for all

    % Map entrants in order to positions
    idx_all = find(entrants);                % positions in 1..N that entered
    panel_e0(idx_all) = e0_in;
    panel_e2(idx_all) = e2_in;
    panel_c(idx_all)  = c_in;
    panel_s(idx_all)  = s_in;
    panel_b0(idx_all) = b0_in;
    panel_b2(idx_all) = b2_in;

    % Winner flag across all N (only one true among entrants)
    panel_winner = false(N,1);
    panel_winner(idx_all(w_local)) = true;

    out.panel = struct('entered', panel_entered, ...
                    'k', panel_k, ...
                    'e0', panel_e0, 'e2', panel_e2, 'c', panel_c, ...
                    's', panel_s, 'b0', panel_b0, 'b2', panel_b2, ...
                    'winner', panel_winner);

    if cfg.trace.enabled
        out.trace = struct('riskTag', riskTag, 'U_enter', U_enter, 'CE_enter', CE_enter, ...
                        'delta', delta, 'kbar', kbar, 'n', n, ...
                        'c_stats', [min(c_in) mean(c_in) max(c_in)]);
    else
        out.trace = [];
    end
end
    