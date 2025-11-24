% sim/simulate_one_auction_full_fp.m
function out = simulate_one_auction_full_fp(cfg, riskTag, sigma, Fh, fh, supp, get_b_handle, entry)
    % One FP auction under the same environment: entry -> FP bidding -> implementation
    % FP has fixed payment, so "overrun for the agency" = 0 by definition.
    
    N   = cfg.sim.N_potential;
    alp = cfg.eff.alpha;
    t0  = cfg.eff.theta0;
    t1  = cfg.eff.theta1;
    
    % (1) Draw entry costs only
    k = draw_entry_costs_norm(N, cfg.entry.mu_k, cfg.entry.sd_k);
    
    % (2) Entry fixed point: precomputed vs compute now
    if isstruct(entry) && isfield(entry,'precomputed') && entry.precomputed
        delta = entry.delta; kbar = entry.kbar;
        CE_enter = entry.CE; U_enter = entry.U;
    else
        [delta, kbar, CE_enter, U_enter] = compute_entry_probability_bisect_fp( ...
            cfg, sigma, Fh, fh, supp);
    end
    
    % (3) Threshold entry
    entrants = (k <= kbar);
    n = sum(entrants);
    if n < 1
        % Return complete structure with default/NaN values for consistency
        out = struct();
        out.skipped   = true;
        out.n         = n;
        out.delta     = delta;
        out.kbar      = kbar;
        out.win_s     = NaN;       % No winner
        out.win_b0    = NaN;       
        out.win_b2    = 0.0;       
        out.e0        = NaN;
        out.e2        = NaN;
        out.epsilon   = 0.0;       
        out.pay_final = NaN;
        out.overrun   = 0.0;
        out.w_idx     = NaN;
        
        % Create empty panel with correct structure
        out.panel = struct('entered', false(N,1), ...
                          'k', k(:), ...
                          'e0', nan(N,1), 'e2', nan(N,1), 'c', nan(N,1), ...
                          's', nan(N,1), 'b0', nan(N,1), 'b2', zeros(N,1), ...
                          'winner', false(N,1));
        
        if cfg.trace.enabled
            out.trace = struct('riskTag', riskTag, 'U_enter', U_enter, 'CE_enter', CE_enter, ...
                               'delta', delta, 'kbar', kbar, 'n', n);
        else
            out.trace = [];
        end
        return;
    end
    
    % (4) Draw types only for entrants
    [e0_in, e2_in] = draw_types_e0e2(n, cfg.type.sd_e0, cfg.type.sd_e2, cfg.type.rho_e);
    
    % FP risk-adjusted cost type: c^FP = t0*e0 + t1*e2 + (alpha/2)*t1^2*sigma^2
    RP  = 0.5 * alp * (t1^2) * (sigma^2);
    c_in = t0.*e0_in + t1.*e2_in + RP;
    
    % (5) FP bidding function b(c;n) from cache or build (analytic integral)
    b_handle = get_b_handle(n);
    b_i      = b_handle(c_in);
    
    % (6) Winner & implementation
    [win_b, w] = min(b_i);
    epsi      = sigma * randn();                         % realized shock (no payment change)
    pay_final = win_b;                                    % agency pays fixed price in FP
    overrun   = 0.0;                                      % score = bid, no ex-post adjustment
    
    % pack outputs (auction-level)
    out.skipped   = false;
    out.n         = n;
    out.delta     = delta;
    out.kbar      = kbar;
    out.win_s     = win_b;     % for consistency, treat FP "score" as the bid
    out.win_b0    = win_b;     % b0 = bid
    out.win_b2    = 0.0;       % no skew term in FP
    out.e0        = e0_in(w);
    out.e2        = e2_in(w);
    out.epsilon   = epsi;      % realized shock (not used to pay)
    out.pay_final = pay_final;
    out.overrun   = overrun;
    out.w_idx     = w;
    
    % bidder-level panel (N rows for consistency with UP)
    panel_e0 = nan(N,1); panel_e2 = nan(N,1); panel_c = nan(N,1);
    panel_s  = nan(N,1); panel_b0 = nan(N,1); panel_b2 = nan(N,1);
    panel_entered = false(N,1);
    
    idx_all = find(entrants);
    panel_entered(idx_all) = true;
    panel_e0(idx_all) = e0_in;
    panel_e2(idx_all) = e2_in;
    panel_c(idx_all)  = c_in;
    panel_s(idx_all)  = b_i;       % score := bid
    panel_b0(idx_all) = b_i;       % b0 := bid
    panel_b2(idx_all) = 0.0;       % no skew
    
    panel_winner = false(N,1);
    panel_winner(idx_all(w)) = true;
    
    out.panel = struct('entered', panel_entered, ...
                       'k', k(:), ...
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
    