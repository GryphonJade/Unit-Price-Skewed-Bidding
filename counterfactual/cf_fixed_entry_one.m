function out = cf_fixed_entry_one(cfg, riskTag, sigma, n, ...
    Fh_UP, fh_UP, supp_UP, ...
    Fh_FP, fh_FP, supp_FP, ...
    get_s_handle, get_b_handle)
    % Compare UP vs FP with the SAME entrants (same n and same (e0,e2) draws).
    % Returns expected final payment (E over eps=0) and winner stats for both.
    %
    % Inputs:
    %   cfg        : config (must have cfg.eff.theta*, cfg.eff.alpha)
    %   riskTag    : 'L' or 'H'
    %   sigma      : project risk (std)
    %   n          : fixed number of entrants
    %   Fh_UP,...  : UP pseudo-cost distribution (for s(c;n) ODE)
    %   Fh_FP,...  : FP type distribution (normal) used for building b(c;n) curve
    %   get_s_handle : @(n) -> s_handle for UP (cached builder)
    %   get_b_handle : @(n) -> b_handle for FP (cached builder)
    %
    % Output (struct):
    %   .n, .riskTag
    %   .UP.win_score, .UP.pay_final, .UP.winner_e0, .UP.winner_e2
    %   .FP.win_score, .FP.pay_final, .FP.winner_e0, .FP.winner_e2

    alp = cfg.eff.alpha;
    th0 = cfg.eff.theta0;
    th1 = cfg.eff.theta1;

    % ---- draw the SAME entrants' private types (e0,e2) ----
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
    % ===== quick diagnostics: compare mean types & RP =====
    % (a) unconditional means (UP type vs FP type) from the SAME entrants:
    mu_cu  = mean(c_u);
    RP     = 0.5 * alp * (th1^2) * (sigma^2);       % FP risk premium component
    mu_cfp = mean(c_fp);

    % (b) report the bbid-level averages
    fprintf('[CF-dbg] risk=%s, n=%d | E[c_u]=%.4f, E[c_fp]=%.4f (RP=%.4f)\n', ...
            riskTag, n, mu_cu, mu_cfp, RP);


    % ===================== pack =====================
    out = struct();
    out.n       = n;
    out.riskTag = riskTag;

    out.UP  = struct('win_score',win_s_u, 'pay_final',pay_u, ...
    'winner_e0',e0(wU), 'winner_e2',e2(wU));
    out.FP  = struct('win_score',win_b_fp, 'pay_final',pay_fp, ...
    'winner_e0',e0(wF),  'winner_e2',e2(wF));
    
end
