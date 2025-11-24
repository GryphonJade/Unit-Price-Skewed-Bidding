function b_handle = build_bfp_handle(cfg, Fh, fh, supp, n)
    % Build FP bidding function:
    %   b(c;n) = c + I(c)/W(c), where
    %   W(c) = (1 - F(c))^(n-1),
    %   I(c) = ∫_c^{R} W(t) dt  (compute via forward cumulative integral).
    alpha   = cfg.eff.alpha;
    opts    = cfg.solver.ode;
    
    % Solve s'(c) = 1 + ((n-1) f / (alpha (1-F))) * (exp(alpha (s-c)) - 1), s(cR)=cR
    b_handle = solve_s_via_ode(Fh, fh, supp, n, alpha, opts);
    end