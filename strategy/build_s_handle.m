function s_handle = build_s_handle(cfg, Fh, fh, supp, n)
    % Build s(c;n) by solving ODE (7) with boundary s(c_max)=c_max:
    % ds/dc = 1 + ((n-1) f(c) / (α [1-F(c)])) * (exp(α (s(c)-c)) - 1)
    alpha = cfg.eff.alpha;
    opts  = cfg.solver.ode;
    s_handle = solve_s_via_ode(Fh, fh, supp, n, alpha, opts);
end