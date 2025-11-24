function s_handle = solve_s_via_ode(Fh, fh, supp, n, alpha, opts)
    % Solve s'(c) = 1 + (n-1) f(c) / [α (1-F(c))] * (exp(α (s(c)-c)) - 1)
    % with boundary s(cR)=cR, integrating from cR down to cL.
    
    cL = supp(1); cR = supp(2);
    rhs = @(c,s) 1 + (n-1) .* fh(c) ./ max(alpha*(1 - Fh(c)), 1e-10) .* (exp(alpha*(s - c)) - 1);
    [c_path, s_path] = ode45(@(c,s) rhs(c,s), [cR cL], cR, odeset(opts));
    % ensure monotone c ascending in interpolant
    [c_sorted, I] = sort(c_path,'ascend'); s_sorted = s_path(I);
    s_handle = @(cq) interp1(c_sorted, s_sorted, min(max(cq,c_sorted(1)),c_sorted(end)), 'pchip', 'extrap');
end