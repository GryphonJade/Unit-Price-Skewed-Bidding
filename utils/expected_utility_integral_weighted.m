function EU = expected_utility_integral_weighted(cfg, s_handle, Fh, fh, supp, alpha, n)
    % E[u | enter, n] = ∫ u(s(c;n) - c) * [1 - F(c)]^(n-1) * f(c) dc
    % u(x) = -exp(-alpha*x)
    
    integrand = @(c) (-exp(-alpha * (s_handle(c) - c))) ...
                     .* max(1 - Fh(c), 0).^(n-1) ...
                     .* fh(c);
    
    EU = integral(integrand, supp(1), supp(2), ...
                  'RelTol', cfg.solver.integral.RelTol, ...
                  'AbsTol', cfg.solver.integral.AbsTol);
    % EU < 0 as required by CE = -(1/alpha) log(-U)
    end
    