function c = sample_c_via_sobol(M, sigma, alpha, cfg)
    % Sample pseudo-cost c via Eq. (5) from BVN (e0,e2) with mean 1
    th0 = cfg.eff.theta0; th1 = cfg.eff.theta1;
    sd0 = cfg.type.sd_e0;   sd2 = cfg.type.sd_e2;  rho = cfg.type.rho_e;
    
    S = [sd0^2, rho*sd0*sd2; rho*sd0*sd2, sd2^2];
    % Sobol quasi-random for variance reduction
    p = sobolset(2,'Skip',1e3,'Leap',1e3); p = scramble(p,'MatousekAffineOwen');
    U = net(p, M);
    Z = icdf('Normal', U, 0, 1);
    L = chol(S, 'lower');
    E = (L * Z')';
    e0 = 1 + E(:,1); e2 = 1 + E(:,2);
    e0(e0<=0)=eps; e2(e2<=0)=eps;
    
    c = th0.*e0 + th1 - ((e2 - 1).^2) ./ (2*alpha*sigma^2);
end