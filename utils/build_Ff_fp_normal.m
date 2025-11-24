% utils/build_Ff_fp_normal.m
function [Fh, fh, supp, moments] = build_Ff_fp_normal(cfg, sigma)
    % Build exact Normal F and f for c^FP = theta0*e0 + theta1*e2 + (alpha/2)*theta1^2*sigma^2
    % using EFFECTIVE (scaled) parameters per Eq.(18): theta_j*s, alpha/s.
    
    % effective (scaled) parameters
    alp = cfg.eff.alpha;
    t0  = cfg.eff.theta0;
    t1  = cfg.eff.theta1;
    
    % BVN(e0,e2): means=1,1; sds & corr from cfg.type
    m0  = 1.0; m2 = 1.0;
    sd0 = cfg.type.sd_e0;
    sd2 = cfg.type.sd_e2;
    rho = cfg.type.rho_e;
    
    % mean and variance of linear combo L = t0*e0 + t1*e2
    mu_L  = t0*m0 + t1*m2;
    var_L = (t0^2)*(sd0^2) + (t1^2)*(sd2^2) + 2*t0*t1*rho*sd0*sd2;
    
    % risk premium (constant): (alpha/2)*t1^2*sigma^2
    RP   = 0.5 * alp * (t1^2) * (sigma^2);
    
    mu_c = mu_L + RP;
    sd_c = sqrt(var_L);
    
    % handles
    Fh = @(x) normcdf((x - mu_c)./sd_c);
    fh = @(x) normpdf((x - mu_c)./sd_c) ./ sd_c;
    
    % support (±6σ) 
    L = mu_c - 6*sd_c;
    R = mu_c + 6*sd_c;
    supp = [L, R];
    
    if nargout >= 4
        moments = struct('mu',mu_c,'sd',sd_c,'RP',RP,'mu_L',mu_L,'sd_L',sqrt(var_L));
    end
    end
    