function cfg = config_up_1d()
    % =================== Hyper-parameters & switches =================== %
    cfg.model.alpha   = 7.98;       % CARA risk aversion α
    cfg.model.theta0  = 0.237;      % cost (lump-sum part), θ0
    cfg.model.theta1  = 0.756;      % cost (unit-price part), θ1
    cfg.model.sigma_L = 0.0837;     % project risk (low),  σ_L
    cfg.model.sigma_H = 0.0995;     % project risk (high), σ_H

    % type distribution for (e0,e2): BVN(mean=1,1), sds & corr
    cfg.type.sd_e0 = 0.982;
    cfg.type.sd_e2 = 0.0140;
    cfg.type.rho_e = 0.758;

    % entry cost k ~ Normal(mu_k, sd_k^2)  (NOT lognormal)
    cfg.entry.mu_k = 0.0181;
    cfg.entry.sd_k = 0.0232;

    % population size & simulation control
    cfg.sim.N_potential = 9;     % potential bidders (ex ante)
    cfg.sim.nAuctions   = 200;    % number of auctions

    % ---------- NEW: stochastic project risk mixing ----------
    % If useProb=true, each auction draws risk state H with prob P, else L with 1-P.
    % If useProb=false, we use riskSchedule (e.g., {'L','H'} alternating) as before.
    cfg.sim.useProb = true;       % NEW switch: Bernoulli(P) over {L,H}
    cfg.sim.P       = 0.449;        % NEW: P = Pr(state = H)
    cfg.sim.riskSchedule = {'L','H'};  % used only if useProb=false

    % entry solver (bisection) setup
    cfg.entry.bisect.maxIter = 100;
    cfg.entry.bisect.tolRoot = 1e-8;
    cfg.entry.bisect.bracket = [1e-8, 1-1e-8];

    cfg.scale.useProjectScale = true;
    cfg.scale.s = 0.15;  

    % kernel / grid for F,f of c (pseudo-cost)
    cfg.kernel.M_c    = 60000;   % samples for c to build F,f 
    % in config_up_1d.m
    cfg.strategy.Ngrid_bfp = 2000;   % grid points for FP bidding curve (fast & smooth)
    cfg.kernel.bw_mult = 1.0;    % bandwidth multiplier (1.0 = Silverman's rule)

    % strategy solver / integral tolerances
    cfg.solver.ode.RelTol = 1e-7;
    cfg.solver.ode.AbsTol = 1e-8;
    cfg.solver.integral.RelTol = 1e-7;
    cfg.solver.integral.AbsTol = 1e-8;

    % tracing / debugging
    cfg.trace.enabled     = true;     % store per-auction trace
    cfg.trace.verbose     = true;     % print progress
    cfg.trace.printEvery  = 10;       % print every k auctions
    cfg.trace.store_s_grid = false;   % optionally store s(c) grid for realized n (large)
    
    % caching options
    cfg.cache.precomputeEntry = true;  % precompute entry probabilities for each risk state
end