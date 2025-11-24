function cf_run_fixed_entry(n_target)
    % Run UP vs FP counterfactual with the SAME number of entrants n_target
    % across many auctions, then summarize expected final payments, etc.
    %
    % Usage:
    %   cf_run_fixed_entry(5);  % for example
    
    if nargin<1, n_target = 5; end
    
    addpath(genpath('config')); addpath(genpath('utils'));
    addpath(genpath('strategy')); addpath(genpath('sim')); addpath(genpath('stats'));
    addpath(genpath('counterfactual'));
    
    cfg = config_up_1d();
    
    % derive effective (scaled) parameters per Eq.(18)
    if isfield(cfg,'scale') && isfield(cfg.scale,'useProjectScale') && cfg.scale.useProjectScale
        s = cfg.scale.s;
    else
        s = 1.0;
    end
    cfg.eff.theta0 = cfg.model.theta0 * s;
    cfg.eff.theta1 = cfg.model.theta1 * s;
    cfg.eff.alpha  = cfg.model.alpha  / s;
    cfg.eff.sigma_L = cfg.model.sigma_L; cfg.eff.sigma_H = cfg.model.sigma_H;
    
    % --------- Pre-build F_u (UP) and F_fp (FP) per risk ----------
    riskTags = {'L','H'};
    FF_UP = struct(); FF_FP = struct();
    for r = 1:numel(riskTags)
        tag   = riskTags{r};
        sigma = cfg.eff.(['sigma_' tag]);
        % UP pseudo-cost via sampling + KDE
        c_emp = sample_c_via_sobol(cfg.kernel.M_c, sigma, cfg.eff.alpha, cfg);
        [Fh_u, fh_u, supp_u] = build_Ff_splines(c_emp, cfg.kernel.bw_mult);
        FF_UP.(tag) = struct('Fh',Fh_u,'fh',fh_u,'supp',supp_u,'sigma',sigma);
    
        % FP type exact normal
        [Fh_f, fh_f, supp_f] = build_Ff_fp_normal(cfg, sigma);
        FF_FP.(tag) = struct('Fh',Fh_f,'fh',fh_f,'supp',supp_f,'sigma',sigma);
    end
    
    % --------- strategy handle caches ----------
    SCacheUP.L = containers.Map('KeyType','double','ValueType','any');
    SCacheUP.H = containers.Map('KeyType','double','ValueType','any');
    SCacheFP.L = containers.Map('KeyType','double','ValueType','any');
    SCacheFP.H = containers.Map('KeyType','double','ValueType','any');
    
    get_s_handle = @(tag,n) get_or_build_s_handle(cfg, FF_UP.(tag).Fh, FF_UP.(tag).fh, ...
                                                  FF_UP.(tag).supp, n, tag, SCacheUP);
    get_b_handle = @(tag,n) get_or_build_bfp_handle(cfg, FF_FP.(tag).Fh, FF_FP.(tag).fh, ...
                                                    FF_FP.(tag).supp, n, tag, SCacheFP);
    
    % --------- loop over auctions ----------
    nAuctions = cfg.sim.nAuctions;
    A = table('Size',[0 9], 'VariableTypes', ...
              {'double','string','double','double','double','double','double','double','double'}, ...
              'VariableNames', {'auction_id','risk','sigma','n', ...
                                'up_win_score','up_pay_final','fp_win_score','fp_pay_final','scale_s'});
    
    t0 = tic;
    auction_id = 0;
    for t = 1:nAuctions
        % draw risk by P or schedule
        if cfg.sim.useProb
            tag = (rand < cfg.sim.P) * 'H' + (rand >= cfg.sim.P) * 'L';
            if ~ischar(tag), tag = ternary(rand<cfg.sim.P,'H','L'); end
        else
            sch = cfg.sim.riskSchedule; tag = sch{mod(t-1,numel(sch))+1};
        end
        sigma = cfg.eff.(['sigma_' tag]);
    
        auction_id = auction_id + 1;
        out = cf_fixed_entry_one(cfg, tag, sigma, n_target, ...
                                 FF_UP.(tag).Fh, FF_UP.(tag).fh, FF_UP.(tag).supp, ...
                                 FF_FP.(tag).Fh, FF_FP.(tag).fh, FF_FP.(tag).supp, ...
                                 @(n)get_s_handle(tag,n), @(n)get_b_handle(tag,n));
    
        A = [A; {auction_id, string(tag), sigma, n_target, ...
                 out.UP.win_score, out.UP.pay_final, ...
                 out.FP.win_score, out.FP.pay_final, s}];
    end
    elapsed = toc(t0);
    fprintf('[CF-fixed-entry] done in %.2fs: n_target=%d, auctions=%d\n', elapsed, n_target, nAuctions);
    
    % --------- write outputs ----------
    if ~exist('output','dir'), mkdir('output'); end
    writetable(A, fullfile('output', sprintf('cf_fixed_entry_n%d.csv', n_target)));
    
    % --------- small summaries (homogenized vs original are the same if s=1) ----------
    fprintf('== Fixed-entry (n=%d): mean by risk ==\n', n_target);
    S = groupsummary(A, 'risk', 'mean', {'up_win_score','up_pay_final','fp_win_score','fp_pay_final'});
    disp(S);
    
    % -------- helper functions --------
    function s_handle = get_or_build_s_handle(cfg, Fh, fh, supp, n, tag, SCache)
        % Return cached s(c;n) for (risk=tag, n); otherwise build & store.
        m = SCache.(tag);
        key = n + 0.0;
        if isKey(m, key)
            s_handle = m(key);
        else
            s_handle = build_s_handle(cfg, Fh, fh, supp, n);  % ODE (7) solver
            m(key)   = s_handle;
            SCache.(tag) = m;
        end
    end

    function y = ternary(cond,a,b), if cond, y=a; else, y=b; end; end

    % strategy/get_or_build_bfp_handle.m
    function b_handle = get_or_build_bfp_handle(cfg, Fh, fh, supp, n, tag, SCacheFP)
        m = SCacheFP.(tag);
        key = n + 0.0;
        if isKey(m,key)
            b_handle = m(key);
        else
            b_handle = build_bfp_handle(cfg, Fh, fh, supp, n);
            m(key)   = b_handle;
            SCacheFP.(tag) = m;
        end
        end


end