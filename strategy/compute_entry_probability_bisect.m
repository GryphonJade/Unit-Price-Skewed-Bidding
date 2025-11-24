function [delta, kbar, CE_enter, U_enter, trace] = compute_entry_probability_bisect( ...
    cfg, sigma, Fh, fh, supp)
% UP: Endogenous entry δ via robust bisection on
%   h(δ) = Φ((CE(δ)-μ_k)/σ_k) - δ,
% where CE(δ) is computed from the WIN-PROB-WEIGHTED expected utility:
%   U_n = ∫ u( s(c;n)-c ) * [1-F(c)]^(n-1) * f(c) dc,  u(x)=-exp(-α x),
% and we condition on n >= 2 to avoid the n=1 degeneracy.
%
% Inputs:
%   cfg   : config (must have cfg.eff.alpha from Eq.(18))
%   sigma : project risk (std), kept constant per Eq.(18)
%   Fh,fh : CDF/PDF handles of UP pseudo-cost c (KDE-based)
%   supp  : support [L,R] used for integration
%
% Outputs:
%   delta    : equilibrium entry probability
%   kbar     : entry threshold CE(delta)
%   CE_enter : CE(delta)
%   U_enter  : U(delta)  (negative)
%   trace    : struct with diagnostic info

% ----- primitives -----
N     = cfg.sim.N_potential;
alpha = cfg.eff.alpha;                       % EFFECTIVE α (scaled per Eq.18)

% use calibrated entry costs if present
if isfield(cfg.entry,'mu_k_eff'), muK = cfg.entry.mu_k_eff; else, muK = cfg.entry.mu_k; end
if isfield(cfg.entry,'sd_k_eff'), sdK = cfg.entry.sd_k_eff; else, sdK = cfg.entry.sd_k; end

% cache E[u] per n (independent of δ once s(c;n) fixed)
EU_cache = containers.Map('KeyType','double','ValueType','double');

rootFun = @(d) Phi_of_CE_UP(d, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache) - d;

% ----- bracketing on (eps, 1-eps) -----
epsL = cfg.entry.bisect.bracket(1);
epsR = cfg.entry.bisect.bracket(2);
hTol = 1e-8;                                   % near-zero tolerance

fa = rootFun(epsL);
fb = rootFun(epsR);
if ~isfinite(fa) || ~isfinite(fb)
    error('Non-finite h(delta) at bracket endpoints.');
end

% accept near-zero endpoints
if abs(fa) < hTol
    delta = epsL;
    [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
    if nargout>=5, trace = struct('it',0,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','endpoint-epsL'); end
    return;
end
if abs(fb) < hTol
    delta = epsR;
    [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
    if nargout>=5, trace = struct('it',0,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','endpoint-epsR'); end
    return;
end

% if already sign-changed, set bracket; else do coarse scan
if fa*fb < 0
    a = epsL; b = epsR;
else
    % coarse vectorized scan
    Kscan = 200;
    gridD = linspace(epsL, epsR, Kscan);
    hvals = arrayfun(@(d) safe_eval(rootFun,d), gridD);

    sgn = sign(hvals); sgn(abs(hvals)<hTol) = 0;
    finiteMask = isfinite(hvals);

    iz = find(sgn==0 & finiteMask, 1, 'first');
    if ~isempty(iz)
        delta = gridD(iz);
        [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
        if nargout>=5, trace = struct('it',0,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','grid-zero'); end
        return;
    end

    idx = [];
    for i = 1:Kscan-1
        if finiteMask(i) && finiteMask(i+1) && sgn(i)*sgn(i+1) < 0
            idx = i; break;
        end
    end
    if ~isempty(idx)
        a = gridD(idx); b = gridD(idx+1);
    else
        [~, imin] = min(abs(hvals(finiteMask)));
        idxFinite = find(finiteMask);
        im = idxFinite(imin);
        L = im; R = im; found = false;
        while (~found) && (L>1 || R<Kscan)
            if L>1, L=L-1; end
            if R<Kscan, R=R+1; end
            if finiteMask(L) && finiteMask(R) && sgn(L)*sgn(R) < 0
                a = gridD(L); b = gridD(R); found = true; break;
            end
            if L>1 && finiteMask(L-1) && sgn(L-1)*sgn(L) < 0
                a = gridD(L-1); b = gridD(L); found = true; break;
            end
            if R<Kscan && finiteMask(R+1) && sgn(R)*sgn(R+1) < 0
                a = gridD(R); b = gridD(R+1); found = true; break;
            end
        end
        if ~found
            if hvals(end) >= 0
                delta = epsR;   % δ* ≈ 1-eps
                [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
                if nargout>=5, trace = struct('it',0,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','fallback δ≈1-eps'); end
                return;
            elseif hvals(1) <= 0
                delta = epsL;   % δ* ≈ eps
                [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
                if nargout>=5, trace = struct('it',0,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','fallback δ≈eps'); end
                return;
            else
                error('Entry root not bracketed after robust scan.');
            end
        end
    end
end

% ----- bisection on [a,b] -----
fa = rootFun(a); fb = rootFun(b);
itMax = cfg.entry.bisect.maxIter;
tol   = cfg.entry.bisect.tolRoot;

for it = 1:itMax
    m  = 0.5*(a+b);
    fm = rootFun(m);
    if abs(fm) < tol || (b-a) < 1e-6
        delta = max(min(m, 1-1e-8), 1e-8);
        [Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
        if nargout>=5, trace = struct('it',it,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter); end
        return;
    end
    if fa*fm <= 0
        b = m; fb = fm;
    else
        a = m; fa = fm;
    end
end

% fallback
delta = max(min(0.5*(a+b), 1-1e-8), 1e-8);
[Phi_val, kbar, CE_enter, U_enter] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache);
if nargout>=5, trace = struct('it',itMax,'delta',delta,'Phi',Phi_val,'CE',CE_enter,'U',U_enter,'note','maxIter'); end
end


% =================== helpers =================== %
function [Phi_val, kbar, CE, U] = Phi_of_CE_UP(delta, cfg, sigma, N, alpha, muK, sdK, Fh, fh, supp, EU_cache)
% Compute CE(δ) for UP (with win-prob weighting), then Φ((CE-μ)/σ)

% weights over n>=2 (conditional) in log domain
logw = -inf(1,N);
ld   = log(max(delta,      1e-15));
l1d  = log(max(1 - delta,  1e-15));
for n = 2:N
    % log C(N-1,n-1) = gammaln(N) - gammaln(n) - gammaln(N-n+1)
    logC     = gammaln(N) - gammaln(n) - gammaln(N-n+1);
    logw(n)  = logC + (n-1)*ld + (N-n)*l1d;
end
maxlog = max(logw(2:end));
w      = exp(logw - maxlog);
Z      = sum(w(2:end));
if Z<=0 || ~isfinite(Z)
    U = -1; CE = -(1/alpha)*log(-U); kbar = CE; Phi_val = normcdf((CE - muK)/sdK); return;
end

U = 0.0;
for n = 2:N
    pn = w(n)/Z; if pn < 1e-14, continue; end
    key = n + 0.0;
    if isKey(EU_cache, key)
        EU_n = EU_cache(key);
    else
        % UP scoring strategy for this n (ODE (7) solver)
        s_handle = build_s_handle(cfg, Fh, fh, supp, n);
        % WIN-PROB utility:
        EU_n = expected_utility_integral_weighted(cfg, s_handle, Fh, fh, supp, alpha, n);
        EU_cache(key) = EU_n;
    end
    U = U + pn * EU_n;
end

U   = min(U, -1e-12);                  % ensure in (-∞,0)
CE  = -(1/alpha) * log(-U);            % inverse funtion of CARA (Certainty Equity/Money Value)
kbar = CE;
Phi_val = normcdf((CE - muK)/sdK);
end

function val = safe_eval(fun, x)
% safe wrapper to avoid exceptions/NaNs in coarse scan
try
    val = fun(x);
    if ~isfinite(val), val = NaN; end
catch
    val = NaN;
end
end
