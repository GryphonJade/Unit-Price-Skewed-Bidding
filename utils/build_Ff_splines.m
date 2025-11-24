function [Fh, fh, supp] = build_Ff_splines(c_emp, bw_mult)
    % Smooth PDF/CDF for c using KDE with reflection and spline interpolation.
    c_emp = sort(c_emp(:));
    if numel(c_emp) < 1000, warning('Few samples for F,f; increase M_c.'); end
    % bandwidth
    bw0 = 1.06 * std(c_emp) * numel(c_emp)^(-1/5);
    bw  = max(bw0 * bw_mult, eps);
    
    % extended grid & reflection
    cmin = c_emp(1); cmax = c_emp(end);
    L = cmin - 3*bw; R = cmax + 3*bw;
    xx = linspace(L,R,5000).';
    
    xR = 2*cmin - c_emp; xL = 2*cmax - c_emp;  % reflect
    c_reflect = [xR; c_emp; xL];
    [fpdf, xpdf] = ksdensity(c_reflect, xx, 'Bandwidth', bw, 'Function','pdf');
    fpdf = max(fpdf, 0);
    cdf = cumtrapz(xx, fpdf); cdf = cdf / max(cdf);
    
    Fh = @(x) interp1(xx, cdf, min(max(x,L),R), 'pchip');
    fh = @(x) interp1(xx, fpdf, min(max(x,L),R), 'pchip');
    supp = [L, R];
end