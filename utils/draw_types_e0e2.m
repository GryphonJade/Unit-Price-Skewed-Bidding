function [e0,e2] = draw_types_e0e2(N, sd_e0, sd_e2, rho)
    % Draw (e0,e2) ~ BVN with mean shifted to 1
    Sigma = [sd_e0^2, rho*sd_e0*sd_e2; rho*sd_e0*sd_e2, sd_e2^2];
    X = mvnrnd([0 0], Sigma, N);
    e0 = 1 + X(:,1);
    e2 = 1 + X(:,2);
    e0(e0<=0) = eps;
    e2(e2<=0) = eps;
end