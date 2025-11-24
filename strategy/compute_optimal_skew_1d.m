function b2 = compute_optimal_skew_1d(s_vec, e2_vec, theta1, alpha, sigma)
    % inner-loop (Eq. (4)): b2* = θ1 + (e2-1)/(α σ^2), truncated to [0, s]
    b2 = theta1 + (e2_vec - 1) ./ (alpha * sigma^2);
    b2 = max(0, min(b2, s_vec));
end