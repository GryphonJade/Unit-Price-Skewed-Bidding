function k = draw_entry_costs_norm(N, mu, sd)
    % Entry cost ~ Normal(mu, sd^2), truncated at 0
    k = mu + sd .* randn(N,1);
    k(k<0) = 0;
end