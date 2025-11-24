function ce = cara_ce(U, alpha)
    % CARA certainty equivalent: CE = u^{-1}(U) = -(1/α) * log(-U)
    ce = -(1/alpha) * log(-U);
end