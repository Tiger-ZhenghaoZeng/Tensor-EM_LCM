function y = evaluate_error(z, z_hat, L)
% Evaluate error rates of estimated membbership z_hat when the true
% membership is z with L latent classes.
    n_error = 0;
    for k = 1:L
        index = find(z_hat == k);
        majority = mode(z(index));
        n_error = n_error + sum(z(index)~=majority);
    end
    y = n_error/length(z);
    