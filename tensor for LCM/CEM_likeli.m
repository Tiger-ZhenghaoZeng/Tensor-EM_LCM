function y = CEM_likeli(sample, z, Thetas)
% Calculate the log-likelihood in CEM algorithm
% Input: sample a N*J dataset, z the latent class assignments, Thetas the item
% parameters
% Output: log-likelihood
    R = sample';
    [J,N] = size(R);
    s = 0;
    for i = 1:N
        q = 1;
        for j = 1:J
            q = q * Thetas(j,z(i))^(R(j,i))*(1-Thetas(j,z(i)))^(1-R(j,i)); 
        end
        s = s+log(q);
    end
    y = s;
end




