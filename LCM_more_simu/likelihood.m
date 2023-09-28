function y = likelihood(sample,p,Thetas)
% Calculate the log-likelihood 
% Input: sample a N*J dataset, p the proportion vector, Thetas the item
% parameters
% Output: log-likelihood
    R = sample';
    [~,N] = size(R);
    G = length(p);
    s = 0;
    for i = 1:N
        q = 0;
        for l = 1:G
            q = q + p(l)*prod(Thetas(:,l).^(R(:,i)))*prod((1-Thetas(:,l)).^(1-R(:,i)));
        end
        s = s+log(q);
    end
    y = s;
end

