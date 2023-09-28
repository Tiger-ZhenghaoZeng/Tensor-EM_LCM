function y = GIC(p,theta,sample,type)
% calculate the GIC value 
% Input: p the proportion vector, theta the item parameters, sample a N*J
% dataset, type = 1 is the case of BIC, type = 2 is the case of GIC with 
% aN=log(log(N))*log(N) 
% Output: GIC value
[J,K] = size(theta);
[N,~] = size(sample);
if type == 1
    an = log(N);
else
    an = log(log(N))*log(N);
end
like = likelihood(sample,p,theta);
y = -2*like + an*(J*K+K-1);
end

