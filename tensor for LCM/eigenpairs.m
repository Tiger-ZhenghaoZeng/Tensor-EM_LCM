function [Lambda,V1,T1] = eigenpairs(T,K,n)
% Obtain an eigenpair of T and the deflation
% Input: T the tensor from which we want to obtain an eigenpair. K the
% number of random initial values for tensor power method. n the number of
% iterations.
% Output: Lambda the eigenvalue, V1 the eigenvector, T1 the deflated tensor
d = size(T);
d = d(1);
theta = zeros(d,K);
% rng(314);
for i = 1:K
    theta(:,i) = randn(d,1);
    theta(:,i) = theta(:,i)/norm(theta(:,i));
    for j = 1:n
        a = double(ttm(T,{eye(d),theta(:,i)',theta(:,i)'},[1,2,3]));
        theta(:,i) = a/norm(a);
    end
end
lambda = zeros(1,K);
for i = 1:K
    lambda(i) = double(ttm(T,{theta(:,i)',theta(:,i)',theta(:,i)'},[1,2,3]));
end
max = 1;
for i = 2:K
    if(lambda(i)>lambda(max))
        max = i;
    end
end
theta = theta(:,max);
for j = 1:n
    a = double(ttm(T,{eye(d),theta',theta'},[1,2,3]));
    theta = a/norm(a);
end
Lambda = double(ttm(T,{theta',theta',theta'},[1,2,3]));
V1 = theta;
S = symktensor(Lambda,theta,3);
T1 = T - full(S);
