function [Lambda,V1,T1] = eigenpairs(T,init,n)
% Obtain an eigenpair of T and the deflation
% Input: T the tensor from which we want to obtain an eigenpair. init the
% initial value for tensor power method. n the number of iterations.
% Output: Lambda the eigenvalue, V1 the eigenvector, T1 the deflated tensor
d = size(T);
d = d(1);
theta = init;
theta = theta/norm(theta);

% rng(314);
for j = 1:n
    a = double(ttm(T,{eye(d),theta',theta'},[1,2,3]));
    theta = a/norm(a);
end
for j = 1:n
    a = double(ttm(T,{eye(d),theta',theta'},[1,2,3]));
    theta = a/norm(a);
end
Lambda = double(ttm(T,{theta',theta',theta'},[1,2,3]));
V1 = theta;
S = symktensor(Lambda,theta,3);
T1 = T - full(S);
