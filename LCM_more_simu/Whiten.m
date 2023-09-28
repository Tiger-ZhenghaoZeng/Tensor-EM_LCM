function [omega,Mu]=Whiten(M2,M3,inits,n,L,p)
% Whitening M3 with M2 and obtain decomposition
% Input: K the number of random
% initials of tensor power method, n the number of iterations, L the number of classes, p true proportion vector. M2 a second-order tensor
% M3 a third-order tensor specified in Tensor Decompositions for Learning Latent Variable Models
% Output: omega the proportion vector, Mu the item parameters
d = length(M2);
[U,D,V] = svd(M2);
%if (rank(D)~=k)
 %  error('M2 must have rank %d', k)
%end
U = U(:,1:L);
D = D(1:L,1:L);
W = U*D^(-0.5);
M3tilde = ttm(M3,{W',W',W'},[1,2,3]);
init_transformed = zeros(L,L);
for i=1:L
    init_transformed(:,i) = sqrt(p(i))*W'*inits(:,i);
end
[lambda,evectors] = eigen(M3tilde,init_transformed,n,L);
B = pinv(W');
mu = zeros(d,L);
for i = 1:L
    mu(:,i) = lambda(i)*B*evectors(:,i);
end
omega = 1./(lambda.^2);
Mu = mu;

%{
function [omega,Mu]=Whiten(M2,M3,K,n,L)
% Whitening M3 with M2 and obtain decomposition
% Input: K the number of random
% initials of tensor power method, n the number of iterations, L the number of classes. M2 a second-order tensor
% M3 a third-order tensor specified in Tensor Decompositions for Learning Latent Variable Models
% Output: omega the proportion vector, Mu the item parameters
d = length(M2);
[U,D,V] = svd(M2);
%if (rank(D)~=k)
 %  error('M2 must have rank %d', k)
%end
U = U(:,1:L);
D = D(1:L,1:L);
W = U*D^(-0.5);
M3tilde = ttm(M3,{W',W',W'},[1,2,3]);
[lambda,evectors, inits] = eigen(M3tilde,K,n,L);
B = pinv(W');
mu = zeros(d,L);
Inits = zeros(d,L);
for i = 1:L
    mu(:,i) = lambda(i)*B*evectors(:,i);
    Inits(:,i) = lambda(i)*B*inits(:,i);
end
omega = 1./(lambda.^2);
Mu = mu;
init_values = Inits;
%}