function [lambda,evectors] = eigen(T,K,N,L)
% Obtain all the eigenpairs of orthogonal decomposable tensor T
% Input: T the tensor, K number of random initial values of tensor power
% method, L the number of latent classes, N the number of iterations
% Ouput: lambda the proportion vectors, evectors the eigenvectors
    d = size(T);
    d = d(1);
    vectors = zeros(d,L);
    Lambda = zeros(1,L);
    for i = 1:L
        [l,v,T] = eigenpairs(T,K,N);
        Lambda(i) = l;
        vectors(:,i) = v;
    end
    lambda = Lambda;
    evectors = vectors;