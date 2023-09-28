function [lambda,evectors] = eigen(T,inits,N,L)
% Obtain all the eigenpairs of orthogonal decomposable tensor T
% Input: T the tensor, inits the initial values of tensor power
% method, L the number of latent classes, N the number of iterations
% Ouput: lambda the proportion vectors, evectors the eigenvectors
    d = size(T);
    d = d(1);
    vectors = zeros(d,L);
    Lambda = zeros(1,L);
    for i = 1:L
        [l,v,T] = eigenpairs(T,inits(:,i),N);
        Lambda(i) = l;
        vectors(:,i) = v;
    end
    lambda = Lambda;
    evectors = vectors;
    
    
