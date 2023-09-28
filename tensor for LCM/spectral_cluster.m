function [idx] = spectral_cluster(sample,k)
% (normalized random-walk) Spectral Clustering with k clusters using given sample
    dist_temp = pdist(sample, 'hamming');
    dist = squareform(dist_temp);
    S = exp(-dist.^2);
    D = diag(sum(S));
    L = D-S;
    opt = struct('IsFunctionSymmetric', true, 'isreal', true);
    [V, ~] = eigs(L, D, k, 'smallestabs', opt);
    idx = kmeans(V, k);
end








