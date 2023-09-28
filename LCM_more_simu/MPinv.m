function y = MPinv(A,r)
% Calculate the M-P pseudoinverse of A with only first r singular values
    [U,D0,V]=svd(A);
    [m,n] = size(D0);
    D = zeros(n,m);
    for i = 1:r
        D(i,i) = D0(i,i)^(-1);
    end
    y = V*D*U';
end

