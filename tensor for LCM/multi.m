function [w,X1,X2,X3] = multi(sample,L,d,K,n)
% Use tensor method to learn parameters from latent class model
% Input: sample is a N by J dataset, L is the number of classes, d is the
% division of J items (for example J = 90, then d can be set to
% [30,30,30]), K the number of initial values for tensor power method, n
% the iteration times.
% Output: w is the proportion vector, X1,X2,X3 is the item parameters
% corresponding to the division specified by d.
    [N,~] = size(sample);
    d = [0,d];
    index = cumsum(d);
    sample1 = sample(:,(index(1)+1):index(2));
    sample2 = sample(:,(index(2)+1):index(3));
    sample3 = sample(:,(index(3)+1):index(4));
    pairs13 = 0;
    for i = 1:N
        pairs13 = pairs13 + sample1(i,:)'*sample3(i,:);
    end
    pairs13 = pairs13/N;
    pairs23 = 0;
    for i = 1:N
        pairs23 = pairs23 + sample2(i,:)'*sample3(i,:);
    end
    pairs23 = pairs23/N;
    pairs12 = 0;
    for i = 1:N
        pairs12 = pairs12 + sample1(i,:)'*sample2(i,:);
    end
    pairs12 = pairs12/N;
    C21 = pairs13*MPinv(pairs23,L);
    C12 = pairs23*MPinv(pairs13,L);
    C31 = pairs12*MPinv(pairs23',L);
    C13 = (pairs23')*MPinv(pairs12,L);
    sample21 = (C21*sample2')';
    sample31 = (C31*sample3')';
    M2 = 0;
    for i = 1:N
        M2 = M2 + sample1(i,:)'*sample21(i,:);
    end
    M2 = M2/N;
    %M3 = ktensor({sample1',sample21',sample31'});
    %M3 = full(M3)/N; 
    M3 = tensor(zeros(d(2),d(2),d(2)));
    for i = 1:N
        M3 = M3 + full(ktensor({sample1(i,:)',sample21(i,:)',sample31(i,:)'}));
    end
    M3 = M3/N;
    [w,X1] = Whiten(M2,M3,K,n,L);
    X2 = C12*X1;
    X3 = C13*X1;
end


