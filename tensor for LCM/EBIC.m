function y = EBIC(p,theta,R,gam,input)
[~,N] = size(R);
card = length(p);
like = likelihood(R,p,theta);
y = -2*like + card*log(N) + 2*gam*log(input/card);