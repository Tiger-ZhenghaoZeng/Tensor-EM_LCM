function [w,v] = EM_random(sample,L,n)
% Latent class EM algorithm with n random initials
% Input: sample is a N*J dataset. L the number of latent class, n the
% number of initial points
% Output: w the proportion vectors, v the item parameters
    R = sample';
    [J,N] = size(R);
    like = -Inf;
    p_c = ones(1,L)/L;
    thetas_c = zeros(J,L);
    kk = 1;
    while kk <= n
        kk = kk+1;
        l0 = -Inf;
        p = ones(1,L)/L;
        Thetas = rand(J,L);
        Phi = zeros(N,L);
        Thetas1 = zeros(J,L);
        count = 0;
        while 1 
            count = count+1;
            for i = 1:N
                de = 0;
                for m = 1:L
                    de = de + p(m)*prod(Thetas(:,m).^(R(:,i)))*prod((1-Thetas(:,m)).^(1-R(:,i)));
                end
                for l = 1:L
                    Phi(i,l) = p(l)*prod(Thetas(:,l).^(R(:,i)))*prod((1-Thetas(:,l)).^(1-R(:,i)))/de;
                end
            end
            for l = 1:L
                p(l) = sum(Phi(:,l));
            end
            p = p./sum(p);
            for j = 1:J
                for g = 1:L
                    Thetas1(j,g) = sum(Phi(:,g).*(R(j,:)'))/(sum(Phi(:,g)));
                end
            end
            l1 = likelihood(sample,p,Thetas1);
            if abs(l1-l0) < 0.1
                break
            end
            if count > 1000
                kk = kk-1;
                break
            end
            l0 = l1;
            Thetas = Thetas1;
        end
        if count <=1000 && l0 > like
            like = l0;
            p_c = p;
            thetas_c = Thetas;
        end
    end
    w = p_c;
    v = thetas_c;
end

