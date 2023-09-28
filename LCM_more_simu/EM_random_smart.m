function [w,v] = EM_random_smart(sample,L,K,n)
% Latent class EM algorithm with K random initials, refined tolerance
% Input: sample is a N*J dataset. L the number of latent class, K the
% number of initial points, n the number of iterations in the first stage
% Output: w the proportion vectors, v the item parameters
    R = sample';
    [J,N] = size(R);
    like = -Inf;
    p_c = ones(1,L)/L;
    thetas_c = zeros(J,L);
    % First K trial
    for kk = 1:K
        p = ones(1,L)/L;
        Thetas = rand(J,L);
        Phi = zeros(N,L);
        Thetas1 = zeros(J,L);
        for a = 1:n
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
        end
        l1 = likelihood(sample,p,Thetas1);
        if  l1 > like
            like = l1;
            p_c = p;
            thetas_c = Thetas;  
        end
    end
    
    % find the best one and do more iterations
    Thetas = thetas_c;
    p = p_c;
    Phi = zeros(N,L);
    Thetas1 = zeros(J,L);
    l0 = l1;
    for a = 1:1000
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
        l0 = l1;
        Thetas = Thetas1;
    end
    w = p;
    v = Thetas;
end

