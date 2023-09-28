function [w,v] = EM2(sample,Theta,p_init)
% Latent class EM algorithm with specific initial value
% Input: sample is a N*J dataset. Theta the initial values for item
% parameters. L the initial values for proportion vectors.
% Output: w the proportion vectors, v the item parameters
    R = sample';
    [J,L] = size(Theta);
    [~,N] = size(R);
    p = p_init;
    Thetas = Theta;
    Phi = zeros(N,L);
    Thetas1 = zeros(J,L);
    l0 = likelihood(sample,p_init,Theta);
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
        if abs(l1-l0) < 0.1 || count > 1000
            break
        end
        l0 = l1;
        Thetas = Thetas1;
    end
    w = p;
    v = Thetas;
end

