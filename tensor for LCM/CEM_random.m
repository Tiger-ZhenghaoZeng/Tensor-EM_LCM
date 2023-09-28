function [w,v,z_hat] = CEM_random(sample,L,n)
% Latent class classification CEM algorithm with n random initials
% Input: sample is a N*J dataset. L the number of latent class, n the
% number of initial points
% Output: w the proportion vectors, v the item parameters, z_hat estimated
% latent class membership
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
        Z = zeros(N,L);
        z = zeros(1,N);
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
                [~, index] = max(Phi(i, :));
                Z(i,:) = 0;
                Z(i, index) = 1;
                z(i) = index;
            end
            for l = 1:L
               p(l) = sum(Z(:,l))/N;
            end
            for j = 1:J
                for g = 1:L
                    Thetas1(j,g) = sum(Z(:,g).*(R(j,:)'))/(sum(Z(:,g)));
                end
            end
            l1 = CEM_likeli(sample,z,Thetas1);
            if abs(l1-l0) < 0.1 || count > 1000
                break
            end
%             if count > 1000
%                 kk = kk-1;
%                 break
%             end
            l0 = l1;
            Thetas = Thetas1;
        end
        if l0 > like
            like = l0;
            p_c = p;
            thetas_c = Thetas;
        end
    end
    w = p_c;
    v = thetas_c;
    z_hat = z;
end

