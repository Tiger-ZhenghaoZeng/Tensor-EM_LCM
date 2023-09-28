function [w,v,z_hat] = CEM(sample,Theta,L)
% Latent class classification EM algorithm with specific initial value
% Input: sample is a N*J dataset. Theta the initial values for item
% parameters. L the initial values for proportion vectors.
% Output: w the proportion vectors, v the item parameters, z_hat the
% estimated latent class membership
    R = sample';
    [J,G] = size(Theta);
    [~,N] = size(R);
    p = L;
    Thetas = Theta;
    Phi = zeros(N,G);
    Z = zeros(N,G);
    z = zeros(1,N);
    Thetas1 = zeros(J,G);
    l0 = -Inf;
    while 1 
        for i = 1:N
            de = 0;
            for m = 1:G
                de = de + p(m)*prod(Thetas(:,m).^(R(:,i)))*prod((1-Thetas(:,m)).^(1-R(:,i)));
            end
            for l = 1:G
                Phi(i,l) = p(l)*prod(Thetas(:,l).^(R(:,i)))*prod((1-Thetas(:,l)).^(1-R(:,i)))/de;
            end
            [~, index] = max(Phi(i, :));
            Z(i,:) = 0;
            Z(i, index) = 1;
            z(i) = index;
        end
        for l = 1:G
            p(l) = sum(Z(:,l))/N;
        end
        for j = 1:J
            for g = 1:G
                Thetas1(j,g) = sum(Z(:,g).*(R(j,:)'))/(sum(Z(:,g)));
            end
        end
        l1 = CEM_likeli(sample,z,Thetas1);
        if abs(l1-l0) < 0.1
            break
        end
        l0 = l1;
        Thetas = Thetas1;
    end
    w = p;
    v = Thetas;
    z_hat = z;
end

