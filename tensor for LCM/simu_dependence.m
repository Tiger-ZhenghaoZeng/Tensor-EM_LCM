% Section 5.2 performance under local dependence
M = 100;
EMtrue_MSE1 = zeros(1,M);
EMtrue_MSE2 = zeros(1,M);
EMrand_MSE1 = zeros(1,M);
EMrand_MSE2 = zeros(1,M);
EMtrue_time = zeros(1,M);
EMrand_time = zeros(1,M);
EM_tensor_MSE1 = zeros(1,M);
EM_tensor_MSE2 = zeros(1,M);
EM_tensor_time = zeros(1,M);
J = 200;
K = 5;
d1 = 66;
d2 = 67;
d3 = 67;
d = [0,d1,d2,d3];
N = 1000;
theta = zeros(J,K);
rho=0.3;
cova_mat = zeros(J,J);
for i = 1:J
    for j = 1:J
        cova_mat(i,j) = rho^(abs(i-j));
    end
end
mu = zeros(1,J);
for i = 1:J
    for j = 1:K
        a = binornd(1,0.5);
        b = binornd(1,0.5);
        if b ==1
            theta(i,j) = a*0.1+(1-a)*0.2;
        else
            theta(i,j) = a*0.8+(1-a)*0.9;
        end
    end
end
while 1
    p = rand(1,K);
    p = p/sum(p);
    if all(p>0.08)
        break
    end
end
cdf = cumsum(p);
parfor mm = 1:M     
    alphas = rand(1,N);
    sample = zeros(J,N);
    for j = 1:N
        % Find class membership
        m = 1;
        while alphas(j)>cdf(m)
            m = m+1;
        end
        % The multivariate normal distribution to generate correlated
        % components
        mv_samples = mvnrnd(mu, cova_mat, 1);
        for i = 1:J
            if mv_samples(i) < norminv(theta(i,m),0,1)
                sample(i,j) = 1;
            end
        end
    end
    tic;
    [w1,v1] = EM2(sample',theta,p);
    [~,w1,~,v1] = Match(p,w1,theta,v1);
    EMtrue_time(mm) = toc; 
    EMtrue_MSE1(mm) = mean((w1-p).^2);
    e1 = (theta-v1).^2;
    EMtrue_MSE2(mm) = mean(e1(:));
    tic;
    [w3,v3] = EM_random(sample',K,5);
    [~,w3,~,v3] = Match(p,w3,theta,v3);
    EMrand_time(mm) = toc;
    EMrand_MSE1(mm) = mean((w3-p).^2);
    e3 = (theta-v3).^2;
    EMrand_MSE2(mm) = mean(e3(:));
    sample = sample';
    tic;
    [w,X1,X2,X3] =  multi(sample,K,[d1,d2,d3],10,10);
    v = [X1;X2;X3];
    [w4,v4] = EM2(sample,v,w);
    [~,w4,~,v4] = Match(p,w4,theta,v4);
    EM_tensor_time(mm) = toc; 
    EM_tensor_MSE1(mm) = mean((w4-p).^2);
    e4 = (theta-v4).^2;
    EM_tensor_MSE2(mm) = mean(e4(:));
end

% MSE plot
boxplot([EMtrue_MSE2',EMrand_MSE2',EM_tensor_MSE2'], 'Labels', {'EM True', 'EM Random', 'EM tensor'});
title('MSE of item parameter','FontSize',25);
ylabel('MSE','FontSize',25)
set(gca,'FontSize',14);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Time plot
boxplot([EMtrue_time',EMrand_time',EM_tensor_time'], 'Labels', {'EM True', 'EM Random','EM tensor'});
title('Running time of algorithms','FontSize',25);
ylabel('Seconds','FontSize',25)
set(gca,'FontSize',14);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

