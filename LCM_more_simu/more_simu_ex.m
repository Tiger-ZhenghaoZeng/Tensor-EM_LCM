% Section 5.1 EM-random and tensor-EM from the same initial values
% EM-random with refined tolerance level

M = 100;
EMtrue_MSE1 = zeros(1,M);
EMtrue_MSE2 = zeros(1,M);
EMrand_MSE1 = zeros(1,M);
EMrand_MSE2 = zeros(1,M);
EMrandsmart_MSE1 = zeros(1,M);
EMrandsmart_MSE2 = zeros(1,M);
EMtrue_time = zeros(1,M);
EMrand_time = zeros(1,M);
EMrandsmart_time = zeros(1,M);
EM_tensor_MSE1 = zeros(1,M);
EM_tensor_MSE2 = zeros(1,M);
EM_tensor_time = zeros(1,M);
J = 100;
L = 5;
d1 = 33;
d2 = 33;
d3 = 34;
d = [0,d1,d2,d3];
N = 1000;
theta = zeros(J,L);
for i = 1:J
    for j = 1:L
        a = binornd(1,0.5);
        b = binornd(1,0.5);
        if b ==1
            theta(i,j) = a*0.2+(1-a)*0.4;
        else
            theta(i,j) = a*0.6+(1-a)*0.8;
        end
    end
end
while 1
    p = rand(1,L);
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
        m = 1;
        while alphas(j)>cdf(m)
            m = m+1;
        end
        for i = 1:J
            sample(i,j) = binornd(1,theta(i,m));
        end
    end
    % Generate theta_1
    theta1_init = rand(d1, L);
    index = cumsum(d);
    sample = sample';
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
    C12 = pairs23*MPinv(pairs13,L);
    C13 = (pairs23')*MPinv(pairs12,L);
    % Obtain initial values for theta2 and theta3
    theta2_init = C12*theta1_init;
    theta3_init = C13*theta1_init;
    for i = 1:d2
        for j = 1:L
            if theta2_init(i,j) <= 0
                theta2_init(i,j) = 0.001;
            elseif theta2_init(i,j) >= 1
                theta2_init(i,j) = 0.999;
            end
        end
    end
    for i = 1:d3
        for j = 1:L
            if theta3_init(i,j) <= 0
                theta3_init(i,j) = 0.001;
            elseif theta3_init(i,j) >= 1
                theta3_init(i,j) = 0.999;
            end
        end
    end
    tic;
    [w1,v1] = EM2(sample,[theta1_init; theta2_init; theta3_init],ones(1,L)/L);
    [~,w1,~,v1] = Match(p,w1,theta,v1);
    EMrand_time(mm) = toc; 
    EMrand_MSE1(mm) = mean((w1-p).^2);
    e1 = (theta-v1).^2;
    EMrand_MSE2(mm) = mean(e1(:));
    tic;
    [w,X1,X2,X3] = multi(sample,L,[d1,d2,d3],theta1_init, 20, p);
    v = [X1;X2;X3];
    [w2,v2] = EM2(sample,v,w);
    [~,w2,~,v2] = Match(p,w2,theta,v2);
    EM_tensor_time(mm) = toc; 
    EM_tensor_MSE1(mm) = mean((w2-p).^2);
    e2 = (theta-v2).^2;
    EM_tensor_MSE2(mm) = mean(e2(:)); 
    tic;
    [w3,v3] = EM_random_smart(sample,L,10,20);
    [~,w3,~,v3] = Match(p,w3,theta,v3);
    EMrandsmart_time(mm) = toc;
    EMrandsmart_MSE1(mm) = mean((w3-p).^2);
    e3 = (theta-v3).^2;
    EMrandsmart_MSE2(mm) = mean(e3(:));
    tic;
    [w4,v4] = EM2(sample,theta,p);
    [~,w4,~,v4] = Match(p,w4,theta,v4);
    EMtrue_time(mm) = toc; 
    EMtrue_MSE1(mm) = mean((w4-p).^2);
    e4 = (theta-v4).^2;
    EMtrue_MSE2(mm) = mean(e4(:)); 
end




% MSE plot
boxplot([EMtrue_MSE2',EMrand_MSE2',EMrandsmart_MSE2',EM_tensor_MSE2'], 'Labels', {'EM True', 'EM Random', 'EM Random(refined)', 'EM tensor'});
title('MSE of item parameter','FontSize',25);
ylabel('MSE','FontSize',25)
set(gca,'FontSize',14);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


% Time plot
boxplot([EMtrue_time',EMrand_time',EMrandsmart_time',EM_tensor_time'], 'Labels', {'EM True', 'EM Random', 'EM Random(refined)','EM tensor'});
title('Running time of algorithms','FontSize',25);
ylabel('Seconds','FontSize',25)
set(gca,'FontSize',14);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
