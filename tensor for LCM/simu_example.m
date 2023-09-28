%Fixed-effect LCM replications N = 1000,[0.1,0.2,0.8,0.9],J=100, L=5
M = 100;
tensor_MSE1 = zeros(1,M);
tensor_MSE2 = zeros(1,M);
EMtrue_MSE1 = zeros(1,M);
EMtrue_MSE2 = zeros(1,M);
EMrand_MSE1 = zeros(1,M);
EMrand_MSE2 = zeros(1,M);
tensor_time = zeros(1,M);
EMtrue_time = zeros(1,M);
EMrand_time = zeros(1,M);
EM_tensor_MSE1 = zeros(1,M);
EM_tensor_MSE2 = zeros(1,M);
EM_tensor_time = zeros(1,M);
J = 100;
L = 5; 
d1 = 33;
d2 = 33;
d3 = 34;
N = 1000;
repeat = 5;
number = zeros(1,M);
theta = zeros(J,L);
p = ones(1,L)/L;
z = randi(L,1,N);
% Generate item parameters
for i = 1:J
    for j = 1:L
        a = binornd(1,0.5);
        b = binornd(1,0.5);
        if b ==1
            theta(i,j) = a*0.1+(1-a)*0.2;
        else
            theta(i,j) = a*0.8+(1-a)*0.9;
        end
    end
end
for mm = 1:M
    % Generate samples
    sample = zeros(N,J);
    for i = 1:N
        for j = 1:J
            sample(i,j) = binornd(1, theta(j, z(i)));
        end
    end    
    % true parameters as initial values
    tic;
    [w1,v1,z1] = CEM(sample,theta,ones(1,L)/L);
    [~,w1,~,v1] = Match(p,w1,theta,v1);
    EMtrue_time(mm) = toc;
    EMtrue_MSE1(mm) = mean((w1-p).^2);
    e1 = (theta-v1).^2;
    EMtrue_MSE2(mm) = mean(e1(:));
    % Random initial values
    tic;
    [w3,v3,z3] = CEM_random(sample,L,5);
    [~,w3,~,v3] = Match(p,w3,theta,v3);
    EMrand_time(mm) = toc;
    EMrand_MSE1(mm) = mean((w3-p).^2);
    e3 = (theta-v3).^2;
    EMrand_MSE2(mm) = mean(e3(:));
    % Tensor estimates
    tic;
    wc = 0;
    Xc = 0;
    for i = 1:repeat
        r = randperm(size(sample,2));
        samples = sample(:,r);
        [w,X1,X2,X3] =  multi(samples,L,[d1,d2,d3],10,10);
        v = [X1;X2;X3];
        X_hat = zeros(size(v));
        for j = 1:J
            X_hat(r(j),:)=v(j,:);
        end
        [~,w,~,X_hat] = Match(p,w,theta,X_hat);
        wc = wc + w;
        Xc = Xc + X_hat;
    end
    wc = wc/repeat;
    Xc = Xc/repeat;
    tensor_time(mm) = toc;
    tensor_MSE1(mm) = mean((wc-p).^2);
    e2 = (Xc-theta).^2;
    tensor_MSE2(mm) = mean(e2(:));
    % Tensor estimator as initial values
    tic;
    [w,X1,X2,X3] =  multi(sample,L,[d1,d2,d3],10,10);
    v = [X1;X2;X3];
    [w4,v4,z4] = CEM(sample,v,w);
    [~,w4,~,v4] = Match(p,w4,theta,v4);
    EM_tensor_time(mm) = toc;
    EM_tensor_MSE1(mm) = mean((w4-p).^2);
    e4 = (theta-v4).^2;
    EM_tensor_MSE2(mm) = mean(e4(:));
end

%replications N = 1000,[0.1,0.2,0.8,0.9],J=100, K=5
M = 100;
tensor_MSE1 = zeros(1,M);
tensor_MSE2 = zeros(1,M);
EMtrue_MSE1 = zeros(1,M);
EMtrue_MSE2 = zeros(1,M);
EMrand_MSE1 = zeros(1,M);
EMrand_MSE2 = zeros(1,M);
tensor_time = zeros(1,M);
EMtrue_time = zeros(1,M);
EMrand_time = zeros(1,M);
EM_tensor_MSE1 = zeros(1,M);
EM_tensor_MSE2 = zeros(1,M);
EM_tensor_time = zeros(1,M);
class = 2:7;
J = 100;
L = 5;
d1 = 33;
d2 = 33;
d3 = 34;
d = [0,d1,d2,d3];
N = 1000;
repeat = 5;
number1 = zeros(1,M);
number2 = zeros(1,M);
theta = zeros(J,L);
for i = 1:J
    for j = 1:L
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
    tic;
    [w1,v1] = EM2(sample',theta,p);
    [~,w1,~,v1] = Match(p,w1,theta,v1);
    EMtrue_time(mm) = toc; 
    EMtrue_MSE1(mm) = mean((w1-p).^2);
    e1 = (theta-v1).^2;
    EMtrue_MSE2(mm) = mean(e1(:));
    tic;
    [w3,v3] = EM_random(sample',L,5);
    [~,w3,~,v3] = Match(p,w3,theta,v3);
    EMrand_time(mm) = toc;
    EMrand_MSE1(mm) = mean((w3-p).^2);
    e3 = (theta-v3).^2;
    EMrand_MSE2(mm) = mean(e3(:));
    sample = sample';
    tic;
    wc = 0;
    Xc = 0;
    for i = 1:repeat
        r = randperm(size(sample,2));
        samples = sample(:,r);
        [w,X1,X2,X3] =  multi(samples,L,[d1, d2, d3],10,10);
        v = [X1;X2;X3];
        X_hat = zeros(size(v));
        for j = 1:J
            X_hat(r(j),:)=v(j,:);
        end
        [~,w,~,X_hat] = Match(p,w,theta,X_hat);
        wc = wc + w;
        Xc = Xc + X_hat;
    end
    wc = wc/repeat;
    Xc = Xc/repeat;
    tensor_time(mm) = toc;
    tensor_MSE1(mm) = mean((wc-p).^2);
    e2 = (Xc-theta).^2;
    tensor_MSE2(mm) = mean(e2(:));
    tic;
    [w,X1,X2,X3] =  multi(sample,L,[d1,d2,d3],10,10);
    v = [X1;X2;X3];
    [w4,v4] = EM2(sample,v,w);
    [~,w4,~,v4] = Match(p,w4,theta,v4);
    EM_tensor_time(mm) = toc; 
    EM_tensor_MSE1(mm) = mean((w4-p).^2);
    e4 = (theta-v4).^2;
    EM_tensor_MSE2(mm) = mean(e4(:));
    GIC1 = zeros(1,length(class));
    GIC2 = zeros(1,length(class));
    for i = 1:length(class)
        [w,v1,v2,v3] = multi(sample,class(i),[d1,d2,d3],20,20); 
        v = [v1;v2;v3];
        for k = 1:J
            for j = 1:class(i)
                if v(k,j)<0
                    v(k,j) = 0.001;
                elseif v(k,j) > 1
                    v(k,j) = 0.999;
                end
            end
        end        
        [w,v] = EM2(sample,v,w/sum(w));
        GIC1(i) = GIC(w,v,sample,1);
        GIC2(i) = GIC(w,v,sample,2);
    end
    number1(mm) = find(GIC1 == min(GIC1));
    number2(mm) = find(GIC2 == min(GIC2));
end

% MSE plot
boxplot([EMtrue_MSE2',EMrand_MSE2',tensor_MSE2',EM_tensor_MSE2'], 'Labels', {'EM True', 'EM Random', 'Tensor', 'EM tensor'});
title('MSE of item parameter','FontSize',25);
ylabel('MSE','FontSize',25)
set(gca,'FontSize',20);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% MSE plot without EM-random
boxplot([EMtrue_MSE2',tensor_MSE2',EM_tensor_MSE2'], 'Labels', {'EM True',  'Tensor', 'EM tensor'});
title('MSE of item parameter','FontSize',25);
ylabel('MSE','FontSize',25)
set(gca,'FontSize',20);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Running time plot
boxplot([EMtrue_time',EMrand_time',tensor_time',EM_tensor_time'], 'Labels', {'EM True', 'EM Random', 'Tensor','EM tensor'});
title('Running time of algorithms','FontSize',25);
ylabel('Seconds','FontSize',25)
set(gca,'FontSize',20);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

boxplot([EMtrue_time',tensor_time',EM_tensor_time'], 'Labels', {'EM True', 'Tensor','EM tensor'});
title('Running time of algorithms','FontSize',25);
ylabel('seconds','FontSize',25)
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


% Error Rate simulations
M = 100;
Js = (2:10)*10;
Errors = zeros(9, M);
rng(521);
parfor ii = 1:9
    J = Js(ii);
    N = 10*J;
    L = 5;
    d1 = floor(J/3);
    d2 = floor(J/3);
    d3 = J-2*floor(J/3);
    theta = zeros(J,L);
    % L = ones(1,K)/K;
    z = randi(L,1,N);
    for i = 1:J
        for j = 1:L
            a = binornd(1,0.5);
            b = binornd(1,0.5);
            if b ==1
                theta(i,j) = a*0.1+(1-a)*0.2;
            else
                theta(i,j) = a*0.8+(1-a)*0.9;
            end
        end
    end
    for mm = 1:M
        sample = zeros(N,J);
        for i = 1:N
            for j = 1:J
                sample(i,j) = binornd(1, theta(j, z(i)));
            end
        end
        [w,X1,X2,X3] =  multi(sample,L,[d1,d2,d3],10,10);
        v = [X1;X2;X3];
        [~,~,zhat] = CEM(sample,v,w);
        Errors(ii, mm) = evaluate_error(z,zhat,L);
    end
end

% Generate plots in Section 5.3
boxplot([Errors(2,:)',Errors(3,:)',Errors(4,:)',Errors(5,:)',Errors(6,:)',Errors(7,:)',Errors(8,:)',Errors(9,:)'], 'Labels', {'J=30','J=40','J=50','J=60','J=70','J=80','J=90','J=100'});
title('Error rate of joint MLE','FontSize',25);
ylabel('Fraction of errors','FontSize',25)
xlabel('Number of items J','FontSize',25)
set(gca,'FontSize',20);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'title','-dpdf')

