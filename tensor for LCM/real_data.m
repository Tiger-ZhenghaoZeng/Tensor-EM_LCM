% Read the five datasets and apply tensor-EM 

%timss1 = table2array(timss1);
class = 2:10;
d1 = 15;
d2 = 16;
d3 = 16;
J = 47;
GIC1_1 = zeros(1,length(class));
GIC2_1 = zeros(1,length(class));
rng(520);
for i = 1:length(class)
    [w,v1,v2,v3] = multi(timss1,class(i),[d1,d2,d3],10,100); 
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
    [w,v] = EM2(timss1,v,w/sum(w));
    GIC1_1(i) = GIC(w,v,timss1,1);
    GIC2_1(i) = GIC(w,v,timss1,2);
end

%timss2 = table2array(timss2);
GIC1_2 = zeros(1,length(class));
GIC2_2 = zeros(1,length(class));
for i = 1:length(class)
    [w,v1,v2,v3] = multi(timss2,class(i),[d1,d2,d3],10,100); 
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
    [w,v] = EM2(timss2,v,w/sum(w));
    GIC1_2(i) = GIC(w,v,timss2,1);
    GIC2_2(i) = GIC(w,v,timss2,2);
end

%timss3 = table2array(timss3);
GIC1_3 = zeros(1,length(class));
GIC2_3 = zeros(1,length(class));
for i = 1:length(class)
    [w,v1,v2,v3] = multi(timss3,class(i),[d1,d2,d3],10,100); 
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
    [w,v] = EM2(timss3,v,w/sum(w));
    GIC1_3(i) = GIC(w,v,timss3,1);
    GIC2_3(i) = GIC(w,v,timss3,2);
end

%timss4 = table2array(timss4);
GIC1_4 = zeros(1,length(class));
GIC2_4 = zeros(1,length(class));
for i = 1:length(class)
    [w,v1,v2,v3] = multi(timss4,class(i),[d1,d2,d3],10,100); 
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
    [w,v] = EM2(timss4,v,w/sum(w));
    GIC1_4(i) = GIC(w,v,timss4,1);
    GIC2_4(i) = GIC(w,v,timss4,2);
end

%timss5 = table2array(timss5);
GIC1_5 = zeros(1,length(class));
GIC2_5 = zeros(1,length(class));
for i = 1:length(class)
    [w,v1,v2,v3] = multi(timss5,class(i),[d1,d2,d3],10,100); 
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
    [w,v] = EM2(timss5,v,w/sum(w));
    GIC1_5(i) = GIC(w,v,timss5,1);
    GIC2_5(i) = GIC(w,v,timss5,2);
end

GIC1_pool = (GIC1_1+GIC1_2+GIC1_3+GIC1_4+GIC1_5)./5;
GIC2_pool = (GIC2_1+GIC2_2+GIC2_3+GIC2_4+GIC2_5)./5;

% Generate the GIC plot
figure
plot(class,GIC1_pool,'-ro','MarkerSize',10,'LineWidth',2);hold on;
plot(class,GIC2_pool,'-b+','MarkerSize',10,'LineWidth',2);hold on;
lgd=legend({'GIC1','GIC2'});
set(gca,'FontSize',20);
title('GIC of TIMSS data');
xlabel('Number of classes');
ylabel('GIC value');
lgd.Location = 'northwest';
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


% Estimation L=3
L = 3;
rng(520);
[w,v1,v2,v3] = multi(timss1,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p1,theta1] = EM2(timss1,v,w/sum(w));

[w,v1,v2,v3] = multi(timss2,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p2,theta2] = EM2(timss2,v,w/sum(w));


[w,v1,v2,v3] = multi(timss3,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p3,theta3] = EM2(timss3,v,w/sum(w));

[w,v1,v2,v3] = multi(timss4,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p4,theta4] = EM2(timss4,v,w/sum(w));

[w,v1,v2,v3] = multi(timss5,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p5,theta5] = EM2(timss5,v,w/sum(w));


% Match 2,3,4,5 with 1
[~,p2,~,theta2] = Match(p1,p2,theta1,theta2);
[~,p3,~,theta3] = Match(p1,p3,theta1,theta3);
[~,p4,~,theta4] = Match(p1,p4,theta1,theta4);
[~,p5,~,theta5] = Match(p1,p5,theta1,theta5);

% Average pooling
p3_hat = (p1+p2+p3+p4+p5)./5;
theta3_hat = (theta1 + theta2 + theta3 + theta4 + theta5)./5;

% Visualize the matrix theta
imagesc(theta3_hat)
xlabel('Different Class');
ylabel('Items','Rotation',90);
colorbar;
set(gca,'FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Estimation L=6
rng(521)
L = 6;
[w,v1,v2,v3] = multi(timss1,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p1,theta1] = EM2(timss1,v,w/sum(w));

[w,v1,v2,v3] = multi(timss2,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p2,theta2] = EM2(timss2,v,w/sum(w));


[w,v1,v2,v3] = multi(timss3,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p3,theta3] = EM2(timss3,v,w/sum(w));

[w,v1,v2,v3] = multi(timss4,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p4,theta4] = EM2(timss4,v,w/sum(w));

[w,v1,v2,v3] = multi(timss5,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p5,theta5] = EM2(timss5,v,w/sum(w));


% Match 2,3,4,5 with 1
[~,p2,~,theta2] = Match(p1,p2,theta1,theta2);
[~,p3,~,theta3] = Match(p1,p3,theta1,theta3);
[~,p4,~,theta4] = Match(p1,p4,theta1,theta4);
[~,p5,~,theta5] = Match(p1,p5,theta1,theta5);

p_hat = (p1+p2+p3+p4+p5)./5;
theta_hat = (theta1 + theta2 + theta3 + theta4 + theta5)./5;

% Estimation L=5
rng(520)
L = 5;
[w,v1,v2,v3] = multi(timss1,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p1,theta1] = EM2(timss1,v,w/sum(w));

[w,v1,v2,v3] = multi(timss2,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p2,theta2] = EM2(timss2,v,w/sum(w));


[w,v1,v2,v3] = multi(timss3,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p3,theta3] = EM2(timss3,v,w/sum(w));

[w,v1,v2,v3] = multi(timss4,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p4,theta4] = EM2(timss4,v,w/sum(w));

[w,v1,v2,v3] = multi(timss5,L,[d1,d2,d3],10,100); 
v = [v1;v2;v3];
for k = 1:J
    for j = 1:L
        if v(k,j)<0
            v(k,j) = 0.001;
        elseif v(k,j) > 1
            v(k,j) = 0.999;
        end
    end
end        
[p5,theta5] = EM2(timss5,v,w/sum(w));


% Match 2,3,4,5 with 1
[~,p2,~,theta2] = Match(p1,p2,theta1,theta2);
[~,p3,~,theta3] = Match(p1,p3,theta1,theta3);
[~,p4,~,theta4] = Match(p1,p4,theta1,theta4);
[~,p5,~,theta5] = Match(p1,p5,theta1,theta5);

% Average pooling
p5_hat = (p1+p2+p3+p4+p5)./5;
theta5_hat = (theta1 + theta2 + theta3 + theta4 + theta5)./5;


% Visualize the matrix
imagesc(theta5_hat)
xlabel('Different Class');
ylabel('Items','Rotation',90);
colorbar;
set(gca,'FontSize',15);
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% Partial order
gamma5 = partial_order(theta5_hat, 0.25);


