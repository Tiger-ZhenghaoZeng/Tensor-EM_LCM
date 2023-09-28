% Compare different clustering methods

M = 100;
Js = (3:10)*10;
Errors_EM = zeros(8, M);
Errors_kmedoids = zeros(8, M);
Errors_kmeans = zeros(8, M);
Errors_link = zeros(8, M);
Errors_spectral = zeros(8, M);
RDs_EM = zeros(8, M);
RDs_kmedoids = zeros(8, M);
RDs_kmeans = zeros(8, M);
RDs_link = zeros(8, M);
RDs_spectral = zeros(8, M);
rng(521);
L = 5;
parfor ii = 1:8
    J = Js(ii);
    N = 10*J;  
    d1 = floor(J/3);
    d2 = floor(J/3);
    d3 = J-2*floor(J/3);
    theta = zeros(J,L);
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
        % Tensor-EM
        [w,X1,X2,X3] =  multi(sample,L,[d1,d2,d3],10,10);
        v = [X1;X2;X3];
        [~,~,z_em] = CEM(sample,v,w);
        Errors_EM(ii, mm) = evaluate_error(z,z_em,L);
        [~,RDs_EM(ii, mm),~] = AccMeasure(z, z_em);
        % K-medoids
        [z_kmedoids, C_kmedoids] = kmedoids(sample,L,'distance','hamming','replicates',10);
        Errors_kmedoids(ii, mm) = evaluate_error(z, z_kmedoids, L);
        [~,RDs_kmedoids(ii, mm),~] = AccMeasure(z, z_kmedoids);
        % K-means
        [z_kmeans,C_kmeans] = kmeans(sample, L, 'Distance','hamming', 'replicates', 10);
        Errors_kmeans(ii, mm) = evaluate_error(z, z_kmeans, L);
        [~,RDs_kmeans(ii, mm),~] = AccMeasure(z, z_kmeans); 
        % Max linkage
        z_link = clusterdata(sample,'Distance','hamming','Linkage','complete','MaxClust',L)';
        if length(categories(categorical(z_link)))== L
            Errors_link(ii, mm) = evaluate_error(z, z_link, L);
            [~,RDs_link(ii, mm),~] = AccMeasure(z, z_link);
        else
            Errors_link(ii, mm) = -1;
            RDs_link(ii, mm) = -1;
        end
        % Spectral Clustering
        z_spectral = spectral_cluster(sample, L);
        Errors_spectral(ii, mm) = evaluate_error(z, z_spectral, L);
        [~,RDs_spectral(ii, mm),~] = AccMeasure(z, z_spectral);
    end
end

Error_EM = mean(Errors_EM, 2);
Error_kmedoids = mean(Errors_kmedoids, 2);
Error_kmeans = mean(Errors_kmeans, 2);
Error_link = mean(Errors_link, 2);
Error_spectral = mean(Errors_spectral, 2);
RD_EM = mean(RDs_EM, 2);
RD_kmedoids = mean(RDs_kmedoids, 2);
RD_kmeans = mean(RDs_kmeans, 2);
RD_link = mean(RDs_link, 2);
RD_spectral = mean(RDs_spectral, 2);


% item parameters in {0.2, 0.4, 0.6, 0.8}
p1 = plot(Js, Error_link,'-b+','LineWidth',1,'MarkerSize',8);hold on;
p2 = plot(Js, Error_kmedoids,'-mx','LineWidth',1,'MarkerSize',8);hold on;
p3 = plot(Js, Error_kmeans,'-bo','LineWidth',1,'MarkerSize',6);hold on; 
p4 = plot(Js, Error_EM,'-r^','LineWidth',1,'MarkerSize',6);hold on; 
p5 = plot(Js, Error_spectral,'-b*','LineWidth',1,'MarkerSize',6);hold on; 
set(p3,'Color',[1 1/2 0])
set(p5,'Color',[0.8 1/3 0])
str = legend([p1 p2 p3 p4 p5],{'Max-Linkage','K-medoids','K-means', 'Tensor-EM', 'Spectral'},'FontSize',20);
set(str,'Interpreter','latex','Location','northeast')
xlabel('Number of Items J','FontSize',25);
ylabel('Fraction of errors','FontSize',25);
set(gca,'XLim',[25,105]); % X-axis range
set(gca,'YLim',[-0.015,0.6]); % Y-axis range
set(gca,'YTick',[0,0.1,0.2,0.3,0.4,0.5, 0.6, 0.7]); % Location of Ylabel Ticks
set(gca,'YTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7'},'FontSize',22); % Label of Those ticks
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% item parameters in {0.1, 0.2, 0.8, 0.9}
p1 = plot(Js, Error_link,'-b+','LineWidth',1,'MarkerSize',8);hold on;
p2 = plot(Js, Error_kmedoids,'-mx','LineWidth',1,'MarkerSize',8);hold on;
p3 = plot(Js, Error_kmeans,'-bo','LineWidth',1,'MarkerSize',6);hold on; 
p4 = plot(Js, Error_EM,'-r^','LineWidth',1,'MarkerSize',6);hold on; 
p5 = plot(Js, Error_spectral,'-b*','LineWidth',1,'MarkerSize',6);hold on; 
set(p3,'Color',[1 1/2 0])
set(p5,'Color',[0.8 1/3 0])
str = legend([p1 p2 p3 p4 p5],{'Max-Linkage','K-medoids','K-means', 'Tensor-EM', 'Spectral'},'FontSize',20);
set(str,'Interpreter','latex','Location','northeast')
xlabel('Number of Items J','FontSize',25);
ylabel('Fraction of errors','FontSize',25);
set(gca,'XLim',[25,105]); % X-axis range
set(gca,'YLim',[-0.015,0.06]); % Y-axis range
set(gca,'YTick',[0,0.005,0.01,0.02, 0.03, 0.04, 0.05]); % Location of Ylabel Ticks
set(gca,'YTickLabel',{'0','0.005','0.01','0.02','0.03','0.04','0.05'},'FontSize',22); % Label of Those ticks
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


