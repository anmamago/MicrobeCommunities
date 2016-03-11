close all;
clear;
threshold_vals = (0.025:0.025:0.975);

noise = 0;
ploton = 0;
en = length(threshold_vals)
avg_nodes_mx_between = zeros(1,en);
avg_edges_mx_between = zeros(1,en);
avg_mod_mx_between = zeros(1,en);
avg_cc_mx_between = zeros(1,en);
avg_diam_mx_between = zeros(1,en);
avg_deg_mx_between = zeros(1,en);
avg_path_mx_between = zeros(1,en);
%avg_c_mx_between = zeros(1,en);

for i = 1:en
    i
    x = threshold_vals(i);
    [ n_nodes, n_edges, modular, cc, avgdeg, apl, diams ] = BarberanData_v16_thresholdrange_nopval_cval_noise(x,noise); 
    avg_nodes_mx_between(i) = mean(n_nodes);
    avg_edges_mx_between(i) = mean(n_edges);
    avg_mod_mx_between(i) = mean(modular);
    avg_cc_mx_between(i) = mean(cc);
    avg_diam_mx_between(i) = mean(diams);
    avg_deg_mx_between(i) = mean(avgdeg);
    avg_path_mx_between(i) = mean(apl);
    %avg_c_mx_between(i) = mean(cval);
end;

%For average degree:
degfile = sprintf('ThreshRangeAvgDegForPlot_noise_%d.csv', noise);
%csvwrite(degfile, avg_deg_mx_between)
if ploton == 1
    figure();
    avg_deg_mx_between_noise_0 = load('ThreshRangeAvgDegForPlot_noise_0.csv');
    semilogy(threshold_vals,avg_deg_mx_between_noise_0,'bs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
    avg_deg_mx_bw_noise_1 = load('ThreshRangeAvgDegForPlot.csv');
    semilogy(threshold_vals,avg_deg_mx_bw_noise_1, 'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average degree')
    linex = zeros(1,8000);
    for i = 1:8000
        linex(i) = 0.36;
    end;
    liney=linspace(0,8000,8000);

    plot(linex,liney,'k--','LineWidth',5)
    xlim([0 1])
    legend('Unaltered', 'Added noise')
end;

%For diameter:
diamfile = sprintf('ThreshRangeAvgDiamForPlot_noise_%d.csv', noise);
%csvwrite(diamfile, avg_diam_mx_between)
if ploton ==1
    figure();
    avg_diam_mx_between_noise_0 = load('ThreshRangeAvgDiamForPlot_noise_0.csv');
    plot(threshold_vals,avg_diam_mx_between_noise_0,'bs-','MarkerFace',[1 1 1], 'LineWidth',2)
    hold on;
    avg_diam_mx_between_noise_1 = load('ThreshRangeAvgDiamForPlot.csv');
    plot(threshold_vals,avg_diam_mx_between_noise_1,'rs-','MarkerFace',[1 1 1], 'LineWidth',2)
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average diameter')
    linex = zeros(1,100);
    for i=1:100
        linex(i)=0.36;
    end
    liney=linspace(0,100,100);
    plot(linex,liney,'k--','LineWidth',5)
    ylim([0 20])
end

%For average path length
aplfile = sprintf('ThreshRangeAplForPlot_noise_%d.csv', noise);
%csvwrite(aplfile, avg_path_mx_between)
if ploton==1
    figure();
    avg_path_mx_between_noise_0 = load('ThreshRangeAplForPlot_noise_0.csv');
    plot(threshold_vals,avg_path_mx_between_noise_0,'bs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
    avg_path_mx_between_noise_1 = load('ThreshRangeAvgPathForPlot.csv');
    plot(threshold_vals,avg_path_mx_between_noise_1,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average path length')
    linex=zeros(1,100);
    for i=1:100
        linex(i)=0.36;
    end
    liney=linspace(0,100,100);
    plot(linex,liney,'k--','LineWidth',5)
    xlim([0 1])
    ylim([1 7])
end

%For modularity
modfile = sprintf('ThreshRangeAvgModForPlot_noise_%d.csv', noise);
%csvwrite(modfile, avg_mod_mx_between)
if ploton==1
    figure();
    avg_mod_mx_between_noise_0 = load('ThreshRangeAvgModForPlot_noise_0.csv');
    plot(threshold_vals,avg_mod_mx_between_noise_0,'bs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
    avg_mod_mx_between_noise_1 = load('ThreshRangeAvgModForPlot_noise_1.csv');
    plot(threshold_vals,avg_mod_mx_between_noise_1,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average modularity')
    linex=zeros(1,100);
    for i=1:100
        linex(i)=0.36;
    end
    liney=linspace(0,100,100);
    plot(linex,liney,'k--','LineWidth',5)
    xlim([0 1])
    ylim([0 1])
end

% For Clustering Coefficient
ccfile = sprintf('ThreshRangeAvgCCForPlot_noise_%d.csv', noise);
%csvwrite(ccfile, avg_cc_mx_between)
if ploton ==1
    figure();
    avg_cc_mx_between_noise_0 = load('ThreshRangeAvgCCForPlot_noise_0.csv');
    plot(threshold_vals,avg_cc_mx_between_noise_0,'bs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
    avg_cc_mx_between_noise_1 = load('ThreshRangeAvgCCForPlot.csv');
    plot(threshold_vals,avg_cc_mx_between_noise_1,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average clustering coeff.')
    linex=zeros(1,100);
    for i=1:100
        linex(i)=0.36;
    end
    liney=linspace(0,100,100);
    plot(linex,liney,'k--','LineWidth',5)
    ylim([0 1])
    xlim([0 1])
end;

%csvwrite('Fig2_avgnodesmx', avg_nodes_mx_between);
%csvwrite('Fig2_avgedgesmx', avg_edges_mx_between);
