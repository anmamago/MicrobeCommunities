close all;
clear;
%threshold_vals = [ 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975];
threshold_vals = (0.025:0.025:0.975);
%threshold_vals = [0.4, 0.45];
doingdegrees=0;
doingdiams = 0;
doingpath = 0;
doingmod = 1;
noise = 0;
en = length(threshold_vals)
avg_nodes_mx_between = zeros(1,en);
avg_edges_mx_between = zeros(1,en);
avg_mod_mx_between = zeros(1,en);
%avg_cc_mx_between = zeros(1,en);
%avg_diam_mx_between = zeros(1,en);
%avg_deg_mx_between = zeros(1,en);
%avg_path_mx_between = zeros(1,en);
%avg_c_mx_between = zeros(1,en);

for i = 1:en
    i
    x = threshold_vals(i);
    if (noise == 1)
        [ n_nodes, n_edges, modular ] = BarberanData_v16_thresholdrange_nopval_cval(x);
    else
        [ n_nodes, n_edges, modular] = BarberanData_v16_thresholdrange_nopval_cval_noise(x,noise);
    end;
        %avg_nodes_mx_between(i) = mean(n_nodes);
    %avg_edges_mx_between(i) = mean(n_edges);
    avg_mod_mx_between(i) = mean(modular);
    %avg_cc_mx_between(i) = mean(cc);
    %avg_diam_mx_between(i) = mean(diams);
    %avg_deg_mx_between(i) = mean(avgdeg);
    %avg_path_mx_between(i) = mean(apl);
    %avg_c_mx_between(i) = mean(cval);
end;


if doingdegrees
    csvwrite('ThreshRangeAvgDegForPlot.csv',avg_deg_mx_between)
    semilogy(threshold_vals,avg_deg_mx_between,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
    set(gca,'FontSize',28)
    xlabel('Threshold')
    ylabel('Average degree')
    linex = zeros(1,1000);
    for i = 1:1000
        linex(i) = 0.36;
    end;
    liney=linspace(0,1000,1000);

    plot(linex,liney,'k--','LineWidth',5)
    xlim([0 1])
end;

if doingdiams
    csvwrite('ThreshRangeAvgDiamForPlot.csv',avg_diam_mx_between)
    plot(threshold_vals,avg_diam_mx_between,'rs-','MarkerFace',[1 1 1], 'LineWidth',2)
    hold on;
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

if doingpath
    csvwrite('ThreshRangeAvgPathForPlot.csv',avg_path_mx_between)
    plot(threshold_vals,avg_path_mx_between,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
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
    ylim([0 5])
end

if doingmod
    csvwrite('ThreshRangeAvgModForPlot.csv',avg_mod_mx_between);

    plot(threshold_vals,avg_mod_mx_between,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
    hold on;
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

% 
% csvwrite('ThreshRangeAvgCCForPlot.csv',avg_cc_mx_between)
% plot(threshold_vals,avg_cc_mx_between,'rs-','MarkerFace',[1 1 1],'LineWidth',2)
% hold on;
% set(gca,'FontSize',28)
% xlabel('Threshold')
% ylabel('Average clustering coeff.')
% linex=zeros(1,100);
% for i=1:100
%     linex(i)=0.36;
% end
% liney=linspace(0,100,100);
% plot(linex,liney,'k--','LineWidth',5)
% ylim([0 1])
% xlim([0 1])

