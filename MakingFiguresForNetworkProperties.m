close all;
clear;
plotperm = 0;

%%%MODULARITY%%%

modularity1 = load('Method1_036thresh_2000reps_mod.csv');
mod_sigma1 = sqrt(var(modularity1));
modularity2 = load('Method2_026thresh_replace0_mod.csv');
mod_sigma2 = sqrt(var(modularity2));
modularity3 = load('ErdosRenyiMethod_mod.csv');
mod_sigma3 = sqrt(var(modularity3));
modularity4 = load('ChungLuMethod_mod.csv'); %Something is weird here.
mod_sigma4 = sqrt(var(modularity4));

y1 = pdfsmooth(modularity1');
plot(y1(:,1),y1(:,2),'b','LineWidth',3)
h = fill([y1(:,1)' y1(1,1)], [y1(:,2)' y1(1,2)],'b');
set(h,'facealpha',0.9)
hold on;
if plotperm
    mod2 = reshape(modularity2, 1, 2000);
    y2 = pdfsmooth(mod2');
    plot(y2(:,1),y2(:,2),'r','LineWidth',3)
end
y3 = pdfsmooth(modularity3');
plot(y3(:,1),y3(:,2),'k','LineWidth',3)
%mod4 = reshape(modularity4,1,2000);
y4 = pdfsmooth(modularity4');
plot(y4(:,1),y4(:,2),'g','LineWidth',3)
set(gca,'FontSize',24)
xlabel('Modularity')
ylabel('Frequency')
xlim([0 1])
legend('Noise-added data','Erdos-Renyi model','Chung-Lu model')




%%%Average path length%%%
apl1 = load('Method1_036thresh_2000reps_apl.csv');
apl_sigma1 = sqrt(var(apl1));
apl2 = load('Method2_026thresh_replace0_apl.csv');
apl_sigma2 = sqrt(var(apl2));
apl3 = load('ErdosRenyiMethod_apl.csv');
apl_sigma3 = sqrt(var(apl3));
apl4 = load('ChungLuMethod_apl.csv');
apl_sigma4 = sqrt(var(apl4));

figure();
z1 = pdfsmooth(apl1');
plot(z1(:,1),z1(:,2),'b','LineWidth',3);
h = fill([z1(:,1)' z1(1,1)], [z1(:,2)' z1(1,2)],'b');
set(h,'facealpha',0.9)
hold on;
if plotperm
    apl2a = reshape(apl2,1,2000);
    z2 = pdfsmooth(apl2a');
    plot(z2(:,1),z2(:,2),'r','Linewidth',3);
end
z3 = pdfsmooth(apl3');
plot(z3(:,1),z3(:,2),'k','LineWidth',3);
%apl4a = reshape(apl4,1,2000);
z4 = pdfsmooth(apl4');
plot(z4(:,1),z4(:,2),'g','LineWidth',3);
xlabel('Average path length')
set(gca,'FontSize',24)
ylabel('Frequency')
legend('Noise-added data','Erdos-Renyi model','Chung-Lu model')


%%%Clustering coefficient%%%
cc1 = load('Method1_036thresh_2000reps_cc.csv');
cc2 = load('Method2_026thresh_replace0_cc.csv');
cc3 = load('ErdosRenyiMethod_cc.csv');
cc4 = load('ChungLuMethod_cc.csv');
cc_sigma1 = sqrt(var(cc1));
cc_sigma2 = sqrt(var(cc2));
cc_sigma3 = sqrt(var(cc3));
cc_sigma4 = sqrt(var(cc4));


figure();
x1 = pdfsmooth(cc1');
plot(x1(:,1),x1(:,2),'b','LineWidth',2);
h=fill([x1(:,1)' x1(1,1)], [x1(:,2)' x1(1,2)],'b');
set(h,'facealpha',0.9)
hold on;
if plotperm
    cc2a = reshape(cc2,1,2000);
    x2 = pdfsmooth(cc2a');
    plot(x2(:,1),x2(:,2),'r','LineWidth',3);
end
x3 = pdfsmooth(cc3');
plot(x3(:,1),x3(:,2),'k','LineWidth',3);
x4 = pdfsmooth(cc4');
plot(x4(:,1),x4(:,2),'g','LineWidth',3);
xlabel('Clustering coefficient')
set(gca,'FontSize',24)
xlabel('Clustering coefficient')
ylabel('Frequency')
xlim([0 1])
legend('Noise-added data','Erdos-Renyi model','Chung-Lu model')


%%%Diameter%%%
diam1 = load('Method1_036thresh_2000reps_diams.csv');
diam2 = load('Method2_026thresh_replace0_diams.csv');
diam3 = load('ErdosRenyiMethod_diams.csv');
diam4 = load('ChungLuMethod_diams.csv');
diam_sigma1 = sqrt(var(diam1));
diam_sigma2 = sqrt(var(diam2));
diam_sigma3 = sqrt(var(diam3));
diam_sigma4 = sqrt(var(diam4));

lowest = 1;
highest = 32;
figure();
v1a = min(diam1);
v1b = max(diam1);
w1 = hist(diam1, (v1b-v1a+1));
s1 = sum(w1);
begin1 = zeros(1,v1a-lowest);
ending1 = zeros(1,highest-v1b);
w1 = [begin1 w1 ending1];
%bar(v1a:v1b, w1./s1,'b')
%xlabel('Diameter')
%hold on;

%There is something wrong with this data set?
if plotperm
    diam2 = reshape(diam2,1,2000);
    v2a = min(diam2);
    v2b = max(diam2);
    w2 = hist(diam2, (v2b-v2a+1));
    s2 = sum(w2);
    begin2 = zeros(1,v2a-lowest);
    ending2 = zeros(1,highest-v2b);
    w2 = [begin2 w2 ending2];
    %bar(v2a:v2b, w2./s2,'r')
end

v3a = min(diam3);
v3b = max(diam3);
w3 = hist(diam3, (v3b-v3a+1));
s3 = sum(w3);
begin3 = zeros(1,v3a-lowest);
ending3 = zeros(1,highest-v3b);
w3 = [begin3 w3 ending3];
%bar(v3a:v3b, w3./s3,'g')

%diam4 = reshape(diam4,1,2000);
v4a = min(diam4);
v4b = max(diam4);
w4 = hist(diam4, (v4b-v4a+1));
s4 = sum(w4);
begin4 = zeros(1,v4a-lowest);
ending4 = zeros(1,highest-v4b);
w4 = [begin4 w4 ending4];
%bar(v4a:v4b, w4./s4,'k')

if plotperm
    u = [(w1./s1)',(w2./s2)',(w3./s3)', (w4./s4)'];
else
    u = [(w1./s1)',(w3./s3)', (w4./s4)'];
end
bar(u,2)
xlabel('Diameter')
set(gca,'FontSize',24)
ylabel('Frequency')
xlim([3 13])
legend('Noise-added data','Erdos-Renyi model','Chung-Lu model')


%%%Average Degree%%%
%Something wrong happened here.
deg1 = load('Method1_036thresh_2000reps_avgdeg.csv');
%diam2 = load('Method2_026thresh_replace0_diams.csv');
deg3 = load('ErdosRenyiMethod_deg.csv');
deg4 = load('ChungLuMethod_avgdeg.csv');
deg5 = load('ChungLuMethod_avgdeg_withsingletons.csv');
deg_sigma1 = sqrt(var(deg1));
deg_sigma3 = sqrt(var(deg3));
deg_sigma4 = sqrt(var(deg4));
deg_sigma4_w_singles = sqrt(var(deg5));

figure();
u1 = pdfsmooth(deg1');
plot(u1(:,1),u1(:,2),'b','LineWidth',2);
h=fill([u1(:,1)' u1(1,1)], [u1(:,2)' u1(1,2)],'b');
set(h,'facealpha',0.9)
hold on;
u3 = pdfsmooth(deg3');
plot(u3(:,1),u3(:,2),'k','LineWidth',2);
u4 = pdfsmooth(deg4');
plot(u4(:,1),u4(:,2),'g','LineWidth',2);
u5 = pdfsmooth(deg5');
plot(u5(:,1),u5(:,2),'c','LineWidth',2);
set(gca,'FontSize',24)
xlabel('Average degree')
ylabel('Frequency')
legend('Noise-added data','Erdos-Renyi model','Chung-Lu model','Chung-Lu with singletons')
