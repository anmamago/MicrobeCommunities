hold off;
%markers = ['x','+','o','*'];

data1 = load('ConnCompAverages_Noise_0_Perm_0.csv')./1577;
%data2 = load('ConnCompAverages_Noise_0_Perm_1_Tics_101.csv')./1577;
data3 = load('ConnCompAverages_Noise_1_Perm_0.csv')./1577;
data4 = load('ConnCompAverages_Noise_1_Perm_1_Tics_101.csv')./1577;
%data2t = load('ConnCompAverages_Noise_0_Perm_1_Tics_301.csv')./1577;
data3t = load('ConnCompAverages_Noise_1_Perm_0_Tics_301.csv')./1577;
data4t = load('ConnCompAverages_Noise_1_Perm_1_Tics_301.csv')./1577;
datastraight = load('ConnCompAverages_Noise_0_Perm_0_Tics_301.csv')./1577;


ThresholdTics = 101; %this should be one more than the value you want: 10 -> 11
ThreshRange = linspace(0,1,ThresholdTics);
%ThresholdTics301 = 301; %this should be one more than the value you want: 10 -> 11
ThreshRange301 = linspace(0,1,301);

%This section of plotting is to make sure that the legend is easily
%readable.

plot(ThreshRange(51),data4(51),'g^-','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
hold on;
%plot(ThreshRange(51),data2(51),'mv-','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(ThreshRange(51),data3(51),'rs-','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(ThreshRange(51),data1(51),'bo-','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)


plot(ThreshRange,data1,'b','LineWidth',3)
%plot(ThreshRange301,data2t,'g','LineWidth',3)
plot(ThreshRange,data3,'r','LineWidth',3)
plot(ThreshRange301,data4t,'g','LineWidth',3)

%This section plots multiple markers along the lines

a = spaced_data(data1, ThreshRange);
%b = spaced_data(data2, ThreshRange);
c = spaced_data(data3, ThreshRange);
d = spaced_data(data4, ThreshRange);

hold on;


ThresholdTics301 = 301; %this should be one more than the value you want: 10 -> 11
ThreshRange301 = linspace(0,1,301);

viz_dropoff = zeros(2,10);
counter = 1;
for i = 75:90
    viz_dropoff(1,counter) = ThreshRange301(i);
    %viz_dropoff(2,counter) = data2t(i);
    counter = counter + 1;
end;


%making things more dense when they aren't otherwise dense.

viz_dropoff2 = zeros(2,10);
counter = 1;
for i = 75:85
    viz_dropoff2(1,counter) = ThreshRange301(i);
    viz_dropoff2(2,counter) = data4t(i);
    counter = counter + 1;
end;

viz_dropoff3 = zeros(2,6);
counter = 1;
for i = 80:85
    viz_dropoff3(1,counter) = ThreshRange301(i);
    viz_dropoff3(2,counter) = data3t(i);
    counter = counter+1;
end;


plot(d(1,:),d(2,:),'g^','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
%plot(b(1,:),b(2,:),'mv','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(c(1,:),c(2,:),'rs','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(a(1,:),a(2,:),'bo','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
%plot(viz_dropoff(1,1:end-3),viz_dropoff(2,1:end-3),'g^','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
%plot(viz_dropoff(1,end-1:end),viz_dropoff(2,end-1:end),'g^','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(viz_dropoff2(1,:),viz_dropoff2(2,:),'g^','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)
plot(viz_dropoff3(1,:),viz_dropoff3(2,:),'rs','MarkerSize',10,'MarkerFaceColor',[1,1,1],'LineWidth',3)

hold on;
x = zeros(1,10);
for i = 1:10
    x(i) = 0.36;
end;
y=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0];
plot(x,y,'k--','LineWidth',4)

set(gca,'FontSize',16)
xlabel('Threshold value')
ylabel('% of OTUs in largest component')
legend('Added noise, with permutation','Added noise only','Unaltered: no noise, no permutation' ) 

