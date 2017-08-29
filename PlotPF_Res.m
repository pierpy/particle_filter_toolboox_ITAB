function PlotPF_Res(t,x,xh,xp,wp,y,yh)





nx=size(xh,1);
%% plot of the state vs estimated state by the particle filter vs particle paths
close all
for i=1:nx
figure
hold on;
h1 = plot(t,squeeze(xp(i,:,:)),'o','Color',[0.4 0.4 0.4],'MarkerSize',1);
h2 = plot(t,x(i,:),'r','LineWidth',3);
h3 = plot(t,xh(i,:),'b','LineWidth',3);
legend([h2 h3 h1(1)],'state','mean of estimated state','particle paths');
title(['State',num2str(i), ' vs estimated state by the particle filter vs particle paths'],'FontSize',14);
xlabel('Time (s)')
end
%% plot of the observation vs filtered observation by the particle filter
figure
plot(t,y,'b', t,yh,'r');
legend('observation','filtered observation');
title('Observation vs filtered observation by the particle filter','FontSize',14);
