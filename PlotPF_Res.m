function PlotPF_Res(t,x,xh,xp,wp,y,yh)

    nx=size(xh,1);
    ny=size(yh,1);
    % plot of the state vs estimated state by the particle filter vs particle paths
    close all
    for i=1:nx
    figure
    hold on;
    h1 = plot(t,squeeze(xp(i,:,:)),'o','Color',[0.4 0.4 0.4],'MarkerSize',3);
    h2 = plot(t,x(i,:),'r','LineWidth',2);
    h3 = plot(t,xh(i,:),'b','LineWidth',2);
    legend([h2 h3 h1(1)],'State','Estimated State','Particle Paths');
    title(['State ',num2str(i)],'FontSize',14);
    xlabel('Time (s)')
    end
    % plot of the observation vs filtered observation by the particle filter
    
    for i=1:ny
    figure
    h1 = plot(t,y(i,:),'r','LineWidth',2);
    hold on
    h2 = plot(t,yh(i,:),'b','LineWidth',2);
    legend([h1 h2],'Observation','Estimated observation');
    title(['Observation ',num2str(i)],'FontSize',14);
    xlabel('Time (s)')
    end
end