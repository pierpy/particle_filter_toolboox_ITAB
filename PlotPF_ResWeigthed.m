function PlotPF_ResWeigthed(t, p, ph, pp, wp, wpp, x, xh, xp, y, yh)
    nx=size(xh,1);
    % plot of the state vs estimated state by the particle filter vs particle paths
    close all
    for i=1:nx
        figure
        hold on;
        for j = 1:size(x,2)
            ww = wp(j,:);
            if j ~=1
                wwn = (ww - min(ww))/(max(ww)-min(ww));
                wwn = wwn + 5*ones(1,length(wwn));
            else
                wwn = ww;
            end
            for p = 1:length(ww)              
                hold on
                h1 = plot(t(j), xp(i,j,p),'o','Color',[0.4 0.4 0.4],'MarkerSize',wwn(p));
            end
        end
        h2 = plot(t,x(i,:),'r','LineWidth',2);
        h3 = plot(t,xh(i,:),'b','LineWidth',2);

        legend([h2 h3 h1(1)],'state','mean of estimated state','particle paths');
        title(['State',num2str(i), ' vs estimated state by the particle filter vs particle paths'],'FontSize',14);
        xlabel('Time (s)')
    end
    for i=1:npar
        figure
        hold on;
        for j = 1:size(x,2)
            ww = wpp(j,:);
            if j ~=1
                wwn = (ww - min(ww))/(max(ww)-min(ww));
                wwn = wwn + 5*ones(1,length(wwn));
            else
                wwn = ww;
            end
            for p = 1:length(ww)
                hold on
                h1 = plot(t(j), pp(i,j,p),'o','Color',[0.4 0.4 0.4],'MarkerSize',wwn(p));
            end
        end
        h2 = line([t(1), t(end)], [p(i), p(i)], 'Color', 'r', 'LineWidth', 2);
        h3 = plot(t, ph(i,:),'b','LineWidth',2);

        legend([h2 h3 h1(1)],'state','mean of estimated state','particle paths');
        title(['State',num2str(i), ' vs estimated state by the particle filter vs particle paths'],'FontSize',14);
        xlabel('Time (s)')
    end
    % plot of the observation vs filtered observation by the particle filter
    figure
    plot(t,y,'b', t, yh,'r');
    legend('observation','filtered observation');
    title('Observation vs filtered observation by the particle filter','FontSize',14);

end