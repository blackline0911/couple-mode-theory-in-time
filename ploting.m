function ploting(x,y,title_name,fontsize)
    % create a figure and plot data
    figure('color','w')
    for i=1:length(y)
        plot(x, y{i},'Linewidth',2)
        hold on
    end
    hold off
%     xlabel(x_label);
    set(gca,'fontsize',fontsize);
    title(title_name);
end