function [] = plotGraph(numerical, xFromTo, xmsg, ymsg, legendMsg, dir)
    figure('visible','off');
    hold on;
    plot(numerical,'--')
    legend({legendMsg});
    set(gca, 'XtickLabel', xFromTo);
    xlabel(xmsg)
    ylabel(ymsg)
    hold off;
    saveas(gcf,[pwd '/' dir '/' ymsg '.png'],'png')
end