function [] = plotGraph(analytical, numerical, xFromTo, xmsg, ymsg, dir)
    figure('visible','off');
    hold on;
    plot(analytical,'-')
    plot(numerical,'--')
    legend({'Analytical', 'Numerical'});
    set(gca, 'XtickLabel', xFromTo);
    xlabel(xmsg)
    ylabel(ymsg)
    hold off;
    saveas(gcf,[pwd '/' dir '/' ymsg '.png'],'png')
    saveas(gcf,[pwd '/' dir '/' ymsg '.fig'],'fig')
end