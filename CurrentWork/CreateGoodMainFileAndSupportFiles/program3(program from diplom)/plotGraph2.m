function [] = plotGraph2(analytical, numerical, xFromTo, xmsg, ymsg, dir)
    figure('visible','off');
    hold on;
    plot(analytical,'-')
    plot(numerical,'--')
    legend({'Our Model', 'Tahvonen'});
    set(gca, 'XtickLabel', xFromTo);
    xlabel(xmsg)
    ylabel(ymsg)
    hold off;
    saveas(gcf,[pwd '/' dir '/' ymsg '.png'],'png')
    saveas(gcf,[pwd '/' dir '/' ymsg '.fig'],'fig')
end