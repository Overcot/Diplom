function [] = plotGraph2(analytical, numerical, xFromTo, xmsg, ymsg, dir, legend1, legend2)
    global precision
    figure('visible','off');
    hold on;
    plot(analytical,'-')
    plot(numerical,'--')
    legend({legend1, legend2});
    set(gca, 'XtickLabel', xFromTo);
    xlabel(xmsg)
    ylabel(ymsg)
    title(['our model precision ', num2str(precision)])
    hold off;
    saveas(gcf,[pwd '/' dir '/' ymsg '.png'],'png')
    saveas(gcf,[pwd '/' dir '/' ymsg '.fig'],'fig')
end