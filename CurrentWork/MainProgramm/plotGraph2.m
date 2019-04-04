function [] = plotGraph2(analytical, numerical, xFromTo, xmsg, ymsg,legend1, legend2, dir, titleTxt)
    figure('visible','off');
    hold on;
    plot(analytical,'-')
    plot(numerical,'--')
    legend({legend1, legend2});
    xlim([1 xFromTo])
    xlabel(xmsg)
    ylabel(ymsg)
    title(titleTxt)
    hold off;
    saveas(gcf,[pwd '/' dir '/' ymsg '.png'],'png')
end