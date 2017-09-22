function [] = plotGraph(analytical, numerical, xFromTo, xmsg, ymsg)
    figure
    hold on;
    plot(analytical,'-')
    plot(numerical,'--')
    legend({'Analytical', 'Numerical'});
    set(gca, 'XtickLabel', xFromTo);
    xlabel(xmsg)
    ylabel(ymsg)
    hold off;
end