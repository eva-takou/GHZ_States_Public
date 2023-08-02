function fig_props(fntsize,xlab,ylab)


set(gca,'fontname','microsoft sans serif','fontsize',fntsize)
set(gcf,'color','w')

xlabel(xlab)
ylabel(ylab,'interpreter','latex')


set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  ...
    'LineWidth', 0.5)

end