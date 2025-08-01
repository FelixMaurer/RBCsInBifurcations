
box on
hold off

scaleFac = 2;
FontSize = 7;
ax = gca;
ax.FontSize = FontSize*scaleFac;
ax.FontName = 'cmss10';
ax.LineWidth = 1.4;
ax.Layer = 'top';

axis square
axis equal
xlim(axLims)
ylim(axLims)
xticks(linspace(axLims(1),axLims(2),3));
yticks(linspace(axLims(1),axLims(2),3));
ax.Units = 'centimeters';
ax.Position = [2 1.5 3.05 3.05]*scaleFac;
%ax.Position = [0.1400    0.16   0.7750    0.7750];
ax.Units = 'normalized';

pR = signrank(xData,yData,'tail','right');
pL = signrank(xData,yData,'tail','left');

signRankString = sprintf('p_y = %.3f',pL);
if (pR<0.001)
    signRankString = sprintf('p_y > 0.999');
end
if (pL<0.001)
    signRankString = sprintf('p_y < 0.001');
end
annotation('textbox',[ax.Position(1)+0.01 ax.Position(2)+0.55 0.5 0],'EdgeColor','none','String',signRankString,'FontSize',FontSize*scaleFac);

signRankString = sprintf('p_x = %.3f',pR);
if (pR<0.001)
    signRankString = sprintf('p_x < 0.001');
end
if (pL<0.001)
    signRankString = sprintf('p_x > 0.999');
end
annotation('textbox',[ax.Position(1)+0.22 ax.Position(2)+0.09 0.5 0],'EdgeColor','none','String',signRankString,'FontSize',FontSize*scaleFac);

fig.Position = [200,200,fig.Position(3),fig.Position(4)];

