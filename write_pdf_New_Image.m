function write_pdf_New_Image(pdfname,fig,resw,resh)
            
fig.Position = [0 0 resw resh];

ax = gca;
% outerpos = ax.OuterPosition;
% ti = 0;
% left = outerpos(1);% + ti(1);
% bottom = outerpos(2);% + ti(2);
% ax_width = outerpos(3);% - ti(1) - ti(3);
% ax_height = outerpos(4);% - ti(2) - ti(4)-0.04;
% ax.Position = [left bottom ax_width ax_height];

set(ax,'Color','w');

fig.PaperPositionMode = 'manual';
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 resw resh];
fig.PaperSize = [resw resh]; 

% fig.InvertHardcopy = 'off';
curpath = cd;
print(fig,'-dpdf',[curpath ,'\', pdfname '.pdf'],'-r0');