function FigChar(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW,LineStyle)

%   FIGCHAR(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW,LineStyle)
%
%   FigChar sets parameters of the figure (FigHandle) with NICHOLS, PZMAP, STEP, etc.
%       - size of axes (FigAxesSize),
%       - FigXLabel, FigYLabel and their sizes (FigAxesSize),
%       - FigTitle and its size (FigTitleSize),
%       - line width of the plots (PlotWidth),
%       - if you want black and white plots, set BaW=1 and LineStyle
%         default LineStyle = {'-','-.','--',':','.'}.         
%
%   See also FIGBODE, FIGNYQUIST, FIGPLOT, FIGSURF.
%
%   Author(s): J. Roubal (roubalj@control.felk.cvut.cz) [04-01-2007]



if nargin < 7
    disp('usage: FigChar(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth)')
    return;
elseif nargin < 8
    BaW = 0; LineStyle = {'-','-.','--',':','.'};
elseif nargin < 9
    LineStyle = {'-','-.','--',':','.'};
end



CurrentAxesHandle = get(FigHandle,'CurrentAxes');

set(CurrentAxesHandle,'FontSize',FigAxesSize,'FontWeight','bold');

XLabelHandle = xlabel(''); YLabelHandle = ylabel(''); TitleHandle = title('');
set(XLabelHandle,'FontSize',FigAxesSize,'FontWeight','bold','String',FigXLabel);
set(YLabelHandle,'FontSize',FigAxesSize,'FontWeight','bold','String',FigYLabel);
set(TitleHandle,'FontSize',FigTitleSize,'FontWeight','bold','String',FigTitle);

PlotsHandle = get(CurrentAxesHandle,'Children');
for i = 1 : 1 : size(PlotsHandle,1)-1
	switch get(PlotsHandle(i),'Type')
        case 'hggroup'
            PlotsHandleChildren = get(PlotsHandle(i),'Children');
            for j = 1 : 1 : size(PlotsHandleChildren,1)
                if BaW == 1, set(PlotsHandleChildren(j),'LineWidth',PlotWidth,'Color','k','LineStyle',LineStyle{size(PlotsHandle,1)-i});
                else set(PlotsHandleChildren(j),'LineWidth',PlotWidth);
                end
            end
    end
end


LegendHandle = legend;
set(LegendHandle,'FontSize',FigAxesSize,'FontWeight','bold');