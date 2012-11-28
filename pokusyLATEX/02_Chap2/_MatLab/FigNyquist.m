function FigNyquist(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW)

%   FIGNYQUIST(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW)
%
%   FigNyquist sets parameters of the figure (FigHandle) with NYQUIST, etc.
%       - size of axes (FigAxesSize),
%       - FigXLabel, FigYLabel and their sizes (FigAxesSize),
%       - FigTitle and its size (FigTitleSize),
%       - line width of the plots (PlotWidth),
%       - if you want black and white plots, set BaW=1
%         (if you want to change line style you must use nyquist(SYS,'--'), etc.)
%
%   See also FIGBODE, FIGCHAR, FIGPLOT, FIGSURF.
%
%   Author(s): J. Roubal (roubalj@control.felk.cvut.cz) [04-01-2007]



if nargin < 7
    disp('usage: MyNyquist(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth)')
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
for i = 1 : 1 : size(PlotsHandle,1)
    switch get(PlotsHandle(i),'Type')
        case 'hggroup'
            PlotsHandleChildren = get(PlotsHandle(i),'Children');
            for j = 1 : 1 : size(PlotsHandleChildren,1)
                switch get(PlotsHandleChildren(j),'Type')
                    case 'line'
                        if BaW == 1, set(PlotsHandleChildren(j),'LineWidth',PlotWidth,'Color',[0 0 0]);
                        else set(PlotsHandleChildren(j),'LineWidth',PlotWidth);
                        end
                    case 'patch'
                        if BaW == 1, set(PlotsHandleChildren(j),'EdgeColor',[0 0 0],'FaceColor',[0 0 0]); end
                end
            end
    end
end

LegendHandle = legend;
set(LegendHandle,'FontSize',FigAxesSize,'FontWeight','bold');