function FigPlot(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW,LineStyle)

%   FigPlot(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW,LineStyle)
%
%   FigPlot sets parameters of the figure (FigHandle) with PLOT-graph
%       - size of axes (FigAxesSize),
%       - FigXLabel, FigYLabel and their sizes (FigAxesSize),
%       - FigTitle and its size (FigTitleSize),
%       - line width of the plots (PlotWidth),
%       - if you want black and white plots, set BaW=1 and LineStyle
%         default LineStyle = {'-','-.','--',':','.'}.         
%
%   See also FIGBODE, FIGCHAR, FIGNYQUIST, FIGSURF.
%
%   Author(s): J. Roubal (roubalj@control.felk.cvut.cz) [19-05-2006]



if nargin < 7
    disp('usage: FigPlot(FigHandle,FigXLabel,FigYLabel,FigAxesSize,FigTitle,FigTitleSize,PlotWidth)')
    return;
elseif nargin < 8 
    BaW = 0; LineStyle = {'-','-.','--',':','.'}; 
elseif nargin < 9 
    LineStyle = {'-','-.','--',':','.'};
end



CurrentAxesHandle = get(FigHandle,'CurrentAxes');

set(CurrentAxesHandle,'FontSize',FigAxesSize,'FontWeight','bold','XColor',[0.38 0.38 0.38],'YColor',[0.38 0.38 0.38]);

set(get(CurrentAxesHandle,'XLabel'),'String',FigXLabel,'FontSize',FigAxesSize,'FontWeight','bold','Color',[0 0 0]);
set(get(CurrentAxesHandle,'YLabel'),'String',FigYLabel,'FontSize',FigAxesSize,'FontWeight','bold','Color',[0 0 0]);
set(get(CurrentAxesHandle,'Title'),'String',FigTitle,'FontSize',FigTitleSize,'FontWeight','bold','Color',[0 0 0]);

PlotsHandle = get(CurrentAxesHandle,'Children');
for i = 1 : 1 : size(PlotsHandle,1)
    if BaW == 1, set(PlotsHandle(i),'LineWidth',PlotWidth,'Color','k','LineStyle',LineStyle{i});
    else set(PlotsHandle(i),'LineWidth',PlotWidth);
    end
end