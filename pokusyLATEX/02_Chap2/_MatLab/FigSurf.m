function FigSurf(FigHandle,FigXLabel,FigYLabel,FigZLabel,FigAxesSize,FigTitle,FigTitleSize)

%   FigSurf(FigHandle,FigXLabel,FigYLabel,FigZLabel,FigAxesSize,FigTitle,FigTitleSize)
%
%   FigSurf sets parameters of the figure (FigHandle) with SURF, MESH, etc.
%       - size of axes (FigAxesSize),
%       - FigXLabel, FigYLabel, FigZLabel and their sizes (FigAxesSize),
%       - FigTitle and its size (FigTitleSize),
%       - line width of the plots (PlotWidth).
%
%   See also FIGBODE, FIGCHAR, FIGNYQUIST, FIGPLOT.
%
%   Author(s): J. Roubal (roubalj@control.felk.cvut.cz) [25-07-2005]



if nargin < 7
    disp('usage: FigSurf(FigHandle,FigXLabel,FigYLabel,FigZLabel,FigAxesSize,FigTitle,FigTitleSize)')
    return;
end



CurrentAxesHandle = get(FigHandle,'CurrentAxes');

set(CurrentAxesHandle,'FontSize',FigAxesSize);

set(get(CurrentAxesHandle,'XLabel'),'String',FigXLabel,'FontSize',FigAxesSize,'FontWeight','bold');
set(get(CurrentAxesHandle,'YLabel'),'String',FigYLabel,'FontSize',FigAxesSize,'FontWeight','bold');
set(get(CurrentAxesHandle,'ZLabel'),'String',FigZLabel,'FontSize',FigAxesSize,'FontWeight','bold');
set(get(CurrentAxesHandle,'Title'),'String',FigTitle,'FontSize',FigTitleSize,'FontWeight','bold');

