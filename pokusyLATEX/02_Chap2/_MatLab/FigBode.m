function FigBode(FigHandle,FigXLabel,FigYLabelG,FigYLabelP,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW)

%   FIGBODE(FigHandle,FigXLabel,FigYLabelG,FigYLabelP,FigAxesSize,FigTitle,FigTitleSize,PlotWidth,BaW)
%       for Matlab R2006b (7.3.0.2)
%   FigBode sets parameters of the figure (FigHandle) with BODE-graph
%       - size of axes (FigAxesSize),
%       - FigXLabel, FigYLabelG, FigYLabelP and their sizes (FigAxesSize),
%       - FigTitle and its size (FigTitleSize),
%       - line width of the plots (PlotWidth),
%       - if you want black and white plots, set BaW=1
%         (if you want to change line style you must use bode(SYS,'--'), etc.)
%
%   Note: legend must be shown.
%
%   See also FIGCHAR, FIGNYQUIST, FIGPLOT, FIGSURF.
%
%   Author(s): J. Roubal (roubalj@control.felk.cvut.cz) [26-12-2007]



if nargin < 8
    disp('usage: FigBode(FigHandle,FigXLabel,FigYLabelG,FigYLabelP,FigAxesSize,FigTitle,FigTitleSize,PlotWidth)')
    return;
elseif nargin < 9
    BaW = 0;
end



XLabelHandle = xlabel(''); TitleHandle = title('');
set(XLabelHandle,'FontSize',FigAxesSize,'FontWeight','bold','String',FigXLabel);
set(TitleHandle,'FontSize',FigTitleSize,'FontWeight','bold','String',FigTitle);


% legend
HandleLeg=legend;
HandlePom=get(HandleLeg,'Children');
for i = 1 : 2 : size(HandlePom,1)
    set(get(HandlePom(i),'Children'),'LineWidth',PlotWidth);
end


FigChildren = get(FigHandle,'Children');
set(FigChildren(1),'FontSize',FigAxesSize,'FontWeight','bold'); %legenda
set(FigChildren(2),'FontSize',FigAxesSize,'FontWeight','bold');
set(get(FigChildren(2),'YLabel'),'String',FigYLabelP,'FontSize',FigAxesSize,'FontWeight','bold');
for i = 3 : 1 : size(FigChildren,1)
    switch get(FigChildren(i),'Type')
        case 'axes'
            set(FigChildren(i),'FontSize',FigAxesSize,'FontWeight','bold');
            set(get(FigChildren(i),'YLabel'),'String',FigYLabelG,'FontSize',FigAxesSize,'FontWeight','bold');
            j = i;
    end
end


MagnitudeHandleH = get(gcf,'Children');
MagnitudeHandle = get(MagnitudeHandleH(end),'Children');
PhaseHandle = get(gca,'Children');
if BaW == 1
    for i = 1 : 1 : size(PhaseHandle,1)
        set(get(MagnitudeHandle(i),'Children'),'LineWidth',PlotWidth,'Color',[0 0 0]);
        set(get(PhaseHandle(i),'Children'),'LineWidth',PlotWidth,'Color',[0 0 0]);
    end
else
    for i = 1 : 1 : size(PhaseHandle,1)
        set(get(MagnitudeHandle(i),'Children'),'LineWidth',PlotWidth);
        set(get(PhaseHandle(i),'Children'),'LineWidth',PlotWidth);
    end
end