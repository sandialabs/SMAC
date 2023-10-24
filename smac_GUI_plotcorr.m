function out=smac_GUI_plotcorr(action,varargin)
% SMAC_GUI_PLOTCORR  Calculate correlation coefficients
%
% OUT=SMAC_GUI_PLOTCORR
%
% SMAC_GUI_PLOTCORR shows a plot of the correlation coefficient plot.  The
% user can select a threshhold above which to select peaks, either by
% typing it in or by selecting off the graph.  A display on the right shows
% the currently selected peak frequencies and their associated correlation
% coefficient value.  These peaks will be used for the automatic
% optimization.  The plot above the correlation coefficient plot lets you
% toggle between viewing the NMIF or CMIF for these FRF.
%
% OUT will be 1 if the user selects the Initiate Auto SMAC button and the
% calculation is successful.  A -1 is returned if the user hits Backup, and
% is empty otherwise.

%==========================================================================
%
% ********************************************************************
% ***  Copyright (C) 2008 Sandia Corporation.                      ***
% ***                                                              ***
% ***  Under the terms of Contract DE-AC04-94AL85000 with Sandia   ***
% ***  Corporation, the U.S. Government retains certain rights in  ***
% ***  this software.                                              ***
% ********************************************************************
%
% The contents of this file are subject to the Mozilla Public License
% Version 1.1 (the "License"); you may not use this file except in
% compliance with the License. You may obtain a copy of the License at
% 
% http://www.mozilla.org/MPL/
% 
% Software distributed under the License is distributed on an "AS IS"
% basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
% License for the specific language governing rights and limitations under
% the License.
% 
% The Original Code is SMAC.
% 
% The Initial Developer of the Original Code is Randy Mayes.
% All Rights Reserved.
% 
% Contributor(s): Dan Hensley / ATA Engineering.

%==========================================================================
% HISTORY
%  22-May-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  15-Jun-2004 / ATA Engineering / Dan Hensley
%    o Minor wording change in tooltipstrings
%    o Fix bug with peak asterisks appearing on the wrong plot when
%      entering a value in the edit box
%    o Better error checking on edit value of correlation coefficient
%
%  27-Jul-2004 / ATA Engineering / Dan Hensley
%    o Fix minor GUI bug when only one root is selected
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Compute NMIF for each reference
%    o Add MMIF to available plot types
%    o Put correlation coefficient plots from individual references on
%      the plot
%    o Add button for individual peak picking
%
%  13-Aug-2004 / ATA Engineering / Dan Hensley
%    o Correlation coefficient plot adjustments
%
%  20-Aug-2004 / ATA Engineering / Dan Hensley
%    o Add Zoom and Zoom All buttons
%    o Add plot highlights when selecting roots from the list
%
%  03-Sep-2004 / ATA Engineering / Dan Hensley
%    o Add number of roots to the GUI list
%    o If more than 1 reference and we are fitting real normal modes, make
%      MMIF the default
%    o Add button to recompute # of repeated roots for selected freqs
%    o Add button to manually override # of repeated roots
%
%  10-Sep-2004 / ATA Engineering / Dan Hensley
%    o Fix some minor bugs
%
%  23-Sep-2004 / ATA Engineering / Dan Hensley
%    o Change Y axis scale for CMIF to log
%    o Make sure MIF plot X axis labels match the CC plot
%    o If zoomed in, make sure MIF plot is zoomed in when switching between
%      them
%
%  01-Oct-2004 / ATA Engineering / Dan Hensley
%    o Draw peak markers and selected roots on the MIF plot as well
%    o Rearrange root management buttons and add select All/None
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make GUI non-blocking (don't use uiwait/uiresume), and change
%      forward/backup calls to call the next/previous GUI explicitly
%    o Update for new .mif fields
%
%  09-Nov-2004 / ATA Engineering / Dan Hensley
%    o Set progress when user hits Backup button
%
%  17-Feb-2005 / ATA Engineering / Dan Hensley
%    o Shrink some of the GUI elements so it better fits on smaller screens
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Tweaks so this runs properly with Matlab 7.0.x and is still
%      compatible with Matlab 6.5.x
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Add resize figure callback to work around new Matlab 7 "feature"
%      that resizes axes differently than Matlab 6.5.
%
%  07-Mar-2005 / ATA Engineering / Dan Hensley
%    o Make sure colors on NMIF plot match CC plot for individual
%      references
%    o Add legend to NMIF plot
%
%  05-May-2005 / ATA Engineering / Dan Hensley
%    o Force axis Position after replotting to make sure the plots don't
%      automatically resize.
%
%  01-Jun-2005 / ATA Engineering / Dan Hensley
%    o Fix some MIF axis display problems when redrawing or resizing the GUI
%    o Fix listbox display problem when deleting all roots
%
%  07-Jun-2005 / ATA Engineering / Dan Hensley
%    o Use new GUI functions in uiWidgets
%
%  10-Aug-2005 / ATA Engineering / Dan Hensley
%    o If fitting real modes, calculate CMIF from imaginary part of FRF
%
%  03-Mar-2006 / ATA Engineering / Dan Hensley
%    o Allow the user to specify multiple roots for single reference data
%
%  07-Mar-2006 / ATA Engineering / Dan Hensley
%    o Resize the figure if it's too big to fit on the screen
%
%  10-Mar-2006 / ATA Engineering / Dan Hensley
%    o Redo figure resize and shrink the listbox a bit
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
    action='';
end

% Register valid callbacks

bname='uicbPlotCorr';
funcs.execute='Execute';
funcs.updateplot='UpdatePlot';
funcs.editcoef='EditCoef';
funcs.selectcoef='SelectCoef';
funcs.deletecoef='DeleteCoef';
funcs.pickpeak='PickPeak';
funcs.zoom='Zoom';
funcs.selectlist='SelectList';

funcs.roots_all='SelectRootsAll';
funcs.roots_none='SelectRootsNone';
funcs.countroots='CountRoots';
funcs.setroots='SetRoots';

funcs.plotnmif='PlotNMIF';
funcs.plotcmif='PlotCMIF';
funcs.plotmmif='PlotMMIF';

funcs.save='Save';
funcs.quit='Quit';
funcs.backup='Backup';

funcs.resizefig='ResizeFigure';

% Evaluate the callback if we have one

if ~isempty(action)
    if isfield(funcs,action),
        fh=str2func(strcat(bname,funcs.(action)));
        feval(fh,varargin{:});
        return;
    else
        error('Unknown action ''%s''',action);
    end
end

% Build the GUI

uicbPlotCorrCreateFigure(varargin);


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbPlotCorrCreateFigure(varargin)
% Create the Correlation Coefficients form

global ss;

scr=get(0,'ScreenSize');
scr=scr(3:4);

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'ResizeFcn','smac_GUI_plotcorr(''resizefig'',gcbf)', ...
    'Name','SMAC Correlation Plot', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_PlotCorrForm', ...
    'Visible','off');
hf=handle(hf);

%--------------------------------------------------------------------------
% Add Execute button

hbe=ui.CreateButton(hf,'Initiate Auto SMAC',[170 40],uip.Gap*[1 1]);
hbe.Callback='smac_GUI_plotcorr(''execute'',gcbf)';
hbe.TooltipString='Perform initial curve-fit of the roots selected in this form';
%uiFitText(hb);

%--------------------------------------------------------------------------
% Add Roots management options

listheight=800*scr(2)/1200;	% Adjust based on 1600x1200 height

hfr=ui.CreatePanel(hf,' Roots',[200 100],[uip.Gap uip.Gap]);

%----------------------------------
% Add Roots management buttons

hb(1)=ui.CreateButton(hfr,'None',[],[uip.Gap uip.Gap]);
hb(1).Callback='smac_GUI_plotcorr(''roots_none'',gcbf)';
hb(1).TooltipString='Deselect all roots';

hb(2)=ui.CreateButton(hfr,'All',[],'top+gaptop',hb(1));
hb(2).Callback='smac_GUI_plotcorr(''roots_all'',gcbf)';
hb(2).TooltipString='Select all roots';

hb(3)=ui.CreateButton(hfr,'Delete',[],'top+gaptop+gaptop',hb(2));
hb(3).Callback='smac_GUI_plotcorr(''deletecoef'',gcbf)';
hb(3).TooltipString='Delete selected roots';
ui.FitText(hb(3),[0 inf]);

hb(4)=ui.CreateButton(hfr,'Set # Roots',[],'top+gaptop',hb(3));
hb(4).Callback='smac_GUI_plotcorr(''setroots'',gcbf)';
hb(4).TooltipString=sprintf('Manually set the number of repeated roots to use for the selected frequencies');
ui.FitText(hb(4),[0 inf]);

hb(5)=ui.CreateButton(hfr,'Recount',[],'top+gaptop',hb(4));
hb(5).Callback='smac_GUI_plotcorr(''countroots'',gcbf)';
hb(5).TooltipString=sprintf('Recount the number of roots at each selected frequency\nusing the MIF function currently selected');
ui.FitText(hb(5),[0 inf]);

ui.ResizeControl(hfr);

% Add Roots listbox

listfmt='%3d: %8.2f  %5.3f  %4d';
titlfmt=sprintf('%3s: %8s  %5s  %4s','Num','Freq(Hz)','Corr','#Rts');

hl=ui.CreateList(hf,sprintf(listfmt,[0 0 0 0]),[170 listheight],'top+gaptop+gaptop',hbe);
hl.Callback='smac_GUI_plotcorr(''selectlist'',gcbf)';
hl.FontSize=10;
hl.Position(2)=hl.Position(2)+uip.Editbox.Size(2)+uip.Gap;
hl.Position(1)=hfr.Position(1) + hfr.Position(3) + uip.Gap;

ui.FitText(hl,[0 inf]);	% Only fit non-infinite dimension
hl.String='';
hl.Value=[];
setappdata(hf,'HandleList',hl);
setappdata(hf,'ListFormat',listfmt);

ht=ui.CreateText(hf,titlfmt,[],'top',hl);
pos=ht.Position; pos(1)=pos(1)+4; ht.Position=pos;
ht.FontName='Courier New';
ht.FontSize=10;
ui.FitText(ht);

% Move the Roots panel and Initiate button around

hfr.Position(2)=hl.Position(2)+hl.Position(4)-hfr.Position(4);
hbe.Position(1)=hl.Position(1);

% Make the buttons all the same size

width=0;
for k=1:length(hb), width=max(width,hb(k).Position(3)); end
for k=1:length(hb), hb(k).Position(3)=width; end

%--------------------------------------------------------------------------
% Create the MIF plot axis and radiobuttons

% axiswidth=900*scr(1)/1600;		% Adjust based on 1600x1200 resolution
% axisheight=250*scr(2)/1200;		% Adjust based on 1600x1200 resolution
asize=[800 250].*scr./[1600 1200];	% Resize based on 1600x1200 resolution

pos=hl.Position;
ha=ui.CreateAxis(hf,titlfmt,asize,'right+gapright+aligntop',hl);
ha.Position(1)=ha.Position(1)+50;
setappdata(hf,'AxisHandleMIF',ha);
setappdata(hf,'PositionAxisMIF',ha.Position);

hr=ui.CreateRadio(hf,'NMIF',[],'right+gapright+aligntop',ha);
hr.Callback='smac_GUI_plotcorr(''plotnmif'',gcbf)';
setappdata(hf,'HandleNMIFRadio',hr);

hr=ui.CreateRadio(hf,'CMIF',[],'bottom+gapbottom',hr);
hr.Callback='smac_GUI_plotcorr(''plotcmif'',gcbf)';
setappdata(hf,'HandleCMIFRadio',hr);

hr=ui.CreateRadio(hf,'MMIF',[],'bottom+gapbottom',hr);
hr.Callback='smac_GUI_plotcorr(''plotmmif'',gcbf)';
setappdata(hf,'HandleMMIFRadio',hr);

% Create the plot axis

pos=hl.Position;
ha=ui.CreateAxis(hf,titlfmt,[asize(1) pos(4)-50-asize(2)-30],'bottom+gapbottom',ha);
ha.Position(2)=ha.Position(2)-30;
setappdata(hf,'AxisHandle',ha);
setappdata(hf,'PositionAxis',ha.Position);
co=get(ha,'ColorOrder');
set(ha,'ColorOrder',co([2:end 1],:));

% Find the peaks

defcc=0.9;
if isempty(ss.corr.corrind),
    [~,ss.corr.corrind]=smac_find_peaks(ss.corr.corr(:,1),defcc);
end

% Update the plot in the axis and resize it

uicbPlotCorrUpdatePlot(hf,defcc);

% Update the MIF plot

if ss.realcomplex==1
    if length(ss.ref_coords)>1,
        smac_GUI_plotcorr('plotmmif',hf);
    else
        smac_GUI_plotcorr('plotnmif',hf);
    end
else
    smac_GUI_plotcorr('plotcmif',hf);
end

uicbPlotCorrUpdatePeaks(hf);

% Determine which peaks have multiple roots and update the list

if isempty(ss.corr.nroots),
    xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
    ss.corr.nroots=smac_getnumroots(ss,xf(ss.corr.corrind),uicbPlotCorrGetSelectedMIFType(hf));
end
uicbPlotCorrUpdateList(hf);

% Store the axis limits

setappdata(hf,'ZoomAllLimits',xlim);

% Add the minimum coefficient editbox

ht=ui.CreateText(hf,'Minimum Coefficient',[300 uip.Editbox.Size(2)],'right+alignbottom+gapright',hl);
ui.FitText(ht,[0 inf]);	% Only fit non-infinite dimension
ht.Position(2)=ht.Position(2)-ht.Position(4);

he=ui.CreateEdit(hf,num2str(defcc),[],'right+gapright',ht);
he.Callback='smac_GUI_plotcorr(''editcoef'',gcbf,gcbo)';
he.TooltipString='Enter a correlation coefficient above which to select peaks';
setappdata(hf,'HandleMaxCCEdit',he);
setappdata(hf,'MinCCOld',defcc);

hb=ui.CreateButton(hf,'Select',[uip.Button.Size(1) uip.Editbox.Size(2)],'right+gapright',he);
hb.Callback='smac_GUI_plotcorr(''selectcoef'',gcbf)';
hb.TooltipString='Select minimum coefficient by clicking on the graph';

hb=ui.CreateButton(hf,'Pick Peak',[uip.Button.Size(1) uip.Editbox.Size(2)],'right+gapright+gapright',hb);
hb.Callback='smac_GUI_plotcorr(''pickpeak'',gcbf)';
hb.TooltipString='Manually pick a peak by clicking on the graph';
ui.FitText(hb,[0 inf]);

hb=ui.CreateButton(hf,'Zoom All',[uip.Button.Size(1) uip.Editbox.Size(2)],'right+gapright',hb);
hb.Callback='smac_GUI_plotcorr(''zoom'',gcbf,''all'')';
hb.TooltipString='Zoom out to show the full frequency range';
ui.FitText(hb,[0 inf]);
hb.Position(1)=ha.Position(1)+ha.Position(3)-hb.Position(3);

hb=ui.CreateButton(hf,'Zoom',[uip.Button.Size(1) uip.Editbox.Size(2)],'left+gapleft',hb);
hb.Callback='smac_GUI_plotcorr(''zoom'',gcbf,''window'')';
hb.TooltipString='Zoom in on the selected frequency range';
ui.FitText(hb,[0 inf]);

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit,backup)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',ha);
hb.Callback='smac_GUI_plotcorr(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';
hb.Position(2)=uip.Gap;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''corre_plot'')';
hb.TooltipString='Get help on the Correlation Coefficient calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_plotcorr(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';

% Backup button

hb=ui.CreateButton(hf,'Backup',[],'left+gapleft',hb);
hb.Callback='smac_GUI_plotcorr(''backup'',gcbf)';
hb.TooltipString='Return to the correlation coefficient calcuation form';


%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

% Shrink the figure if it's bigger than the screen

pos=hf.Position(3:4);
if any(pos>.98)
    pos=pos.*(.98/max(pos));
end
hf.Position(3:4)=pos;

% Now make sure it's not off the screen

hf.Position(hf.Position<0.01)=0.01;

% Make the figure visible

hf.Visible='on';


%==========================================================================
function hf=uicbPlotCorrExecute(varargin)
% Initiate auto-SMAC

global ss;

hf=varargin{1};

% First make sure we have some peaks to process

if isempty(ss.corr.corrind),
    uiwait(errordlg('You must have some peaks selected!'));
    return;
end

smac_setprogress('selroots');

% Delete the figure and move to the next one
if ishandle(hf), delete(hf); end
smac_GUI_autofit;


%==========================================================================
function hf=uicbPlotCorrUpdatePlot(varargin)
% Update the plot

global ss;

hf=varargin{1};
ha=getappdata(hf,'AxisHandle');

% Get data to plot

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));

% Update the plot

set(hf,'CurrentAxes',ha);

% Plot coefficients from individual references

hold on
hl=plot(xf,ss.corr.corr_ref,'-','LineWidth',1);
legend(cellstr(char(ss.ref_coords)),'Location','NorthEast');
setappdata(hf,'HandleCCLines',hl);

% Plot the overall coefficient

plot(xf,ss.corr.corr(:,1),'b-','LineWidth',1)   %  Change from abs(Mcc) to Mcc aug 1 03 rlm

% Add labels

grid on
xlabel('Frequency (Hz)');
ylabel('Correlation Value');

xlim([xf(1) xf(end)]);
ylim([0 1]);

% Draw the max CC line if necessary

if nargin>1,
    ccval=varargin{2};
    hold on;
    hl=plot([xf(1) xf(end)],ccval*[1 1],'m');
    setappdata(hf,'HandleMaxCCLine',hl);
    hold off;
end


%==========================================================================
function hf=uicbPlotCorrSelectRootsAll(varargin)
% Select all of the roots in the list

hf=varargin{1};
hl=getappdata(hf,'HandleList');

hl.Value=1:size(hl.String,1);

% Update the plots

uicbPlotCorrSelectList(hf);


%==========================================================================
function hf=uicbPlotCorrSelectRootsNone(varargin)
% Deselect all of the roots in the list

hf=varargin{1};
hl=getappdata(hf,'HandleList');

hl.Value=[];

% Update the plots

uicbPlotCorrSelectList(hf);


%==========================================================================
function hf=uicbPlotCorrSelectList(varargin)
% Select peaks from the list and highlight on the graph

global ss;

hf=varargin{1};
ha=getappdata(hf,'AxisHandle');
hl=getappdata(hf,'HandleList');

% First delete the highlight peaks

hp=getappdata(hf,'HandleHighlights');
delete(hp(ishandle(hp)));
hp=[];

% Highlight the selected peaks

ind=get(hl,'Value');
if ~isempty(ind),

    % Redraw the peaks

    if ~isempty(ss.corr.corrind),
        xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
        set(hf,'CurrentAxes',ha);
        hold on;
        hp(1)=plot(xf(ss.corr.corrind(ind)),ss.corr.corr(ss.corr.corrind(ind),1),'go ','MarkerSize',7,'Color',[0 .85 0]);
        hold off;

        % Update the MIF plot as well

        ha=getappdata(hf,'AxisHandleMIF');
        set(hf,'CurrentAxes',ha);
        hpp=getappdata(hf,'HandlePeaks');
        xx=get(hpp(2),'XData');
        yy=get(hpp(2),'YData');
        hold on;
        hp(2)=plot(xx(ind),yy(ind),'go ','MarkerSize',7,'Color',[0 .85 0]);
        units=get(ha,'Units');
        set(ha,'Units','pixels','Position',getappdata(hf,'PositionAxisMIF'),'Units',units);
        hold off;

        % Save the handles

        setappdata(hf,'HandleHighlights',hp);
    end

end


%==========================================================================
function uicbPlotCorrEditCoef(varargin)
% Update the max CC line

global ss;

hf=varargin{1};
he=varargin{2};
ha=getappdata(hf,'AxisHandle');

ccval=str2num(get(he,'String'));

% Make sure a valid number was entered

if isempty(ccval),
    uiwait(errordlg('Value must be a numeric scalar'));
    defcc=getappdata(hf,'MinCCOld');
    set(he,'String',num2str(defcc));
    return;
end

% Make sure the number is within range

ccval=ccval(1);
yl=get(ha,'YLim');
if ccval<yl(1) || ccval>yl(2)
    uiwait(errordlg(sprintf('Value must be between %g and %g',yl)));
    return;
end

% Find the peaks again

[~,ss.corr.corrind]=smac_find_peaks(ss.corr.corr(:,1),ccval);

% Look at the peaks for multiple indices

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
ss.corr.nroots=smac_getnumroots(ss,xf(ss.corr.corrind),uicbPlotCorrGetSelectedMIFType(hf));

% Update the line, edit box, and list

uicbPlotCorrUpdateMaxLine(hf,ccval);
uicbPlotCorrUpdatePeaks(hf);

set(he,'String',num2str(ccval));
uicbPlotCorrUpdateList(hf);
uicbPlotCorrSelectList(hf);
setappdata(hf,'MinCCOld',ccval);

ss.done.selroots=false;


%==========================================================================
function uicbPlotCorrSelectCoef(varargin)
% Select the coefficient threshhold from the plot

global ss;

hf=varargin{1};
he=getappdata(hf,'HandleMaxCCEdit');
ha=getappdata(hf,'AxisHandle');

% Get the input from the user's click

[~,ccval]=ginput(1);

% Check against the valid limits

yl=get(ha,'YLim');

if ccval<yl(1) || ccval>yl(2)
    uiwait(errordlg(sprintf('Value must be between %g and %g',yl)));
    return;
end

% Find the peaks again

[~,ss.corr.corrind]=smac_find_peaks(ss.corr.corr(:,1),ccval);

% Look at the peaks for multiple indices

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
ss.corr.nroots=smac_getnumroots(ss,xf(ss.corr.corrind),uicbPlotCorrGetSelectedMIFType(hf));

% Update the line, edit box, and list

uicbPlotCorrUpdateMaxLine(hf,ccval);
uicbPlotCorrUpdatePeaks(hf);

set(he,'String',num2str(ccval));
uicbPlotCorrUpdateList(hf);
uicbPlotCorrSelectList(hf);
setappdata(hf,'MinCCOld',ccval);

ss.done.selroots=false;


%==========================================================================
function uicbPlotCorrUpdateList(varargin)
% Update the max CC line

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');
listfmt=getappdata(hf,'ListFormat');

% Get the data for the list

freqs=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
ind=ss.corr.corrind;
corr=ss.corr.corr(:,1);

% Update the listbox

str=cell(length(ind),1);
for k=1:length(ind),
    str{k}=sprintf(listfmt,k,freqs(ind(k)),corr(ind(k)),ss.corr.nroots(k));
end

if length(ind)==1, val=1; else val=[]; end
%drawnow;
set(hl,'Min',1,'Max',3,'String',str,'Value',val);
%drawnow;


%==========================================================================
function uicbPlotCorrPickPeak(varargin)
% Pick a peak from the graph

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

% Get the point picked from the graph

ha=getappdata(hf,'AxisHandle');
set(hf,'CurrentAxes',ha);

[x,y]=ginput(1);

% Figure out the closest point to the one picked

xx=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
yy=ss.corr.corr(:,1);

xr=diff(xlim);
yr=diff(ylim);

dist=sqrt(((xx-x)/xr).^2 + ((yy-y)/yr).^2);	% Distance, normalized

[~,ind]=min(dist);

% hold on;
% plot(xx(ind)*[1 1],yy(ind)*[1 1],'ro')

% Insert the point into the index listing

[ss.corr.corrind,ind2]=sort([ss.corr.corrind; ind]);
ss.corr.nroots(end+1)=0;
ss.corr.nroots=ss.corr.nroots(ind2);

[~,ind2]=max(ind2);


% Determine the number of roots at the selected index

ss.corr.nroots(ind2)=smac_getnumroots(ss,xx(ind),uicbPlotCorrGetSelectedMIFType(hf));

% Auto select the root we just added

set(hl,'Value',ind2);

% Update the listbox

uicbPlotCorrUpdateList(hf);
uicbPlotCorrUpdatePeaks(hf);
uicbPlotCorrSelectList(hf);


%==========================================================================
function uicbPlotCorrDeleteCoef(varargin)
% Delete the selected coefficients

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

% Remove the selected roots

val=get(hl,'Value');
if ~isempty(val),
    ss.corr.corrind(val)=[];
    ss.corr.nroots(val)=[];

    % Delete the highlighted peaks

    hp=getappdata(hf,'HandleHighlights');
    delete(hp);
    hp=[];
    setappdata(hf,'HandleHighlights',hp);

    % Update the listbox

    uicbPlotCorrUpdateList(hf);
    uicbPlotCorrUpdatePeaks(hf);
end


%==========================================================================
function hf=uicbPlotCorrCountRoots(varargin)
% Count the number of repeated roots for the selected ones

global ss;

hf=varargin{1};
ha=getappdata(hf,'AxisHandle');
hl=getappdata(hf,'HandleList');

% See what peaks are highlighted

ind=get(hl,'Value');
if ~isempty(ind),

    % Count the repeated roots

    if ~isempty(ss.corr.corrind),
        xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
        ss.corr.nroots(ind)=smac_getnumroots(ss,xf(ss.corr.corrind(ind)),uicbPlotCorrGetSelectedMIFType(hf));

        % Update the listbox

        uicbPlotCorrUpdateList(hf);
        set(hl,'Value',ind)
    end
end


%==========================================================================
function hf=uicbPlotCorrSetRoots(varargin)
% Manually set the number of repeated roots for the selected ones

global ss;

hf=varargin{1};
ha=getappdata(hf,'AxisHandle');
hl=getappdata(hf,'HandleList');

% See what peaks are highlighted

ind=get(hl,'Value');
if isempty(ind),
    uiwait(warndlg('Please select at least one root to override'));
    return;

else

    % Ask for the number to use

    num=inputdlg('Enter number of repeated roots');
    if isempty(num), return; end
    nov=str2num(num{1});

    % Check (integer 1<#<# of refs)

%     if numel(nov)~=1 | ~ismember(nov,1:length(ss.ref_coords)),
%         uiwait(errordlg(sprintf('Number must be a scalar between 1 and %d',length(ss.ref_coords))));
%         return;
%     end

    % Set the number of roots

    if ~isempty(nov),
        ss.corr.nroots(ind)=nov;

        % Update the listbox

        uicbPlotCorrUpdateList(hf);
        set(hl,'Value',ind)
    end
end


%==========================================================================
function uicbPlotCorrUpdateMaxLine(varargin)
% Update the max CC line

global ss;

hf=varargin{1};
ccval=varargin{2};

hl=getappdata(hf,'HandleMaxCCLine');
set(hl,'YData',ccval*[1 1]);


%==========================================================================
function uicbPlotCorrUpdatePeaks(varargin)
% Update the peaks

global ss;

hf=varargin{1};

% Draw or update the peaks if they have been found

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));

% First delete the peaks

hp=getappdata(hf,'HandlePeaks');
delete(hp(ishandle(hp)));
hp=[];

% Redraw the peaks

if ~isempty(ss.corr.corrind),
    set(0,'CurrentFigure',hf);
    ha=getappdata(hf,'AxisHandle');
    set(hf,'CurrentAxes',ha);
    hold on;
    hp(1)=plot(xf(ss.corr.corrind),ss.corr.corr(ss.corr.corrind,1),'r* ');
    hold off;
    units=get(ha,'Units');
    set(ha,'Units','pixels','Position',getappdata(hf,'PositionAxis'),'Units',units);

    % Get Y Axis data based on MIF type

    miftype=uicbPlotCorrGetSelectedMIFType(hf);
    switch miftype,
        case 'nmif',
            ydata=min(ss.mif.exp.nmif.ordinate(ss.corr.corrind+ss.freqrangecc(1,2)-1,:),[],2);
        case 'cmif',
            ydata=ss.mif.exp.cmif(1).ordinate(ss.corr.corrind+ss.freqrangecc(1,2)-1);
        case 'mmif',
            ydata=ss.mif.exp.mmif(1).ordinate(ss.corr.corrind+ss.freqrangecc(1,2)-1);
        otherwise,
            ydata=zeros(size(ss.corr.corrind));
    end

    % Draw the peaks on the MIF plot

    ha=getappdata(hf,'AxisHandleMIF');
    set(hf,'CurrentAxes',ha);

    hold on;
    hp(2)=plot(xf(ss.corr.corrind),ydata,'r* ');
    hold off;
    units=get(ha,'Units');
    set(ha,'Units','pixels','Position',getappdata(hf,'PositionAxisMIF'),'Units',units);

end
setappdata(hf,'HandlePeaks',hp);


%==========================================================================
function uicbPlotCorrPlotNMIF(varargin)
% Plot the NMIF

global ss;

hf=varargin{1};
hr(1)=getappdata(hf,'HandleNMIFRadio');
hr(2)=getappdata(hf,'HandleCMIFRadio');
hr(3)=getappdata(hf,'HandleMMIFRadio');
ha=getappdata(hf,'AxisHandleMIF');

% Toggle the radio buttons

set(hr,'Value',0);
set(hr(1),'Value',1);

% Compute the NMIF if necessary

if isempty(ss.mif.exp.nmif),
    hh=watchon;
    drawnow;
    ss.mif.exp.nmif=nmif(ss.fe);
    ss.mif.exp.nmif=ss.mif.exp.nmif(1:length(ss.ref_coords));
    watchoff(hh);
end

% Now plot it

set(hf,'CurrentAxes',ha);
hl=plot(ss.mif.exp.nmif.abscissa,ss.mif.exp.nmif.ordinate);
ylabel('NMIF');
grid on;
xlim([ss.freqrangecc(:,1)]);
set(ha,'Position',getappdata(hf,'PositionAxisMIF'));

% Use the X axis labels from the other graph

ha2=getappdata(hf,'AxisHandle');
set(ha,'XTick',get(ha2,'XTick'));
set(ha,'XLim',get(ha2,'XLim'));

% Change colors to match the CC plot

hl2=getappdata(hf,'HandleCCLines');
for k=1:length(hl),
    set(hl(k),'Color',get(hl2(k),'Color'));
end

ha2=getappdata(hf,'AxisHandle');

% Add a legend showing all of the references

legend(cellstr(char(ss.ref_coords)),'Location','NorthEast');

% Update the peak markers

uicbPlotCorrUpdatePeaks(hf);
uicbPlotCorrSelectList(hf);


%==========================================================================
function uicbPlotCorrPlotCMIF(varargin)
% Plot the CMIF

global ss;

hf=varargin{1};
hr(1)=getappdata(hf,'HandleNMIFRadio');
hr(2)=getappdata(hf,'HandleCMIFRadio');
hr(3)=getappdata(hf,'HandleMMIFRadio');
ha=getappdata(hf,'AxisHandleMIF');

% Toggle the radio buttons

set(hr,'Value',0);
set(hr(2),'Value',1);

% Compute the CMIF if necessary

if isempty(ss.mif.exp.cmif),
    hh=watchon;
    drawnow;
	if ss.realcomplex==1		% If fitting real modes, compute CMIF from imag(frf)
		ss.mif.exp.cmif=cmif(imag(ss.fe));
	else
		ss.mif.exp.cmif=cmif(ss.fe);
	end
    watchoff(hh);
end

% Now plot it

set(hf,'CurrentAxes',ha);
cmifo=ss.mif.exp.cmif.ordinate;
semilogy(ss.mif.exp.cmif.abscissa,cmifo);
ylabel('CMIF');
grid on;
xlim([ss.freqrangecc(:,1)]);
set(ha,'Position',getappdata(hf,'PositionAxisMIF'));

% Only display 4 decades at most

yl=[min(cmifo(:)) max(cmifo(:))];
yl(1)=max(yl(1),yl(2)/10^4);
ylim(yl)

% Use the X axis labels from the other graph

ha2=getappdata(hf,'AxisHandle');
set(ha,'XTick',get(ha2,'XTick'));
set(ha,'XLim',get(ha2,'XLim'));

% Update the peak markers

uicbPlotCorrUpdatePeaks(hf);
uicbPlotCorrSelectList(hf);


%==========================================================================
function uicbPlotCorrPlotMMIF(varargin)
% Plot the MMIF

global ss;

hf=varargin{1};
hr(1)=getappdata(hf,'HandleNMIFRadio');
hr(2)=getappdata(hf,'HandleCMIFRadio');
hr(3)=getappdata(hf,'HandleMMIFRadio');
ha=getappdata(hf,'AxisHandleMIF');

% Toggle the radio buttons

set(hr,'Value',0);
set(hr(3),'Value',1);

% Compute the MMIF if necessary

if isempty(ss.mif.exp.mmif),
    hh=watchon;
    drawnow;
    ss.mif.exp.mmif=mmif(ss.fe,'silent','nofp');
    watchoff(hh);
end

% Now plot it

set(hf,'CurrentAxes',ha);
plot(ss.mif.exp.mmif.abscissa,ss.mif.exp.mmif.ordinate);
ylabel('MMIF');
grid on;
xlim([ss.freqrangecc(:,1)]);
set(ha,'Position',getappdata(hf,'PositionAxisMIF'));


% Use the X axis labels from the other graph

ha2=getappdata(hf,'AxisHandle');
set(ha,'XTick',get(ha2,'XTick'));
set(ha,'XLim',get(ha2,'XLim'));

% Update the peak markers

uicbPlotCorrUpdatePeaks(hf);
uicbPlotCorrSelectList(hf);

%==========================================================================
function hf=uicbPlotCorrZoom(varargin)
% Zoom in or out on the plot

global ss;

hf=varargin{1};
zoomtype=varargin{2};
ha(1)=getappdata(hf,'AxisHandle');
ha(2)=getappdata(hf,'AxisHandleMIF');

switch zoomtype,
    case 'all',
        xl=getappdata(hf,'ZoomAllLimits');

    case 'window',
        [xl,~]=ginput(2);
        xl=sort(xl);

    otherwise
        uiwait(errordlg(sprintf('Internal zoom error:  Unknown zoom option ''%s''',zoomtype)));
        return;

end

% Set the axis limits

set(hf,'CurrentAxes',ha(1));
xlim(xl);
set(hf,'CurrentAxes',ha(2));
xlim(xl);
set(ha(2),'XTick',get(ha(1),'XTick'));


%==========================================================================
function miftype=uicbPlotCorrGetSelectedMIFType(varargin)
% Get the selected MIF type

hf=varargin{1};

% Determine the selected MIF type

if get(getappdata(hf,'HandleNMIFRadio'),'Value'),
    miftype='nmif';
end
if get(getappdata(hf,'HandleCMIFRadio'),'Value'),
    miftype='cmif';
end
if get(getappdata(hf,'HandleMMIFRadio'),'Value'),
    miftype='mmif';
end


%==========================================================================
function hf=uicbPlotCorrResizeFigure(varargin)
% Gets called when we resize the figure

hf=gcbf;

ham=getappdata(hf,'AxisHandleMIF');
ha=getappdata(hf,'AxisHandle');

units_ha=ha.Units;
units_ham=ham.Units;

ha.Units='pixels';
ham.Units='pixels';

ha.Position([1 3])=ham.Position([1 3]);

setappdata(hf,'PositionAxisMIF',ham.Position);
setappdata(hf,'PositionAxis',ha.Position);

ha.Units=units_ha;
ham.Units=units_ham;


%==========================================================================
function hf=uicbPlotCorrSave(varargin)
% Save current structure to .mat file

global ss;
global SMAC_SAVE_V6;

[fname,pname]=uiputfile('*.mat','Select .mat file to save to');
if isnumeric(fname), return; end

fname=fullfile(pname,fname);
if isempty(findstr(fname,'.mat')), fname=[fname '.mat']; end

if SMAC_SAVE_V6,
    save('-v6',fname,'ss');
else
    save(fname,'ss');
end


%==========================================================================
function hf=uicbPlotCorrBackup(varargin)
% Backup

hf=varargin{1};

if ishandle(hf), delete(hf); end
smac_setprogress('pinv');
smac_GUI_corrcoef;


%==========================================================================
function hf=uicbPlotCorrQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strmatch(rslt,'Yes'),
    delete(hf);
end
