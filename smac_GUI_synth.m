function out=smac_GUI_synth(action,varargin)
% SMAC_GUI_SYNTH  Final root management and shape synthesis
%
% OUT=SMAC_GUI_SYNTH
%
% SMAC_GUI_SYNTH allows the user to perform final root management and
% generate and animate mode shapes.  The user can manually add roots,
% delete selected roots, and condense roots with 0 frequency.  In addition,
% the user can synthesize NMIF, CMIF, and/or FRF using selected roots to
% compare with the experimental data.
%
% OUT provides the status of the GUI when it is closed.  A -1 is returned
% if the user hits Backup, and is empty otherwise.

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
%  02-Jun-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  08-Jun-2004 / ATA Engineering / Dan Hensley
%    o Connect Add Roots form to this GUI
%    o Residual Terms only valid for complex mode fit
%    o Prompt user when deleting roots
%
%  09-Jul-2004 / ATA Engineering / Dan Hensley
%    o Add MAC button to the form
%    o Autofit text in Create Shapes button
%
%  17-Sep-2004 / ATA Engineering / Dan Hensley
%    o Highlight all roots by default when entering the form
%    o Add MMIF synthesize button
%
%  24-Sep-2004 / ATA Engineering / Dan Hensley
%    o Add initial damping column to the list
%    o Put up error message when trying to synthesize a zeroed root
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Rearrange Synthesize and Close Plots buttons
%    o Fix bugs with plotting MIFs
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make GUI non-blocking (don't use uiwait/uiresume), and change
%      forward/backup calls to call the next/previous GUI explicitly
%    o New function to delete analytical data when needed
%    o Only calculate synthesis when necessary (store it)
%    o New Save Shapes button, moved functionality from Create Shapes
%    o Added Status titleframe for FEM, shape, synthesis
%
%  09-Nov-2004 / ATA Engineering / Dan Hensley
%    o Prompt user for shape selection technique
%    o Set progress when user hits Backup button
%    o Trap for condensed roots when creating shapes
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Add .fit.corr_ref handling for root management
%    o Support 2nd output argument from smac_GUI_addroot
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Show initial damping as % rather than fraction
%    o Make FEM status buttons 'inactive' instead of 'off'
%    o Tweaks so this runs properly with Matlab 7.0.x and is still
%      compatible with Matlab 6.5.x
%    o Add diagnostics display button when creating shapes
%
%  10-Mar-2005 / ATA Engineering / Dan Hensley
%    o Fix bug where residuals are not always calculated
%
%  07-Jun-2005 / ATA Engineering / Dan Hensley
%    o Use new GUI functions in uiWidgets
%    o Add call to smac_GUI_residuals
%
%  08-Jun-2005 / ATA Engineering / Dan Hensley
%    o Change layout and style of residuals GUI
%
%  15-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remember selected roots when entering the form again
%    o Add button for Resynthesize
%    o When saving shapes, also save DOF order shape for residuals
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Write out low and high mode residual frequency ranges
%
%  06-Mar-2006 / ATA Engineering / Dan Hensley
%    o Add ref and DP coef columns to listbox
%    o Resize the figure if it's too big to fit on the screen
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
    action='';
end

% Register valid callbacks

bname='uicbSynth';
funcs.createshapes='CreateShapes';
funcs.plotshapes='PlotShapes';
funcs.saveshapes='SaveShapes';
funcs.plotmac='PlotMAC';
funcs.resynthesize='Resynthesize';

funcs.roots_all='SelectAll';
funcs.roots_none='SelectNone';
funcs.roots_select='SelectFromList';
funcs.roots_add='AddRoots';
funcs.roots_condense='CondenseRoots';
funcs.roots_delete='DeleteRoots';

funcs.toggleresiduals='ToggleResiduals';
funcs.residuals='Residuals';
funcs.synthesize='Synthesize';
funcs.closesynthplots='CloseSynthPlots';
funcs.clearanalysisdata='ClearAnalysisData';

funcs.CMIFtoggle='CMIFtoggle';
funcs.UpdateTracking = 'UpdateTracking';

funcs.save='Save';
funcs.quit='Quit';
funcs.backup='Backup';

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

uicbSynthCreateFigure(varargin);


%==========================================================================
% CALLBACKS
%==========================================================================
%% uicbSynthCreateFigure
function hf=uicbSynthCreateFigure(varargin)
% Create the Synthesis form

global ss;

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Synthesis and Mode Shape Generator', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_SynthForm', ...
    'Visible','off');
hf=handle(hf);

%--------------------------------------------------------------------------
% Add Create Shapes button

hb=ui.CreateButton(hf,'Create Shapes',[170 40],uip.Gap*[1 1]);
hb.Callback='smac_GUI_synth(''createshapes'',gcbf)';
hb.TooltipString='Create mode shapes for the selected roots and save to an ADF';
ui.FitText(hb,[0 inf]);

hb3=ui.CreateCheckbox(hf,'Display diagnostics when creating shapes',[170 40],'top',hb);
hb3.TooltipString='Show diagnostic output when generating shapes';
hb3.Value=1;
ui.FitText(hb3,[0 0]);
setappdata(hf,'HandleShapeDiagnostics',hb3);

%--------------------------------------------------------------------------
% Add Plot Shapes button

hb2=ui.CreateButton(hf,'Plot Shapes',[170 40],'right+gapright',hb);
hb2.Callback='smac_GUI_synth(''plotshapes'',gcbf)';
hb2.TooltipString='Plot and animate mode shapes';
ui.FitText(hb2,[0 inf]);

% Add Save Shapes button

hb2=ui.CreateButton(hf,'Save Shapes',[170 40],'right+gapright',hb2);
hb2.Callback='smac_GUI_synth(''saveshapes'',gcbf)';
hb2.TooltipString='Save mode shapes to a shape ADF (ASH)';
ui.FitText(hb2,[0 inf]);

% Add Plot MAC button

hb2=ui.CreateButton(hf,'Plot MAC',[170 40],'right+gapright+gapright',hb2);
hb2.Callback='smac_GUI_synth(''plotmac'',gcbf)';
hb2.TooltipString='Generate and plot a Modal Assurance Criteria for the generated shapes';
ui.FitText(hb2,[0 inf]);

% Add Resynthesize button

hb2=ui.CreateButton(hf,'Resynthesize',[170 40],'right+gapright+gapright',hb2);
hb2.Callback='smac_GUI_synth(''resynthesize'',gcbf)';
hb2.TooltipString='Resynthesize FRF from mode shapes and residual information';
ui.FitText(hb2,[0 inf]);

% Add List Info button

hb2=ui.CreateButton(hf,'List Info',[170 40],'right+gapright+gapright',hb2);
hb2.Callback='smac_list_ref_info';
hb2.TooltipString='List information on what references produced each correlation coefficient and shape';
ui.FitText(hb2,[0 inf]);

%--------------------------------------------------------------------------
% Add Status panel and options

hfrs=ui.CreatePanel(hf,' Status',[400 100],'top',hb3);

hfrs.Position(2)=hfrs.Position(2)+20;

% Add status boxes

cdata=zeros(17,17,3);
cdata(:,:,1)=1;

hs(1)=ui.CreateButton(hfrs,' ',[20 20],[uip.Gap uip.Gap]);
hs(1).Enable='inactive';
hs(1).CData=cdata;
ht(1)=ui.CreateText(hfrs,'FEM Loaded',[inf 20],'right+gapright',hs(end));
hs(1).TooltipString='Specifies whether FEM has been loaded';
ht(1).TooltipString=hs(1).TooltipString;
setappdata(hf,'StatusFEM',hs(1));

hs(2)=ui.CreateButton(hfrs,' ',[20 20],'top+gaptop',hs(end));
hs(2).Enable='inactive';
hs(2).CData=cdata;
ht(2)=ui.CreateText(hfrs,'Shapes Created',[inf 20],'right+gapright',hs(end));
hs(2).TooltipString='Specifies whether shapes have been generated';
ht(2).TooltipString=hs(end).TooltipString;
setappdata(hf,'StatusShapes',hs(2));

hs(3)=ui.CreateButton(hfrs,' ',[20 20],'top+gaptop',hs(end));
hs(3).Enable='inactive';
hs(3).CData=cdata;
ht(3)=ui.CreateText(hfrs,'FRF Synthesized',[inf 20],'right+gapright',hs(end));
hs(3).TooltipString='Specifies whether shapes have been generated';
ht(3).TooltipString=hs(end).TooltipString;
setappdata(hf,'StatusSynthesis',hs(3));

ui.ResizeControl(hfrs);

uicbSynthUpdateStatus(hf);


%--------------------------------------------------------------------------
% Add Synthesis panel and options

hfr=ui.CreatePanel(hf,' Synthesis',[400 100],'top+gaptop+gaptop',hfrs);

%----------------------------------
% Add Synthesize and Close Plots buttons

hb=ui.CreateButton(hfr,'Close Plots',[],[uip.Gap uip.Gap]);
hb.Callback='smac_GUI_synth(''closesynthplots'',gcbf)';
hb.TooltipString='Close all of the synthesis plots';

hb2=ui.CreateButton(hfr,'Synthesize',[],'top+gaptop',hb);
hb2.Callback='smac_GUI_synth(''synthesize'',gcbf)';
hb2.TooltipString='Compare experimental to analytical data for the selected roots';

%----------------------------------
% Add synthesis options toggles

hc1=ui.CreateCheckbox(hfr,'Use Residuals',[],'top+gaptop+gaptop',hb2);
hc1.Callback='smac_GUI_synth(''toggleresiduals'',gcbf)';
hc1.TooltipString='Toggle for whether to include residual terms';
hc1.Value=ss.residuals.use;
ui.FitText(hc1,[0 inf]);
setappdata(hf,'HandleResidualToggle',hc1);

hc2=ui.CreateButton(hfr,'Specify',[],'right+gapright',hc1);
hc2.Callback='smac_GUI_synth(''residuals'',gcbf)';
hc2.TooltipString='Specify what residual terms to apply';
ui.FitText(hc2,[0 inf]);
setappdata(hf,'HandleResiduals',hc2);

hc=ui.CreateCheckbox(hfr,'FRF',[],'top+gaptop+gaptop+gaptop',hc1);
hc.TooltipString='Synthesize frequency response functions';
setappdata(hf,'ToggleFRF',hc);
ui.FitText(hc);

hc=ui.CreateCheckbox(hfr,'CMIF',[],'top+gaptop',hc);
hc.TooltipString='Synthesize complex mode indicator functions';
hc.Callback='smac_GUI_synth(''CMIFtoggle'',gcbf)';
if ss.realcomplex==2, hc.Value=1; end
setappdata(hf,'ToggleCMIF',hc);
ui.FitText(hc);

hc3=ui.CreateCheckbox(hfr,'Use ref tracking',[],'right+gapright+gapright',hc);
hc3.TooltipString='Synthesize CMIF using reference tracking';
hc3.Value=0;
hc3.Callback='smac_GUI_synth(''UpdateTracking'',gcbf)';
setappdata(hf,'ToggleTracking',hc3);
ui.FitText(hc3);

hc=ui.CreateCheckbox(hfr,'MMIF',[],'top+gaptop',hc);
hc.TooltipString='Synthesize multivariate mode indicator functions';
if ss.realcomplex==1 && length(ss.ref_coords)>1, hc.Value=1; end
setappdata(hf,'ToggleMMIF',hc);
ui.FitText(hc);

hc=ui.CreateCheckbox(hfr,'NMIF',[],'top+gaptop',hc);
hc.TooltipString='Synthesize normal mode indicator functions';
if ss.realcomplex==1 && length(ss.ref_coords)==1, hc.Value=1; end
setappdata(hf,'ToggleNMIF',hc);
ui.FitText(hc);

hfr.Position(3)=hfrs.Position(3);
ui.ResizeControl(hfr);
hb.Position(3)=hfr.Position(3)-2*uip.Gap;
hb2.Position(3)=hfr.Position(3)-2*uip.Gap;

%--------------------------------------------------------------------------
% Add Roots panel and options

hfr2=ui.CreatePanel(hf,' Roots',[400 100],'top+gaptop+gaptop',hfr);

%----------------------------------
% Add Roots management buttons

hba=ui.CreateButton(hfr2,'All',[],[uip.Gap uip.Gap]);
hba.Callback='smac_GUI_synth(''roots_all'',gcbf)';
hba.TooltipString='Select all roots';

hb(1)=ui.CreateButton(hfr2,'None',[],'right+gapright',hba);
hb(1).Callback='smac_GUI_synth(''roots_none'',gcbf)';
hb(1).TooltipString='Deselect all roots';

hb(2)=ui.CreateButton(hfr2,'Condense',[],'top+gaptop+gaptop',hba);
hb(2).Callback='smac_GUI_synth(''roots_condense'',gcbf)';
hb(2).TooltipString='Remove zeroed roots';

hb(3)=ui.CreateButton(hfr2,'Merge',[],'top+gaptop',hb(2));
hb(3).Callback='smac_merge_roots';
hb(3).TooltipString='Merge root in from an ASH file';

hb(4)=ui.CreateButton(hfr2,'Delete',[],'top+gaptop',hb(3));
hb(4).Callback='smac_GUI_synth(''roots_delete'',gcbf)';
hb(4).TooltipString='Delete selected roots';

hb(5)=ui.CreateButton(hfr2,'Add',[],'top+gaptop',hb(4));
hb(5).Callback='smac_GUI_synth(''roots_add'',gcbf)';
hb(5).TooltipString='Manually add roots';

ui.ResizeControl(hfr2);

wide=max([hfrs.Position(3) hfr.Position(3) hfr2.Position(3)]);
hfr.Position(3)=wide;
hfr2.Position(3)=wide;
hfrs.Position(3)=wide;

% Resize the top 4 buttons to be the width of the frame

for k=2:5
    hb(k).Position(3)=wide-2*uip.Gap;
end

% Now move the None button over just in case the lower frame was wider

hb(1).Position(1)=hfr2.Position(3)-hb(1).Position(3)-uip.Gap;

%--------------------------------------------------------------------------
% Add Roots listbox and title

fieldwidth=[4 10 8 6 10 8 6];

listfmt=sprintf('%% %dd  %%%d.3f  %%%d.3f  %%%d.3f    %%%d.3f  %%%d.3f  %%%d.3f',fieldwidth);
listfmt2='  %10s  %10s%c';
ttlfmt =sprintf('%%%ds  %%%ds  %%%ds  %%%ds    %%%ds  %%%ds  %%%ds',fieldwidth);
ttlfmt2=sprintf('%%%ds  %%-%ds  %%%ds  %%%ds    %%-%ds  %%%ds  %%%ds',fieldwidth);
titlfmt=sprintf([ttlfmt listfmt2], ...
                'Root','Freq (Hz)','Damp (%)','Corr', ...
                'Freq (Hz)','Damp (%)','Corr', ...
                'Shp Ref','DP coef',' ' );
titltop=sprintf(ttlfmt2,'','FINAL','','','INITIAL');

hl=ui.CreateList(hf,sprintf([listfmt listfmt2],[1 0 0 0 0 0 0],'','',' '),[],'right+gapright',hfrs);
hl.Callback='smac_GUI_synth(''roots_select'',gcbf)';
hl.FontSize=10;
hl.Position(4)=hfr2.Position(2)+hfr2.Position(4)-hfrs.Position(2);
ui.FitText(hl,[0 inf]);	% Only fit non-infinite dimension

setappdata(hf,'HandleList',hl);
setappdata(hf,'ListFormat',listfmt);
setappdata(hf,'ListFormat2',listfmt2);

ht=ui.CreateText(hf,titlfmt,[],'top',hl);
pos=ht.Position; pos(1)=pos(1)+4; ht.Position=pos;
ht.FontName='Courier New';
ht.FontSize=10;
ui.FitText(ht);

ht=ui.CreateText(hf,titltop,[],'top',ht);
%pos=ht.Position;, pos(1)=pos(1)+4;, ht.Position=pos;
ht.FontName='Courier New';
ht.FontSize=10;
ui.FitText(ht);

% Draw some vertical lines

height=ht.Position(2)+ht.Position(4)-hl.Position(2);
textwidth=ht.Position(3)/length(ht.String);

hv=ui.CreateVertLine(hf,'',height,'',hl);
hv.Position(1)=hv.Position(1)+(fieldwidth(1)+2)*textwidth+2;

hv=ui.CreateVertLine(hf,'',height,'',hl);
hv.Position(1)=hv.Position(1)+(sum(fieldwidth(1:4))+2*5)*textwidth+2;

hv=ui.CreateVertLine(hf,'',height,'',hl);
hv.Position(1)=hv.Position(1)+(sum(fieldwidth(1:7))+2*9)*textwidth+2;

% Fill in the initial text

uicbSynthUpdateList(hf);
hl.Value=ss.shape.rootsel;

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit,backup)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',hl);
hb.Callback='smac_GUI_synth(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';
hb.Position(2)=uip.Gap;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''sh_hlp'')';
hb.TooltipString='Get help on the Synthesis calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_synth(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';

% Backup button

hb=ui.CreateButton(hf,'Backup',[],'left+gapleft',hb);
hb.Callback='smac_GUI_synth(''backup'',gcbf)';
hb.TooltipString='Return to the SMAC autofit form';

%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

% Shrink the figure if it's bigger than the screen; don't keep proportion

pos=hf.Position(3:4);
if pos(1)>.98, pos(1)=.98; end
if pos(2)>.94, pos(2)=.94; end
hf.Position(3:4)=pos;

% Now make sure it's not off the screen

hf.Position(hf.Position<0.03)=0.03;

% Make the figure visible

hf.Visible='on';


%==========================================================================
%% uicbSynthResynthesize
function hf=uicbSynthResynthesize(varargin)
% Resynthesize based on mode shapes

global ss;

hf=varargin{1};
if ishandle(hf), delete(hf); end

% First clear out the analytical FRF

ss.fa=imat_fn([]);

% Open the GUI

smac_GUI_resynth;


%==========================================================================
%% uicbSynthCreateShapes
function hf=uicbSynthCreateShapes(varargin)
% Synthesize mode shapes

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

ind=get(hl,'Value');

diagtoggle=logical(get(getappdata(hf,'HandleShapeDiagnostics'),'Value'));

% If no roots have been selected don't do anything

if isempty(ind),
    uiwait(warndlg('Please select the roots you want to include'));
    return;
end

% Make sure none of the selected roots are condensed

if any(ss.fit.rootlist(ind,1)==0),
    uiwait(warndlg('You cannot perform synthesis using a 0 frequency root (condense first)'));
    return;
end

% Build the root list and other inputs to the synthesis

rootlist=[ind(:) ss.fit.rootlist(ind,:)];
freqrange=ss.freqrange(:,1);

% Select a shape fitting type for multi-reference fits

if length(ss.ref_coords)>1,
    types={'Max DP','Max CC','Best Synth'};
    btn=questdlg('How do you want to select which shapes to store?', ...
                 'Multi-reference shape selection',types{:},types{end});
    if isempty(btn), return; end
    fittype=find(strcmp(btn,types));
    if fittype<1, fittype=3; end
else
    fittype=3;
end

% Get the residual inputs

ht=getappdata(hf,'HandleResidualToggle');
if ht.Value
    residuals=ss.residuals;
else
    tmp=smac_create_data_structure;
    residuals=tmp.residuals;
end

% Generate shapes

hw=watchon;
drawnow;
[shp,refpick,psi]=smac_create_shapes(rootlist,residuals,fittype,diagtoggle);
ss.shape.shape=shp;
ss.shape.refpick=refpick;
ss.shape.psi=psi;
watchoff(hw);

if isempty(shp), return; end

% Update the status area on the form

uicbSynthUpdateStatus(hf);
uicbSynthUpdateList(hf);
set(hl,'Value',ind);


%==========================================================================
%% uicbSynthPlotShapes
function hf=uicbSynthPlotShapes(varargin)
% Plot mode shapes

global ss;

hf=varargin{1};

% See if we have any shapes to plot

if isempty(ss.shape.shape),
    uiwait(errordlg('Create the mode shapes first!'));
    return;
end

% Read in the FEM geometry if necessary

fem=ss.shape.fem;
if isempty(fem),
    fem=readunv('*.unv',[15 781 2411 780 2412 2420 82 2416 2431]);
    if isnumeric(fem), return; end
    ss.shape.fem=fem;
end

% Update the status area on the form

uicbSynthUpdateStatus(hf);

% Plot the shapes

plot(ss.shape.shape,fem);


%==========================================================================
%% uicbSynthSaveShapes
function hf=uicbSynthSaveShapes(varargin)
% Save mode shapes to ASH file

global ss;

hf=varargin{1};

% See if we have any shapes to plot

if isempty(ss.shape.shape),
    uiwait(errordlg('Create the mode shapes first!'));
    return;
end

% Get the ASH file to write to

[fname,pname]=uiputfile('*.ash','Select the ASH file to export shapes to');
if isnumeric(fname), return; end

fname=fullfile(pname,fname);
if isempty(findstr(fname,'.ash')),
    fname=[fname '.ash'];
end

% Delete the file if it exists

if exist(fname,'file'),
    delete(fname);
end

res=imat_shp([]);

% Build up the residuals as shapes for storage

if ss.residuals.use
    
    % Allocate space for the residual shapes
    
    res=imat_shp(6*length(ss.ref_coords)+1);
    nresshp=0;

    % First write out indexing shape for residuals
    
    refc=cellstr(ss.ref_coords);
    
    nresshp=nresshp+1;
    res(nresshp)=build_shape(allres(ss.fe(:,1)),(1:size(ss.fe,1))');
    res(nresshp).idline1='DOF indexing shape (shape coefficients are DOF order)';
    res(nresshp).idline2=['References: ' sprintf('%s ',refc{:})];
    res(nresshp).frequency=1000;

    % Low and high modes

    if ss.residuals.lowmode.active
        for k=1:length(ss.ref_coords)
            nresshp=nresshp+1;
            ctres=allres(ss.fe(:,k));
            nres=length(ctres);
            ind=nres*(k-1);
            res(nresshp)=build_shape(ctres,ss.residuals.lowmode.residual(ind+1:ind+nres).');
            res(nresshp).frequency=100+k;
            res(nresshp).idline1='Low mode (IDLine4=Frequency+Damping)';
            res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
            res(nresshp).idline3=num2str(ss.residuals.lowmode.frange);
            res(nresshp).idline4=num2str([ss.residuals.lowmode.freq ss.residuals.lowmode.damp]);
        end
    end

    if ss.residuals.highmode.active
        for k=1:length(ss.ref_coords)
            nresshp=nresshp+1;
            ctres=allres(ss.fe(:,k));
            nres=length(ctres);
            ind=nres*(k-1);
            res(nresshp)=build_shape(ctres,ss.residuals.highmode.residual(ind+1:ind+nres).');
            res(nresshp).frequency=200+k;
            res(nresshp).idline1='High mode (IDLine4=Frequency+Damping)';
            res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
            res(nresshp).idline3=num2str(ss.residuals.highmode.frange);
            res(nresshp).idline4=num2str([ss.residuals.highmode.freq ss.residuals.highmode.damp]);
        end
    end

%     % Mode inertance and compliance
% 
%     if ss.residuals.modeinertance.active
%         for k=1:length(ss.ref_coords)
%             nresshp=nresshp+1;
%             ctres=allres(ss.fe(:,k));
%             nres=length(ctres);
%             ind=nres*(k-1);
%             res(nresshp)=build_shape(ctres,ss.residuals.modeinertance.residual(ind+1:ind+nres).');
%             res(nresshp).frequency=300+k;
%             res(nresshp).idline1='Mode Inertance';
%             res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
%         end
%     end
% 
%     if ss.residuals.modecompliance.active
%         for k=1:length(ss.ref_coords)
%             nresshp=nresshp+1;
%             ctres=allres(ss.fe(:,k));
%             nres=length(ctres);
%             ind=nres*(k-1);
%             res(nresshp)=build_shape(ctres,ss.residuals.modecompliance.residual(ind+1:ind+nres).');
%             res(nresshp).frequency=400+k;
%             res(nresshp).idline1='Mode Compliance';
%             res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
%         end
%     end

    % General inertance and compliance

    if ss.residuals.inertance.active
        for k=1:length(ss.ref_coords)
            nresshp=nresshp+1;
            ctres=allres(ss.fe(:,k));
            nres=length(ctres);
            ind=nres*(k-1);
            res(nresshp)=build_shape(ctres,ss.residuals.inertance.residual(ind+1:ind+nres).');
            res(nresshp).frequency=500+k;
            res(nresshp).idline1='General Inertance (IDLine4=Frequency Range)';
            res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
            res(nresshp).idline4=num2str(ss.residuals.inertance.freq(:).');
        end
    end

    if ss.residuals.compliance.active
        for k=1:length(ss.ref_coords)
            nresshp=nresshp+1;
            ctres=allres(ss.fe(:,k));
            nres=length(ctres);
            ind=nres*(k-1);
            res(nresshp)=build_shape(ctres,ss.residuals.compliance.residual(ind+1:ind+nres).');
            res(nresshp).frequency=600+k;
            res(nresshp).idline1='General Compliance (IDLine4=Frequency Range)';
            res(nresshp).idline2=sprintf('Reference %d: %s',k,char(ss.ref_coords(k)));
            res(nresshp).idline4=num2str(ss.residuals.compliance.freq(:).');
        end
    end
    res=res(1:nresshp);

    % Tag as SMAC residual shapes
    
    res.damp=-1;
    res.idline5='SMAC Residual information';
    
end

% Write out the ADF

writeadf(fname,[ss.shape.shape(:); res]);


%==========================================================================
%% uicbSynthPlotMAC
function hf=uicbSynthPlotMAC(varargin)
% Plot mode shapes

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

% See if we have any shapes to fit

if isempty(ss.shape.shape),
    uiwait(errordlg('Create the mode shapes first!'));
    return;
end

% Get the indices to the selected shapes

ind=get(hl,'Value');

% If no roots have been selected don't do anything

if isempty(ind),
    uiwait(warndlg('Please select the roots you want to use in the MAC calculation'));
    return;
end

% Match the root indices to the fit mode shapes

ind=ind(:);
ind2=zeros(size(ind));
shfreq=ss.shape.shape.frequency;

for k=1:length(ind),
    [tmp,indt]=min(abs(shfreq-ss.fit.rootlist(ind(k))));
    if abs(tmp)<eps,
        ind2(k)=indt;
    end
end

if ~all(ind2),
    str=sprintf('Warning:  Shape not found for root %d (%g Hz)\n',[ind(~ind2) ss.fit.rootlist(ind(~ind2))]');
    uiwait(errordlg(str));
    if ~any(ind2), return; end
    ind2=ind2(ind2~=0);
end

% Calculate the MAC

o=ortho(ss.shape.shape(ind2));

% Plot the MAC

figure;
bar3o(o);

% Dress up the plot

axis tight
zlim([0 110]);
view(0,90);

set(gca,'XTick',1:length(shfreq),'YTick',1:length(shfreq));
set(gca,'XTickLabel',cellstr(num2str(shfreq(ind2),'%.3f')));
set(gca,'YTickLabel',cellstr(num2str(shfreq(ind2),'%.3f')));

ind=find(o>50);
if ~isempty(ind);
    [ii,jj]=ind2sub(size(o),ind);
    ht=text(ii,jj,o(ind)*1.05,cellstr(num2str(o(ind),'%.0f')));
    set(ht,'Color','green');
end


%==========================================================================
%% uicbSynthSelectAll
function hf=uicbSynthSelectAll(varargin)
% Select all of the roots

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

rootsel=1:size(get(hl,'String'),1);
ss.shape.rootsel=rootsel;

% Clean out data that is no longer valid

uicbSynthClearAnalysisData(hf);
uicbSynthUpdateList(hf);
set(hl,'Value',rootsel);


%==========================================================================
%% uicbSynthSelectNone
function hf=uicbSynthSelectNone(varargin)
% Deselect all of the roots

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

set(hl,'Value',[]);
ss.shape.rootsel=[];

% Clean out data that is no longer valid

uicbSynthClearAnalysisData(hf);
uicbSynthUpdateList(hf);


%==========================================================================
%% uicbSynthSelectFromList
function hf=uicbSynthSelectFromList(varargin)
% User selected roots from the listbox

global ss;
hf=varargin{1};

hl=getappdata(hf,'HandleList');
ss.shape.rootsel=hl.Value;

% Clean out data that is no longer valid

uicbSynthClearAnalysisData(hf);

% Update the list on the screen

uicbSynthUpdateList(hf);
hl.Value=ss.shape.rootsel;


%==========================================================================
%% uicbSynthAddRoots
function hf=uicbSynthAddRoots(varargin)
% Select all of the roots

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

% Open manual fit GUI

[out,cormax,refmax]=smac_GUI_addroot;

% If a new root was returned, squeeze it into the list

if ~isempty(out),
    rootlist=ss.fit.rootlist;
    rootlist(end+1,:)=[out refmax];
    rootlistorig=ss.fit.rootlistorig;
    rootlistorig(end+1,end)=0;
    corr_ref=ss.fit.corr_ref;
    corr_ref(end+1,:)=cormax;

    % Sort it to see where it belongs

    [tmp,ind]=sort(rootlist(:,1));
    ss.fit.rootlist=rootlist(ind,:);
    ss.fit.rootlistorig=rootlistorig(ind,:);
    ss.fit.corr_ref=corr_ref(ind,:);

    % Clean out data that is no longer valid

    uicbSynthClearAnalysisData(hf);

    % Update the list

    uicbSynthUpdateList(hf);

    % Preselect the new root in the list so we know which one was added

    rootsel=find(ind==length(ind));
    set(hl,'Value',rootsel);
    ss.shape.rootsel=rootsel;
end


%==========================================================================
%% uicbSynthDeleteRoots
function hf=uicbSynthDeleteRoots(varargin)
% Delete selected roots

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

ind=get(hl,'Value');

if ~isempty(ind),
    Btn=questdlg('Ok to delete selected root(s)?');
    if ~strcmp(Btn,'Yes'), return; end

end

ss.fit.rootlist(ind,:)=[];
ss.fit.rootlistorig(ind,:)=[];
ss.fit.corr_ref(ind,:)=[];

% Update the list

uicbSynthUpdateList(hf);

% Clean out data that is no longer valid

uicbSynthClearAnalysisData(hf);


%==========================================================================
%% uicbSynthCondenseRoots
function hf=uicbSynthCondenseRoots(varargin)
% Condense zeroed roots

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');
val=get(hl,'Value');

% Find the index of roots to get rid of

ind=ss.fit.rootlist(:,1)==0;

kind=false(size(ss.fit.rootlist(:,1),1),1);
kind(val)=true;
kind(ind)=[];

% Delete these roots

ss.fit.rootlist(ind,:)=[];
ss.fit.rootlistorig(ind,:)=[];
ss.fit.corr_ref(ind,:)=[];

valnew=find(kind);

% Update the list

uicbSynthUpdateList(hf);

% Select the roots previously selected that were not condensed

set(hl,'Value',valnew);
ss.shape.rootsel=valnew;

% Clean out data that is no longer valid

uicbSynthClearAnalysisData(hf);

% Update the list

uicbSynthUpdateList(hf);


%==========================================================================
%% uicbSynthUpdateList
function hf=uicbSynthUpdateList(varargin)
% Update the listbox

global ss

hf=varargin{1};
hl=getappdata(hf,'HandleList');
listfmt=getappdata(hf,'ListFormat');
listfmt2=getappdata(hf,'ListFormat2');

% Update the listbox contents

nroots=size(ss.fit.rootlist,1);
rl=ss.fit.rootlist;
rlo=ss.fit.rootlistorig;

str=cell(nroots,1);
for k=1:length(str)
    str{k}=sprintf(listfmt,[k rl(k,1) rl(k,2)*100 rl(k,3) rlo(k,1) rlo(k,2)*100 rlo(k,3)]);
end

% Now fill in the drive point information

defstr=sprintf(listfmt2,'','');
str2=repmat({defstr},nroots,1);
for k=1:length(ss.shape.refpick)
    %ref=ss.shape.refpick(k);
    %refc=ss.ref_coords(ref);
    refc=imat_ctrace(ss.shape.shape.referencecoord{k});
    dpcoef=ss.shape.shape{refc}(k);
    flag=' ';
    if dpcoef<0, flag='*'; end

    tempstr=sprintf('%10.3f',dpcoef);
    str2{ss.shape.rootsel(k)}=sprintf(listfmt2,char(refc),tempstr,flag);
end

% Append the new strings to the end of this list

for k=1:length(str)
    str{k}=[str{k} str2{k}];
end

% Set attributes

set(hl,'String',str);
set(hl,'Value',[]);
set(hl,'Max',length(ss.fit.rootlist));


%==========================================================================
%% uicbSynthSynthesize
function hf=uicbSynthSynthesize(varargin)
% Synthesize

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

ind=get(hl,'Value');

% If no roots have been selected don't do anything

if isempty(ind),
    uiwait(warndlg('Please select the roots you want to include in the synthesis'));
    return;
end

% Make sure none of the selected roots are condensed

if any(ss.fit.rootlist(ind,1)==0),
    uiwait(warndlg('You cannot perform synthesis using a 0 frequency root (condense first)'));
    return;
end

% If no synthesis options have been selected, return

if ~get(getappdata(hf,'ToggleNMIF'),'Value') && ...
        ~get(getappdata(hf,'ToggleMMIF'),'Value') && ...
        ~get(getappdata(hf,'ToggleCMIF'),'Value') && ...
        ~get(getappdata(hf,'ToggleFRF'),'Value'),
    uiwait(warndlg('Select at least one synthesis comparison!'));
    return;
end

% Build the root list and other inputs to the synthesis

rootlist=[ind(:) ss.fit.rootlist(ind,:)];
freqrange=ss.freqrangecc(:,1);

% Get the residual inputs

ht=getappdata(hf,'HandleResidualToggle');
if ht.Value
    residuals=ss.residuals;
else
    tmp=smac_create_data_structure;
    residuals=tmp.residuals;
end

% Perform the synthesis

hw=watchon;
drawnow;
if isempty(ss.fa),
    [frfsyn,residues,residuals]=smac_synthesize(rootlist,residuals);
    ss.fa=frfsyn;
    ss.shape.residues=residues;
    if ss.residuals.use
        ss.residuals=residuals;
    end

    % Update the status area on the form

    uicbSynthUpdateStatus(hf);
else
    frfsyn=ss.fa;
end

%--------------------------------------------------------------------------
% Plot NMIF if requested

h_nmif=getappdata(hf,'Figure_Handle_NMIF');
if get(getappdata(hf,'ToggleNMIF'),'Value'),
    h_nmif=smac_synth_nmif(h_nmif,frfsyn,ss.fe,rootlist,freqrange);
    setappdata(h_nmif.hf,'FigureType','NMIF');
    setappdata(hf,'Figure_Handle_NMIF',h_nmif);
else

    % Delete figure if we are not plotting it

    if ~isempty(h_nmif) && ishandle(h_nmif.hf),
        if strcmp(getappdata(h_nmif.hf,'FigureType'),'NMIF')
            delete(h_nmif.hf);
            setappdata(hf,'Figure_Handle_NMIF',[]);
        end
    end
end


%--------------------------------------------------------------------------
% Plot CMIF if requested

h_cmif=getappdata(hf,'Figure_Handle_CMIF'); 

if get(getappdata(hf,'ToggleCMIF'),'Value'),
    h_cmif=smac_synth_cmif(h_cmif,frfsyn,ss.fe,rootlist,freqrange);
    setappdata(h_cmif.hf,'FigureType','CMIF');
    setappdata(hf,'Figure_Handle_CMIF',h_cmif);
else

    % Delete figure if we are not plotting it

    if ~isempty(h_cmif) && ishandle(h_cmif.hf),
        if strcmp(getappdata(h_cmif.hf,'FigureType'),'CMIF')
            delete(h_cmif.hf);
            setappdata(hf,'Figure_Handle_CMIF',[]);
        end
    end
end


%--------------------------------------------------------------------------
% Plot MMIF if requested

h_mmif=getappdata(hf,'Figure_Handle_MMIF');
if get(getappdata(hf,'ToggleMMIF'),'Value'),
    h_mmif=smac_synth_mmif(h_mmif,frfsyn,ss.fe,rootlist,freqrange);
    setappdata(h_mmif.hf,'FigureType','MMIF');
    setappdata(hf,'Figure_Handle_MMIF',h_mmif);
else

    % Delete figure if we are not plotting it

    if ~isempty(h_mmif) && ishandle(h_mmif.hf),
        if strcmp(getappdata(h_mmif.hf,'FigureType'),'MMIF')
            delete(h_mmif.hf);
            setappdata(hf,'Figure_Handle_MMIF',[]);
        end
    end
end


%--------------------------------------------------------------------------
% Plot synthesized FRF compared to experimental if requested

h_frf=getappdata(hf,'Figure_Handle_FRF');
if get(getappdata(hf,'ToggleFRF'),'Value'),
    h_frf=smac_synth_frf(h_frf,frfsyn,ss.fe,rootlist,freqrange,ss.realcomplex,ss.residuals.use);
    setappdata(h_frf.hf,'FigureType','FRF');
    setappdata(hf,'Figure_Handle_FRF',h_frf);
else

    % Delete figure if we are not plotting it

    if ~isempty(h_frf) && ishandle(h_frf.hf),
        if strcmp(getappdata(h_frf.hf,'FigureType'),'FRF')
            delete(h_frf.hf);
            setappdata(hf,'Figure_Handle_FRF',[]);
        end
    end
end

watchoff(hw);


%==========================================================================
%% uicbSynthCloseSynthPlots
function hf=uicbSynthCloseSynthPlots(varargin)
% Close synthesis plots

hf=varargin{1};

h=getappdata(hf,'Figure_Handle_NMIF');
if ~isempty(h) && ishandle(h.hf), delete(h.hf); end
setappdata(hf,'Figure_Handle_NMIF',[]);
h=getappdata(hf,'Figure_Handle_CMIF');
if ~isempty(h) && ishandle(h.hf), delete(h.hf); end
setappdata(hf,'Figure_Handle_CMIF',[]);
h=getappdata(hf,'Figure_Handle_MMIF');
if ~isempty(h) && ishandle(h.hf), delete(h.hf); end
setappdata(hf,'Figure_Handle_MMIF',[]);
h=getappdata(hf,'Figure_Handle_FRF');
if ~isempty(h) && ishandle(h.hf), delete(h.hf); end
setappdata(hf,'Figure_Handle_FRF',[]);


%==========================================================================
%% uicbSynthUpdateStatus
function uicbSynthUpdateStatus(varargin)
% Updates the status area on the form

global ss;

hf=varargin{1};
h_fem=getappdata(hf,'StatusFEM');
h_shp=getappdata(hf,'StatusShapes');
h_syn=getappdata(hf,'StatusSynthesis');

cdata_no=zeros(17,17,3);
cdata_no(:,:,1)=1;
cdata_yes=zeros(17,17,3);
cdata_yes(:,:,2)=0.67;

% FEM

if isempty(ss.shape.fem),
    set(h_fem,'CData',cdata_no);
else
    set(h_fem,'CData',cdata_yes);
end

% Shapes

if isempty(ss.shape.shape),
    set(h_shp,'CData',cdata_no);
else
    set(h_shp,'CData',cdata_yes);
end

% Synthesized FRF

if isempty(ss.fa),
    set(h_syn,'CData',cdata_no);
else
    set(h_syn,'CData',cdata_yes);
end


%==========================================================================
%% uicbSynthResiduals
function uicbSynthResiduals(varargin)
% Get residual parameters

global ss;

hf=varargin{1};
hr=getappdata(hf,'HandleResiduals');

out=smac_GUI_residuals([],handle(hr));
if isempty(out)
    return;
else
    if ~isequal(ss.residuals,out)
        uicbSynthClearAnalysisData(hf);
    end
    ss.residuals=out;
end


%==========================================================================
%% uicbSynthToggleResiduals
function uicbSynthToggleResiduals(varargin)
% Toggle whether residual parameters will be added

global ss

hf=varargin{1};

uicbSynthClearAnalysisData(varargin{:});

hc=getappdata(hf,'HandleResidualToggle');
ss.residuals.use=logical(get(hc,'Value'));

ss.residuals.inertance.residual=[];
ss.residuals.compliance.residual=[];

% ss.residuals.modeinertance.residual=[];
% ss.residuals.modecompliance.residual=[];

ss.residuals.lowmode.residual=[];
ss.residuals.highmode.residual=[];

%==========================================================================
%% uicbSynthCMIFtoggle
function uicbSynthCMIFtoggle(varargin)
% Toggle whether residual parameters will be added

hf = varargin{1};

toggleCMIF = getappdata(hf,'ToggleCMIF');
ToggleTracking = getappdata(hf,'ToggleTracking');


%==========================================================================
%% uicbSynthUpdateTracking
function uicbSynthUpdateTracking(varargin)
% Save Tracking variable when altered

global ss;
hf = varargin{1};
ToggleTracking = getappdata(hf,'ToggleTracking');
ss.mif.ana.reftracking = get(ToggleTracking,'Value');
ss.mif.exp.reftracking = get(ToggleTracking,'Value');

% Blank out the CMIF functions since they will need to be recalculated

ss.mif.ana.cmif = [];
ss.mif.exp.cmif = [];


%==========================================================================
%% uicbSynthClearAnalysisData
function uicbSynthClearAnalysisData(varargin)
% Delete all synthesized data from the structure

global ss;

hf=varargin{1};

% Original functions

ss.fa=imat_fn([]);

% MIFs

ss.mif.ana.nmif=imat_fn([]);
ss.mif.ana.cmif=imat_fn([]);
ss.mif.ana.mmif=imat_fn([]);

% Shapes

ss.shape.shape=imat_shp([]);
ss.shape.residues=[];
ss.shape.refpick=[];

% Update the status area on the form

uicbSynthUpdateStatus(hf);


%==========================================================================
%% uicbSynthSave
function hf=uicbSynthSave(varargin)
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
%% uicbSynthBackup
function hf=uicbSynthBackup(varargin)
% Backup

hf=varargin{1};

if ishandle(hf), delete(hf); end
smac_setprogress('selroots');
smac_GUI_autofit;


%==========================================================================
%% uicbSynthQuit
function hf=uicbSynthQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strmatch(rslt,'Yes'),

    % Get rid of synthesis figures

    uicbSynthCloseSynthPlots(hf);

    % Delete the form

    delete(hf);

end
