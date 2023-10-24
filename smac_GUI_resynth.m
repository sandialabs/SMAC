function out=smac_GUI_resynth(action,varargin)
% SMAC_GUI_RESYNTH  Resynthesis based on shapes and residuals
%
% OUT=SMAC_GUI_RESYNTH
%
% SMAC_GUI_RESYNTH allows the user to synthesize FRF (and then NMIF, CMIF,
% and MMIF) based on the curve-fit mode shapes.  The residual effects
% applied during the synthesis can be applied here as well.  The parameters
% cannot be changed, but the individual effects can be toggled on and off.
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
%  15-Jun-2005 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Update residual information shown on the form
%
%  07-Mar-2006 / ATA Engineering / Dan Hensley
%    o Add toggle to recalculate residual terms
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
    action='';
end

% Register valid callbacks

bname='uicbResynth';

funcs.plot='Plot';
funcs.close_plots='ClosePlots';

funcs.shapes_all='SelectAll';
funcs.shapes_none='SelectNone';
funcs.shapes_select='SelectShapes';

funcs.toggleresiduals='ToggleResiduals';
funcs.toggleresidualrecalc='ToggleResidualRecalc';

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
        error(sprintf('Unknown action ''%s''',action));
    end
end

% Build the GUI

hf=uicbResynthCreateFigure(varargin);


%==========================================================================
%% CALLBACKS
%==========================================================================
%% uicbResynthCreateFigure
function hf=uicbResynthCreateFigure(varargin)
% Create the Re-Synthesis form

global ss;

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Re-Synthesis', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_ReSynthForm', ...
    'Visible','off');
hf=handle(hf);

% First look through the residual structure and turn off if necessary

residuals=ss.residuals;
if ~residuals.use
    fn=fieldnames(residuals);
    fn(strcmp(fn,'use'))=[];
    for k=1:length(fn)
        residuals.(fn{k}).active=false;
    end
end

setappdata(hf,'Residuals',residuals);

%--------------------------------------------------------------------------
% Add Plot button

hb=ui.CreateButton(hf,'Plot',[],uip.Gap*[1 1]);
hb.Callback='smac_GUI_resynth(''plot'',gcbf)';
hb.TooltipString='Synthesize the FRF and plot the selected functions against the experimental data';

hb2=ui.CreateButton(hf,'Close Plots',[],'right+gapright',hb);
hb2.Callback='smac_GUI_resynth(''close_plots'',gcbf)';
hb2.TooltipString='Close any open plots';
ui.FitText(hb2,[0 inf]);

%--------------------------------------------------------------------------
% Add Residual information panel and options

hfrs=ui.CreatePanel(hf,' Residual Information',[400 100],'top',hb);

hfrs.Position(2)=hfrs.Position(2)+20;

%----------------------------------
% Add Residual options toggles and disable if residuals were not applied

clear hc ht
enable={'off','on'};

enabletog=@(x) enable{double(x)+1};

% Recalculate residuals?

anyres=residuals.compliance.active || ...
       residuals.inertance.active || ...
       residuals.lowmode.active || ...
       residuals.highmode.active;

hc(1)=ui.CreateCheckbox(hfrs,'Recalculate residuals',[],[uip.Gap uip.Gap]);
hc(end).Callback='smac_GUI_resynth(''toggleresidualrecalc'',gcbf)';
hc(end).Tag='Recalculate residuals instead of using previously stored values';
hc(end).Enable=enabletog(anyres);
hc(end).Value=length(ss.ref_coords)>1;
ui.FitText(hc(end),[0 inf]);
setappdata(hf,'RecalcResiduals',hc(1));

% Horizontal line

hh=ui.CreateHorizLine(hfrs,'HorizLine',100,'top+gaptop',hc(end));

% General Compliance

hc(end+1)=ui.CreateCheckbox(hfrs,'General Compliance',[],'top+gaptop',hh);
hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
hc(end).Tag='Compliance';
hc(end).Enable=enabletog(residuals.compliance.active);
hc(end).Value=residuals.compliance.active;
ui.FitText(hc(end),[0 inf]);

ht(1)=ui.CreateText(hfrs,' ',[],'right+gapright',hc(end));
if residuals.compliance.active
    ht(1).String=sprintf('Frequency range: %g-%g Hz',residuals.compliance.freq);
end
ui.FitText(ht(1),[0 inf]);

% General Inertance

hc(end+1)=ui.CreateCheckbox(hfrs,'General Inertance',[],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
hc(end).Tag='Inertance';
hc(end).Enable=enabletog(residuals.inertance.active);
hc(end).Value=residuals.inertance.active;
ui.FitText(hc(end),[0 inf]);

ht(end+1)=ui.CreateText(hfrs,' ',[],'right+gapright',hc(end));
if residuals.inertance.active
    ht(end).String=sprintf('Frequency range: %g-%g Hz',residuals.inertance.freq);
end
ui.FitText(ht(end),[0 inf]);

% % Mode Compliance
% 
% hc(end+1)=ui.CreateCheckbox(hfrs,'Mode Compliance',[],'top+gaptop',hc(end));
% hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
% hc(end).Tag='ModeCompliance';
% hc(end).Enable=enable{double(residuals.modecompliance.active+1)};
% hc(end).Value=residuals.modecompliance.active;
% ui.FitText(hc(end),[0 inf]);
% 
% % Mode Inertance
% 
% hc(end+1)=ui.CreateCheckbox(hfrs,'Mode Inertance',[],'top+gaptop',hc(end));
% hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
% hc(end).Tag='ModeInertance';
% hc(end).Enable=enable{double(residuals.modeinertance.active+1)};
% hc(end).Value=residuals.modeinertance.active;
% ui.FitText(hc(end),[0 inf]);

% High Mode

hc(end+1)=ui.CreateCheckbox(hfrs,'High Mode',[],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
hc(end).Tag='HighMode';
hc(end).Enable=enabletog(residuals.highmode.active);
hc(end).Value=residuals.highmode.active;
ui.FitText(hc(end),[0 inf]);

ht(end+1)=ui.CreateText(hfrs,' ',[],'right+gapright',hc(end));
if residuals.highmode.active
    ht(end).String=sprintf('%g Hz, %g%% (fit range %g-%g Hz)', ...
                         residuals.highmode.freq, ...
                         residuals.highmode.damp*100, ...
                         residuals.highmode.frange);
end
ui.FitText(ht(end),[0 inf]);

% Low Mode

hc(end+1)=ui.CreateCheckbox(hfrs,'Low Mode',[],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_resynth(''toggleresiduals'',gcbf,gcbo)';
hc(end).Tag='LowMode';
hc(end).Enable=enabletog(residuals.lowmode.active);
hc(end).Value=residuals.lowmode.active;
ui.FitText(hc(end),[0 inf]);

ht(end+1)=ui.CreateText(hfrs,' ',[],'right+gapright',hc(end));
if residuals.lowmode.active
    ht(end).String=sprintf('%g Hz, %g%% (fit range %g-%g Hz)', ...
                         residuals.lowmode.freq, ...
                         residuals.lowmode.damp*100, ...
                         residuals.lowmode.frange);
end
ui.FitText(ht(end),[0 inf]);

% Move all of the text widgets over and resize the frame

maxw=0;
for k=1:length(hc)
    maxw=max(hc(k).Position(3)+hc(k).Position(1)+2*uip.Gap,maxw);
end
for k=1:length(ht)
    ht(k).Position(1)=ht(k).Position(1)+(maxw-ht(k).Position(1));
end

ui.ResizeControl(hfrs);

% Resize the horizontal line

hh.Position(3)=hfrs.Position(3)-hfrs.Position(1)-uip.Gap;

%--------------------------------------------------------------------------
% Add Synthesis Options panel and options

hfr=ui.CreatePanel(hf,' Synthesis Options',[400 100],'top+gaptop+gaptop',hfrs);

hc1=ui.CreateCheckbox(hfr,'FRF',[],[uip.Gap uip.Gap]);
hc1.TooltipString='Synthesize frequency response functions';
hc1.Value=1;
setappdata(hf,'ToggleFRF',hc1);
ui.FitText(hc1);

hc=ui.CreateCheckbox(hfr,'CMIF',[],'top+gaptop',hc1);
hc.TooltipString='Synthesize complex mode indicator functions';
if ss.realcomplex==2, hc.Value=1; end
setappdata(hf,'ToggleCMIF',hc);
ui.FitText(hc);

hc2=ui.CreateCheckbox(hfr,'Ref tracking',[],'right+gapright+gapright',hc);
hc2.TooltipString='Synthesize CMIF using reference tracking';
hc2.Value=double(ss.mif.ana.reftracking);
setappdata(hf,'CMIFRefTracking',hc2);
ui.FitText(hc2);

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

% Resize the controls

ui.ResizeControl(hfr);
hfr.Position(3)=hfrs.Position(3);

%--------------------------------------------------------------------------
% Add Shapes listbox and title

colwidth=[4 10 9 10 35];

listfmt=sprintf('%% %dd  %%%d.3f  %%%d.4f  %%-%ds  %%-%ds',colwidth);
titlfmt=sprintf(sprintf('%%%ds  %%%ds  %%%ds  %%-%ds  %%-%ds',colwidth), ...
                'Root','Freq (Hz)','Damp (%)','Type','IDLine1');

hl=ui.CreateList(hf,sprintf(listfmt,[1 0 0],'a','a'),[],'right+gapright+alignbottom',hfrs);
hl.Callback='smac_GUI_resynth(''shapes_select'',gcbf)';
hl.FontSize=10;
hl.Max=3;
hl.Position(4)=hfr.Position(2)+hfr.Position(4)-hl.Position(2);
ui.FitText(hl,[0 inf]);	% Only fit non-infinite dimension

setappdata(hf,'HandleList',hl);

ht=ui.CreateText(hf,titlfmt,[],'top',hl);
pos=ht.Position; pos(1)=pos(1)+4; ht.Position=pos;
ht.FontName='Courier New';
ht.FontSize=10;
ui.FitText(ht);

% Fill in the text

shp=ss.shape.shape;
str=char(32*ones(length(shp),length(ht.String)));
for k=1:prod(size(shp))
    idline1=shp(k).idline1;
    if iscell(idline1), idline1=idline1{:}; end
    idline1=idline1(1:min(length(idline1),colwidth(end)));
    str(k,:)=sprintf(listfmt,k,shp(k).frequency,shp(k).damping,shp(k).shapetype,idline1);
end
hl.String=str;
hl.Value=1:size(hl.String,1);

%--------------------------------------------------------------------------
% Add All and None buttons

hb=ui.CreateButton(hf,'All',[],'aligntop+right+gapright',hl);
hb.Callback='smac_GUI_resynth(''shapes_all'',gcbf)';
hb.TooltipString='Select all shapes';

hb=ui.CreateButton(hf,'None',[],'bottom+gapbottom',hb);
hb.Callback='smac_GUI_resynth(''shapes_none'',gcbf)';
hb.TooltipString='Deselect all shapes';

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit,backup)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',hl);
hb.Callback='smac_GUI_resynth(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';
hb.Position(2)=uip.Gap;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''res_hlp'')';
hb.TooltipString='Get help on the Resynthesis calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_resynth(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';

% Backup button

if ~isempty(ss.pinv)
    hb=ui.CreateButton(hf,'Backup',[],'left+gapleft',hb);
    hb.Callback='smac_GUI_resynth(''backup'',gcbf)';
    hb.TooltipString='Return to the SMAC synthesis form';
end

%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';


%==========================================================================
%% uicbResynthPlot
function hf=uicbResynthPlot(varargin)
% Synthesize and plot

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');
ind=get(hl,'Value');

% If no shapes have been selected don't do anything

if isempty(ind),
    uiwait(warndlg('Please select the shapes you want to include in the synthesis'));
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

% Compute the analytical functions

hw=watchon;
drawnow;

if ~isempty(ss.fa)
    frfsyn=ss.fa;
else
    residuals=getappdata(hf,'Residuals');
    recalc=logical(get(getappdata(hf,'RecalcResiduals'),'Value'));
    shp=ss.shape.shape(ind);

    frfsyn=smac_resynthesize(ss.fe,shp,residuals,recalc);
    watchoff(hw);

    if isempty(frfsyn),
        uiwait(errordlg('Error generating synthesized functions'));
        return;
    end

    ss.fa=frfsyn;
end

rootlist=[];
freqrange=[];

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
if get(getappdata(hf,'ToggleCMIF'),'Value')
    % Toggle the reference tracking
    hc=getappdata(hf,'CMIFRefTracking');
    ss.mif.ana.reftracking = get(hc,'Value');
    ss.mif.exp.reftracking = get(hc,'Value');

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
%% uicbResynthToggleResiduals
function uicbResynthToggleResiduals(varargin)
% Toggle whether residual parameters will be added

global ss

hf=varargin{1};

residuals=getappdata(hf,'Residuals');
tag=get(gcbo,'Tag');

residuals.(lower(tag)).active=get(gcbo,'Value');

setappdata(hf,'Residuals',residuals);

uicbResynthClearAnalysisData(varargin{:});


%==========================================================================
%% uicbResynthToggleResidualRecalc
function uicbResynthToggleResidualRecalc(varargin)
% Clear out synthesize functions

uicbResynthClearAnalysisData(varargin{:});


%==========================================================================
%% uicbResynthClosePlots
function hf=uicbResynthClosePlots(varargin)
% Close synthesis plots

global ss;

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
%% uicbResynthClearAnalysisData
function uicbResynthClearAnalysisData(varargin)
% Delete all synthesized data from the structure

global ss;

hf=varargin{1};

% Original functions

ss.fa=imat_fn([]);

% MIFs

ss.mif.ana.nmif=imat_fn([]);
ss.mif.ana.cmif=imat_fn([]);
ss.mif.ana.mmif=imat_fn([]);


%==========================================================================
%% uicbResynthSelectAll
function hf=uicbResynthSelectAll(varargin)
% Select all of the shapes

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

rootsel=1:size(get(hl,'String'),1);
set(hl,'Value',rootsel);

% Clean out data that is no longer valid

uicbResynthClearAnalysisData(hf);


%==========================================================================
%% uicbResynthSelectNone
function hf=uicbResynthSelectNone(varargin)
% Deselect all of the shapes

global ss;

hf=varargin{1};
hl=getappdata(hf,'HandleList');

set(hl,'Value',[]);

% Clean out data that is no longer valid

uicbResynthClearAnalysisData(hf);


%==========================================================================
%% uicbResynthSelectShapes
function hf=uicbResynthSelectShapes(varargin)
% Shapes were selected from the list

hf=varargin{1};

% Clean out data that is no longer valid

uicbResynthClearAnalysisData(hf);


%==========================================================================
%% uicbResynthSave
function hf=uicbResynthSave(varargin)
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
%% uicbResynthBackup
function hf=uicbResynthBackup(varargin)
% Backup

hf=varargin{1};

if ishandle(hf), delete(hf); end

% Delete the synthesized FRF since they are no longer valid

ss.fa=imat_fn([]);

% Open the GUI

smac_GUI_synth;


%==========================================================================
%% uicbResynthQuit
function hf=uicbResynthQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strcmpi(rslt,'Yes'),

    % Get rid of synthesis figures

    uicbResynthClosePlots(hf);

    % Delete the form

    delete(hf);

end
