function out=smac_GUI_autofit(action,varargin)
% SMAC_GUI_AUTOFIT  Set up and perform auto fit
%
% OUT=SMAC_GUI_AUTOFIT
%
% SMAC_GUI_AUTOFIT allows the user to select a damping and frequency range
% over which to iterate when performing the automatic peak value detection
% around the initially supplied peaks.  The user can also specify how many
% frequency and damping values are used in the initial optimization.
%
% OUT will be 1 if the user selects the Execute Auto SMAC button and the
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
%    o Changed wording on some controls
%    o Added damping and frequency convergence edit boxes
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Add toggle for convergence plots during autofit
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make GUI non-blocking (don't use uiwait/uiresume), and change
%      forward/backup calls to call the next/previous GUI explicitly
%
%  09-Nov-2004 / ATA Engineering / Dan Hensley
%    o Set progress when user hits Backup button
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Tweaks so this runs properly with Matlab 7.0.x and is still
%      compatible with Matlab 6.5.x
%
%  07-Jun-2005 / ATA Engineering / Dan Hensley
%    o Use new GUI functions in uiWidgets
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
    action='';
end

% Register valid callbacks

bname='uicbAutoFit';
funcs.execute='Execute';

funcs.numdampedit='NumDampEdit';
funcs.dampconvedit='DampConvEdit';
funcs.lowdampedit='LowDampEdit';
funcs.highdampedit='HighDampEdit';
funcs.numfreqedit='NumFreqEdit';
funcs.freqconvedit='FreqConvEdit';
funcs.freqedit='FreqEdit';

funcs.plotfrf='PlotFRF';

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

uicbAutoFitCreateFigure(varargin);


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbAutoFitCreateFigure(varargin)
% Create the SMAC Auto Fit form

global ss;

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Auto Fit', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_AutoFitForm', ...
    'Visible','off');
hf=handle(hf);
% setappdata(hf,'Backup',false);

%--------------------------------------------------------------------------
% Add Execute button

hb=ui.CreateButton(hf,'Execute Auto SMAC',[170 40],uip.Gap*[1 1]);
hb.Callback='smac_GUI_autofit(''execute'',gcbf)';
hb.TooltipString='Perform the curve-fit on the initial roots selected';
%uiFitText(hb);

% Convergence toggle button

ht=ui.CreateCheckbox(hf,'Display convergence graphs',[inf 40],'top+gaptop',hb);
ht.TooltipString='Display convergence graphs during autofit';
ht.Value=0;
setappdata(hf,'GraphToggle',ht);
%uiFitText(ht);

%--------------------------------------------------------------------------
% Add panel for damping parameters information

hfr=ui.CreatePanel(hf,' Damping range',[400 100],'top+gaptop',ht);

%----------------------------------
% Add Damping range information

hr=ui.CreateText(hfr,'Damping convergence percentage',[],[uip.Gap uip.Gap]);
hr.TooltipString=sprintf(['Specifies the damping convergence tolerance used in the autofit routine.\n', ...
                          'PLEASE USE CAUTION WHEN CHANGING THIS VALUE']);
hr.HorizontalAlignment='right';

he=ui.CreateEdit(hfr,num2str(ss.fit.zetaconv*100),[],'right+gapright',hr);
he.Callback='smac_GUI_autofit(''dampconvedit'',gcbf,gcbo)';
he.TooltipString=sprintf('Enter the damping convergence percentage');
setappdata(hf,'HandleEditDampConv',he);
setappdata(hf,'DampConvOld',ss.fit.zetaconv*100);

%

he=ui.CreateEdit(hfr,num2str(ss.fit.zetapts),[],'top+gaptop',he);
he.Callback='smac_GUI_autofit(''numdampedit'',gcbf,gcbo)';
he.TooltipString='Enter the number of damping values to use in the iteration';
setappdata(hf,'HandleEditNumDamp',he);
setappdata(hf,'NumDampOld',ss.fit.zetapts);

hr=ui.CreateText(hfr,'Number of values in each range',[],'left+gapleft',he);
hr.TooltipString='Use this number of damping values in the iteration';
hr.HorizontalAlignment='right';

%

he=ui.CreateEdit(hfr,num2str(ss.fit.zeta(1)*100),[],'top+gaptop',he);
he.Callback='smac_GUI_autofit(''lowdampedit'',gcbf,gcbo)';
he.TooltipString='Enter the low damping value in %';
setappdata(hf,'DampLowOld',ss.fit.zeta(1)*100);

hr=ui.CreateText(hfr,'Damping Range (%)',[],'left+gapleft',he);
hr.TooltipString='Fit this range of damping values to determine the optimal value';
hr.HorizontalAlignment='right';

ht=ui.CreateText(hfr,'-',[],'right+gapright',he);

he=ui.CreateEdit(hfr,num2str(ss.fit.zeta(2)*100),[],'right+gapright',ht);
he.Callback='smac_GUI_autofit(''highdampedit'',gcbf,gcbo)';
he.TooltipString='Enter the high damping value in %';
setappdata(hf,'DampHighOld',ss.fit.zeta(2)*100);

ui.ResizeControl(hfr);

%--------------------------------------------------------------------------
% Add panel for frequency parameters information

hfr2=ui.CreatePanel(hf,' Frequency Range',[400 100],'top+gaptop+gaptop',hfr);

%----------------------------------
% Add Frequency Lines information

hr=ui.CreateText(hfr2,'Frequency convergence percentage',[],[uip.Gap uip.Gap]);
hr.TooltipString=sprintf(['Specifies the frequency convergence tolerance used in the autofit routine.\n', ...
                          'PLEASE USE CAUTION WHEN CHANGING THIS VALUE']);
hr.HorizontalAlignment='right';

he=ui.CreateEdit(hfr2,num2str(ss.fit.freqconv*100),[],'right+gapright',hr);
he.Callback='smac_GUI_autofit(''freqconvedit'',gcbf,gcbo)';
he.TooltipString=sprintf('Enter the frequency convergence percentage');
setappdata(hf,'HandleEditFreqConv',he);
setappdata(hf,'FreqConvOld',ss.fit.freqconv*100);

%

he=ui.CreateEdit(hfr2,num2str(ss.fit.freqpts),[],'top+gaptop',he);
he.Callback='smac_GUI_autofit(''numfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter the number of frequency values to use in the iteration';
setappdata(hf,'HandleEditNumFreq',he);
setappdata(hf,'NumFreqOld',ss.fit.freqpts);

hr=ui.CreateText(hfr2,'Number of values in each range',[],'left+gapleft',he);
hr.TooltipString='Use this number of frequency values in the iteration';
hr.HorizontalAlignment='right';

%

he=ui.CreateEdit(hfr2,num2str(ss.fit.freqp),[],'top+gaptop',he);
he.Callback='smac_GUI_autofit(''freqedit'',gcbf,gcbo)';
he.TooltipString='Enter the frequency range in %';
setappdata(hf,'HandleEditFreq',he);
setappdata(hf,'FreqRangeOld',ss.fit.freqp);

hr=ui.CreateText(hfr2,'Frequency Range (%)',[],'left+gapleft',he);
hr.TooltipString='Fit this range of frequency values to determine the optimal value';
hr.HorizontalAlignment='right';

ui.ResizeControl(hfr2);

% Now make sure the two frames are the same width

pos=[hfr.Position; hfr2.Position];
pos(:,3)=max(pos(:,3));
hfr.Position=pos(1,:);
hfr2.Position=pos(2,:);

%--------------------------------------------------------------------------
% Add plot FRF button

hb=ui.CreateButton(hf,'Plot FRFs',[100 40],'aligntop+right+gapright+gapright',hfr);
hb.Callback='smac_GUI_autofit(''plotfrf'')';
hb.TooltipString='Plot experimental FRF';

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',hb);
hb.Callback='smac_GUI_autofit(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';

% Move this button down

hb.Position(2)=uip.Gap;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''sf_hlp'')';
hb.TooltipString='Get help on the Correlation Coefficient calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_autofit(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';

% Backup button

hb=ui.CreateButton(hf,'Backup',[],'left+gapleft',hb);
hb.Callback='smac_GUI_autofit(''backup'',gcbf)';
hb.TooltipString='Return to the correlation coefficient plot form';


%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';


%==========================================================================
function hf=uicbAutoFitExecute(varargin)
% Calculate correlation coefficients

global ss;

hf=watchon;
drawnow;
ht=getappdata(hf,'GraphToggle');
dograph=logical(get(ht,'Value'));

% Perform the curve-fit

smac_autofit(dograph);

smac_setprogress('autofit');

watchoff(hf);

% Delete the figure and move to the next one

if ishandle(hf), delete(hf); end
smac_GUI_synth;


%==========================================================================
function hf=uicbAutoFitNumDampEdit(varargin)
% Make sure number of damping values is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the damping value from the edit box

nd=str2num(get(he,'String'));
if isempty(nd),
    uiwait(errordlg('Number of damping values must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'NumDampOld')));
    return;
end

nd=fix(nd(1));

% Check it for validity
if nd<1,
    uiwait(errordlg('Invalid number of damping values to use'));
    nd=ss.fit.zetapts;
end

% Update the value and the form

ss.fit.zetapts=nd;
set(he,'String',num2str(nd));
setappdata(hf,'NumDampOld',nd);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitDampConvEdit(varargin)
% Make sure damping convergence is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the damping value from the edit box

dc=str2num(get(he,'String'));
if isempty(dc),
    uiwait(errordlg('Damping convergence fraction must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'DampConvOld')));
    return;
end

dc=dc(1);

% Check it for validity
if dc<=0 || dc>10
    uiwait(errordlg('Damping convergence must be between 0 and 10%'));
    dc=ss.fit.zetaconv*100;
end

% Update the value and the form

ss.fit.zetaconv=dc/100;
set(he,'String',num2str(dc));
setappdata(hf,'DampConvOld',dc);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitLowDampEdit(varargin)
% Make sure low damping value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the damping value from the edit box

damp=str2num(get(he,'String'));
if isempty(damp),
    uiwait(errordlg('Damping value must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'DampLowOld')));
    return;
end

damp=damp(1);

% Check it for validity

if damp<=0 || damp>25
    uiwait(errordlg('Invalid damping value (0<zeta<=25%)'));
    damp=ss.fit.zeta(1)*100;
end
if damp>ss.fit.zeta(2)*100,
    uiwait(errordlg('The low damping value must be less than the high value'));
    damp=ss.fit.zeta(1)*100;
end

% Update the value and the form

ss.fit.zeta(1)=damp/100;
set(he,'String',num2str(damp));
setappdata(hf,'DampLowOld',damp);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitHighDampEdit(varargin)
% Make sure high damping value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the damping value from the edit box

damp=str2num(get(he,'String'));
if isempty(damp),
    uiwait(errordlg('Damping value must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'DampHighOld')));
    return;
end

damp=damp(1);

% Check it for validity

if damp<=0 || damp>25
    uiwait(errordlg('Invalid damping value (0<zeta<=25%)'));
    damp=ss.fit.zeta(2)*100;
end
if damp<ss.fit.zeta(1)*100
    uiwait(errordlg('The high damping value must be greater than the low value'));
    damp=ss.fit.zeta(2)*100;
end

% Update the value and the form

ss.fit.zeta(2)=damp/100;
set(he,'String',num2str(damp));
setappdata(hf,'DampHighOld',damp);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitFreqConvEdit(varargin)
% Make sure frequency convergence is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

fc=str2num(get(he,'String'));
if isempty(fc),
    uiwait(errordlg('Frequency convergence fraction must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'FreqConvOld')));
    return;
end

fc=fc(1);

% Check it for validity
if fc<=0 || fc>2
    uiwait(errordlg('Frequency convergence must be between 0 and 2%'));
    fc=ss.fit.freqconv*100;
end

% Update the value and the form

ss.fit.freqconv=fc/100;
set(he,'String',num2str(fc));
setappdata(hf,'FreqConvOld',fc);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitNumFreqEdit(varargin)
% Make sure number of frequency values is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

freq=str2num(get(he,'String'));
if isempty(freq),
    uiwait(errordlg('Number of frequency values must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'NumFreqOld')));
    return;
end

freq=fix(freq);

% Check it for validity

if freq<1,
    uiwait(errordlg('Invalid number of frequency values to use'));
    freq=ss.fit.freqpts;
end

% Update the value and the form

ss.fit.freqpts=freq;
set(he,'String',num2str(freq));
setappdata(hf,'NumFreqOld',freq);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitFreqEdit(varargin)
% Make sure low frequency value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency percentage from the edit box

freq=str2num(get(he,'String'));
if isempty(freq),
    uiwait(errordlg('Frequency range must be a numeric scalar'));
    set(he,'String',num2str(getappdata(hf,'FreqRangeOld')));
    return;
end

% Check it for validity

if freq<=0 || freq > 10
    uiwait(errordlg('Invalid frequency range (0<freq<=10%)'));
    freq=ss.fit.freqp;
end

% Update the value and the form

ss.fit.freqp=freq;
set(he,'String',num2str(freq));
setappdata(hf,'FreqRangeOld',freq);

% Tag that we are not done here

ss.done.autofit=false;


%==========================================================================
function hf=uicbAutoFitPlotFRF(varargin)
% Plot experimental FRF

global ss;

smac_view_frf(ss.fe,'freqrangecc');


%==========================================================================
function hf=uicbAutoFitSave(varargin)
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
function hf=uicbAutoFitBackup(varargin)
% Backup

hf=varargin{1};

if ishandle(hf), delete(hf); end
smac_setprogress('corrcoef');
smac_GUI_plotcorr;


%==========================================================================
function hf=uicbAutoFitQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strmatch(rslt,'Yes'),
    delete(hf);
end
