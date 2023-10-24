function out=smac_GUI_corrcoef(action,varargin)
% SMAC_GUI_CORRCOEF  Calculate correlation coefficients
%
% OUT=SMAC_GUI_CORRCOEF
%
% SMAC_GUI_CORRCOEF allows the user to specify an initial damping estimate,
% number of frequency lines, and bandwidth of curve-fit for calculating the
% correlation coefficients.  The user can select the frequency bandwidth
% either directly or by picking from a mode indicator function (NMIF for
% real normal mode fit and CMIF for complex fit).
%
% OUT will be 1 if the user selects the Calculate Correlation button and the
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
%  14-Jun-2004 / ATA Engineering / Dan Hensley
%    o Change "deltaF" to "Bandwidth"
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Compute NMIF for each reference
%
%  01-Oct-2004 / ATA Engineering / Dan Hensley
%    o Add error checking to numeric inputs
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make GUI non-blocking (don't use uiwait/uiresume), and change
%      forward/backup calls to call the next/previous GUI explicitly
%    o Update for new .mif fields
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

bname='uicbCorrCoef';
funcs.execute='Execute';
funcs.selectpartialfreqrange='SelectPartialFreqRange';
funcs.getfreqrange='GetFreqRange';
funcs.lowfreqedit='LowFreqEdit';
funcs.highfreqedit='HighFreqEdit';
funcs.dampingedit='DampingEdit';
funcs.numlinesedit='NumLinesEdit';

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

uicbCorrCoefCreateFigure(varargin);


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbCorrCoefCreateFigure(varargin)
% Create the Correlation Coefficients form

global ss;

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Correlation Coefficient Calculation', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_CorrCoefForm', ...
    'Visible','off');
hf=handle(hf);

%--------------------------------------------------------------------------
% Add Execute button

hb=ui.CreateButton(hf,'Calculate Correlation',[170 40],uip.Gap*[1 1]);
hb.Callback='smac_GUI_corrcoef(''execute'',gcbf)';
hb.TooltipString='Calculate correlation coefficients for initial root selections';
%uiFitText(hb);

%--------------------------------------------------------------------------
% Add panel for coefficient parameters information

hfr=ui.CreatePanel(hf,' Correlation coefficient parameters',[400 100],'top+gaptop',hb);

if isempty(ss.freqrangecc)
    ss.freqrangecc=ss.freqrange;
end

%----------------------------------
% Add Partial Frequency radiobutton

hr=ui.CreateText(hfr,'Bandwidth of fit',[150 30],[uip.Gap uip.Gap]);
hr.TooltipString='Use the selected frequency range in the correlation coefficient calculation';
hr.HorizontalAlignment='right';

hb=ui.CreateButton(hfr,'Select',[],'right+gapright',hr);
hb.Callback='smac_GUI_corrcoef(''getfreqrange'',gcbf)';
hb.TooltipString='Select frequency range from NMIF or CMIF function';

he=ui.CreateEdit(hfr,num2str(ss.freqrangecc(1,1)),[],'right+gapright',hb);
he.Callback='smac_GUI_corrcoef(''lowfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter starting frequency';
setappdata(hf,'HandleEditFreqLow',he);

ht=ui.CreateText(hfr,'-',[],'right+gapright',he);

he=ui.CreateEdit(hfr,num2str(ss.freqrangecc(2,1)),[],'right+gapright',ht);
he.Callback='smac_GUI_corrcoef(''highfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter ending frequency';
setappdata(hf,'HandleEditFreqHigh',he);

ht=ui.CreateText(hfr,'Hz',[],'right+gapright',he);

% Frequency Lines

hr=ui.CreateText(hfr,'Frequency Lines',[150 30],'top+gaptop',hr);
hr.TooltipString=sprintf(['Use the selected number of frequency lines around each peak\n', ...
    'in the correlation coefficient calculation']);
hr.HorizontalAlignment='right';

he=ui.CreateEdit(hfr,num2str(ss.corr.nl),[],'right+gapright',hr);
he.Callback='smac_GUI_corrcoef(''numlinesedit'',gcbf,gcbo)';
he.TooltipString='Enter the number of spectral lines to use in the calculation';

ht=ui.CreateText(hfr,sprintf('Bandwidth = %g Hz',ss.corr.nl*ss.fe(1).abscissainc),[],'right+gapright',he);
setappdata(hf,'HandleDeltaFText',ht);

% Damping Estimate

hr=ui.CreateText(hfr,'Damping Estimate',[150 30],'top+gaptop',hr);
hr.TooltipString='Use the selected damping estimate in the correlation coefficient calculation';
hr.HorizontalAlignment='right';

he=ui.CreateEdit(hfr,num2str(ss.corr.zeta),[],'right+gapright',hr);
he.Callback='smac_GUI_corrcoef(''dampingedit'',gcbf,gcbo)';
he.TooltipString='Enter damping fraction to use';

% Text gap
hr=ui.CreateText(hfr,' ',[10 uip.Gap],'top+gaptop',hr);

% FIXME:  Make sure the text is sized appropriately

ui.ResizeControl(hfr);

%--------------------------------------------------------------------------
% Add plot FRF button

hb=ui.CreateButton(hf,'Plot FRFs',[100 40],'aligntop+right+gapright+gapright',hfr);
hb.Callback='smac_GUI_corrcoef(''plotfrf'')';
hb.TooltipString='Plot experimental FRF';

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',hb);
hb.Callback='smac_GUI_corrcoef(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';

% Move this button down

hb.Position(2)=uip.Gap;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''nsmac_hlp'')';
hb.TooltipString='Get help on the Correlation Coefficient calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_corrcoef(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';

% Backup button

hb=ui.CreateButton(hf,'Backup',[],'left+gapleft',hb);
hb.Callback='smac_GUI_corrcoef(''backup'',gcbf)';
hb.TooltipString='Return to the pseudo-inverse form';


%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';


%==========================================================================
function hf=uicbCorrCoefExecute(varargin)
% Calculate correlation coefficients

global ss;

hf=watchon;
drawnow;

smac_corr;

smac_setprogress('corrcoef');

watchoff(hf);

% Delete the figure and move to the next one

if ishandle(hf), delete(hf); end
smac_GUI_plotcorr;


%==========================================================================
function hf=uicbCorrCoefGetFreqRange(varargin)
% Select partial frequency range

global ss;

hf=varargin{1};

hh=watchon;
drawnow;

% Get handles

hp(1)=getappdata(hf,'HandleEditFreqLow');
hp(2)=getappdata(hf,'HandleEditFreqHigh');

% Calculate NMIF or CMIF if we need to

if ss.realcomplex==1	% Real normal
    if isempty(ss.mif.exp.nmif),
        ss.mif.exp.nmif=nmif(ss.fe);
        ss.mif.exp.nmif=ss.mif.exp.nmif(1:length(ss.ref_coords));
    end
    mif=ss.mif.exp.nmif;
else
    if isempty(ss.mif.exp.cmif)
        ss.mif.exp.cmif=cmif(ss.fe);
    end
    mif=ss.mif.exp.cmif;
end

% Plot the MIF and ask for the frequency range

h=plot(mif);

[x,~]=ginput(2);
x=sort(x);
delete(h.hf);

watchoff(hh);
drawnow;

% Make sure the X values are within range

fa=mif(1).abscissa;
for k=1:length(x)
    [~,ind]=min(abs(fa-x(k)));
    ss.freqrangecc(k,:)=[fa(ind) ind];
end

% Update the editboxes

for k=1:2
    set(hp(k),'String',num2str(ss.freqrangecc(k,1)));
end

% Tag that we are not done here

ss.done.corr=false;
ss.corr.corr=[];


%==========================================================================
function hf=uicbCorrCoefLowFreqEdit(varargin)
% Make sure low frequency value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

freq=str2num(get(he,'String'));
if numel(freq)~=1,
    uiwait(errordlg('Invalid entry'));
    return;
end

fa=ss.fe(1).abscissa;
[~,ind]=min(abs(fa-freq));

% Make sure this is sane

if fa(ind)>=ss.freqrangecc(2,1)
    uiwait(errordlg('Low frequency value must be lower than high frequency value','Frequency Range Error'));
    set(he,'String',num2str(ss.freqrangecc(1,1)));
else
    ss.freqrangecc(1,:)=[fa(ind) ind];
    set(he,'String',num2str(fa(ind)));
end

% Tag that we are not done here

ss.done.corr=false;
ss.corr.corr=[];


%==========================================================================
function hf=uicbCorrCoefHighFreqEdit(varargin)
% Make sure high frequency value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

freq=str2num(get(he,'String'));
if numel(freq)~=1
    uiwait(errordlg('Invalid entry'));
    return;
end

fa=ss.fe(1).abscissa;
[~,ind]=min(abs(fa-freq));

% Make sure this is sane

if fa(ind)<=ss.freqrangecc(1,1)
    uiwait(errordlg('High frequency value must be higher than low frequency value','Frequency Range Error'));
    set(he,'String',num2str(ss.freqrangecc(2,1)));
else
    ss.freqrangecc(2,:)=[fa(ind) ind];
    set(he,'String',num2str(fa(ind)));
end

% Tag that we are not done here

ss.done.corr=false;
ss.corr.corr=[];


%==========================================================================
function hf=uicbCorrCoefDampingEdit(varargin)
% Make sure damping value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the damping value from the edit box

damp=str2num(get(he,'String'));
if numel(damp)~=1,
    uiwait(errordlg('Invalid entry'));
    return;
end

% Make sure this is sane

if damp<=0 || damp>1
    uiwait(errordlg('Unreasonable damping fraction','Damping Value Error'));
    set(he,'String',num2str(ss.corr.zeta));
else
    ss.corr.zeta=damp;
    set(he,'String',num2str(damp));
end

% Tag that we are not done here

ss.done.corr=false;
ss.corr.corr=[];


%==========================================================================
function hf=uicbCorrCoefNumLinesEdit(varargin)
% Make sure number of spectral lines is sane

global ss;

hf=varargin{1};
he=varargin{2};
ht=getappdata(hf,'HandleDeltaFText');

% Get the number of lines value from the edit box

nl=str2num(get(he,'String'));
if numel(nl)~=1,
    uiwait(errordlg('Invalid entry'));
    return;
end

% Make sure this is sane

if nl<1 || nl>fix(ss.fe(1).numberelements/2)
    uiwait(errordlg(sprintf('Unreasonable number of spectral lines (1 to %d)',fix(ss.fe(1).numberelements/2)), ...
        'Spectral Line Value Error'));
    set(he,'String',num2str(ss.corr.nl));
else
    ss.corr.nl=fix(nl);
    set(he,'String',num2str(nl));
end

% Update the text description

set(ht,'String',sprintf('deltaF = %g Hz',nl*ss.fe(1).abscissainc));

% Tag that we are not done here

ss.done.corr=false;
ss.corr.corr=[];


%==========================================================================
function hf=uicbCorrCoefPlotFRF(varargin)
% Plot experimental FRF

global ss;

smac_view_frf(ss.fe,'freqrangecc');


%==========================================================================
function hf=uicbCorrCoefSave(varargin)
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
function hf=uicbCorrCoefBackup(varargin)
% Backup

hf=varargin{1};
if ishandle(hf), delete(hf); end
smac_GUI_pinv;


%==========================================================================
function hf=uicbCorrCoefQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strmatch(rslt,'Yes'),
    delete(hf);
end
