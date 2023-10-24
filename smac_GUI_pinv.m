function out=smac_GUI_pinv(action,varargin)
% SMAC_GUI_PINV  Calculate SMAC pseudo-inverse
%
% OUT=SMAC_GUI_PINV
%
% SMAC_GUI_PINV allows the user to select a curve-fitting type (real or
% complex), as well as the frequency range over which the pseudo-inverse of
% the FRF matrix will be calculated.  The user can select the frequency
% range either directly or by picking from a mode indicator function (NMIF
% for real normal mode fit and CMIF for complex fit).
%
% OUT will be 1 if the user selects the Execute PINV button and the
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
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Compute pseudo-inverse for each reference
%    o Compute NMIF for each reference
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

bname='uicbPinv';
funcs.execute='Execute';
funcs.selectpartialfreqrange='SelectPartialFreqRange';
funcs.getfreqrange='GetFreqRange';
funcs.lowfreqedit='LowFreqEdit';
funcs.highfreqedit='HighFreqEdit';
funcs.selectfullfreqrange='SelectFullFreqRange';

funcs.selectreal='SelectReal';
funcs.selectcomplex='SelectComplex';
funcs.plotfrf='PlotFRF';

funcs.save='Save';
funcs.quit='Quit';

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

hf=uicbPinvCreateFigure(varargin);


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbPinvCreateFigure(varargin)
% Create the Pseudo-Inverse form

global ss;

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Pseudo Inverse Calculation', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'Tag','SMAC_PinvForm', ...
    'Visible','off');
hf=handle(hf);

%--------------------------------------------------------------------------
% Add Execute button

hb=ui.CreateButton(hf,'Execute Pseudo-Inverse',[170 40],uip.Gap*[1 1]);
hb.Callback='smac_GUI_pinv(''execute'',gcbf)';
hb.TooltipString='Generate pseudo-inverse of FRF matrix';
%uiFitText(hb);

%--------------------------------------------------------------------------
% Add panel for frequency range information

hfr=ui.CreatePanel(hf,' Frequency Range',[400 100],'top+gaptop',hb);

freqrange=ss.fe(1).abscissa([1 end]);

%----------------------------------
% Add Partial Frequency radiobutton

hr=ui.CreateRadio(hfr,'Partial Range',[],[uip.Gap uip.Gap]);
hr.Callback='smac_GUI_pinv(''selectpartialfreqrange'',gcbf)';
hr.TooltipString='Use the selected frequency range in the pseudo-inverse calculation';
hr.Value=0;
setappdata(hf,'HandlePartFreqRange',hr);

hb=ui.CreateButton(hfr,'Select',[],'right+gapright+gapright',hr);
hb.Callback='smac_GUI_pinv(''getfreqrange'',gcbf)';
hb.TooltipString='Select frequency range from NMIF or CMIF function';
hb.Enable='off';
hp=hb;

he=ui.CreateEdit(hfr,num2str(freqrange(1)),[],'right+gapright',hb);
he.Callback='smac_GUI_pinv(''lowfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter starting frequency';
he.Enable='off';
hp(end+1)=he;
setappdata(hf,'HandleEditFreqLow',he);

ht=ui.CreateText(hfr,'-',[],'right+gapright',he);

he=ui.CreateEdit(hfr,num2str(freqrange(2)),[],'right+gapright',ht);
he.Callback='smac_GUI_pinv(''highfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter ending frequency';
he.Enable='off';
hp(end+1)=he;
setappdata(hf,'HandleEditFreqHigh',he);

ht=ui.CreateText(hfr,'Hz',[],'right+gapright',he);
setappdata(hf,'HandlePartialList',hp);

% Add Full Frequency radiobutton

hr=ui.CreateRadio(hfr,sprintf('Full Range (%g-%g Hz) --recommended',freqrange),[],'top+gaptop+gaptop',hr);
hr.Callback='smac_GUI_pinv(''selectfullfreqrange'',gcbf)';
hr.TooltipString='Use the full frequency range in the pseudo-inverse calculation';
hr.Value=1;
setappdata(hf,'HandleFullFreqRange',hr);

ui.ResizeControl(hfr);

%--------------------------------------------------------------------------
% Add panel for solution type

hfr2=ui.CreatePanel(hf,' Solution Method',[400 100],'top+gaptop+gaptop',hfr);

hr=ui.CreateRadio(hfr2,'Complex modes',[],[uip.Gap uip.Gap]);
hr.Callback='smac_GUI_pinv(''selectcomplex'',gcbf)';
hr.TooltipString='Curve-fit complex modes';
hr.Value=0;
setappdata(hf,'HandleComplex',hr);

hr=ui.CreateRadio(hfr2,'Real Normal modes',[],'top+gaptop',hr);
hr.Callback='smac_GUI_pinv(''selectreal'',gcbf)';
hr.TooltipString='Curve-fit real normal modes';
hr.Value=1;
setappdata(hf,'HandleReal',hr);

ui.ResizeControl(hfr2);

% Now make sure the two frames are the same width

pos=[hfr.Position; hfr2.Position];
pos(:,3)=max(pos(:,3));
hfr.Position=pos(1,:);
hfr2.Position=pos(2,:);

%--------------------------------------------------------------------------
% Add FRF and SVD buttons

hb=ui.CreateButton(hf,'Plot FRFs',[100 40],'right+gapright+gapright',hfr2);
hb.Callback='smac_GUI_pinv(''plotfrf'')';
hb.TooltipString='Plot experimental FRF';

hb=ui.CreateButton(hf,'SVD Reduction',[100 40],'bottom+gapbottom+gapbottom',hb);
hb.Callback='uiwait(warndlg(''Not Implemented''));';
hb.TooltipString='Perform SVD reduction on FRFs';
hb.Enable='off';

%--------------------------------------------------------------------------
% Add standard form buttons (save,help,quit)

% Quit button

hb=ui.CreateButton(hf,'QUIT',[],'alignright+bottom',hb);
hb.Callback='smac_GUI_pinv(''quit'',gcbf)';
hb.TooltipString='Exit SMAC';

% Move this button down

pos=hb.Position;
pos(2)=uip.Gap;
hb.Position=pos;

% Help button

hb=ui.CreateButton(hf,'HELP',[],'left+gapleft',hb);
hb.Callback='smac_helper(''pseudo_inv'')';
hb.TooltipString='Get help on the Pseudo-Inverse calculations';

% Save button

hb=ui.CreateButton(hf,'Save',[],'left+gapleft',hb);
hb.Callback='smac_GUI_pinv(''save'',gcbf)';
hb.TooltipString='Save the current data structure to a .mat file';


%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';


%==========================================================================
function hf=uicbPinvExecute(varargin)
% Execute pseudo-inverse

global ss;

% Calculate pseudo-inverse

h=varargin{1};
hf=watchon;
drawnow;

% Compute pseudo-inverse for each reference

for k=1:length(ss.ref_coords),

    fo=ss.fe(:,k).ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:);

    if ss.realcomplex==2,
        ss.pinv(:,:,k) = pinv(fo);
    else
        ss.pinv(:,:,k) = pinv([real(fo); imag(fo)]);
    end
end

smac_setprogress('pinv');

watchoff(hf);

% Delete the figure and move to the next one
if ishandle(hf),
    delete(hf);
end

smac_GUI_corrcoef;


%==========================================================================
function hf=uicbPinvSelectPartialFreqRange(varargin)
% Select partial frequency range

hf=varargin{1};

% Get handles

hr(1)=getappdata(hf,'HandlePartFreqRange');
hr(2)=getappdata(hf,'HandleFullFreqRange');
hp=getappdata(hf,'HandlePartialList');

% Set radio buttons

set(hr(1:2),'Value',0);
set(hr(1),'Value',1);

% Enable partial controls

set(hp,'Enable','on');

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvGetFreqRange(varargin)
% Select partial frequency range

global ss;

hf=varargin{1};

hh=watchon;
drawnow;

% Get handles

hr(1)=getappdata(hf,'HandleReal');
hr(2)=getappdata(hf,'HandleComplex');
hp(1)=getappdata(hf,'HandleEditFreqLow');
hp(2)=getappdata(hf,'HandleEditFreqHigh');

% Calculate NMIF or CMIF if we need to

if hr(1).Value,	% Real normal
    if isempty(ss.mif.exp.nmif),
        ss.mif.exp.nmif=nmif(ss.fe);
        ss.mif.exp.nmif=ss.mif.exp.nmif(1:length(ss.ref_coords));
    end
    mif=ss.mif.exp.nmif;
else
    if isempty(ss.mif.exp.cmif),
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
for k=1:length(x),
    [~,ind]=min(abs(fa-x(k)));
    ss.freqrange(k,:)=[fa(ind) ind];
end

% Update the editboxes

for k=1:2,
    set(hp(k),'String',num2str(ss.freqrange(k,1)));
end

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvLowFreqEdit(varargin)
% Make sure low frequency value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

freq=str2double(get(he,'String'));

fa=ss.fe(1).abscissa;
[freq,ind]=min(abs(fa-freq));

% Make sure this is sane

if fa(ind)>=ss.freqrange(2,1),
    uiwait(errordlg('Low frequency value must be lower than high frequency value', ...
                    'Frequency Range Error'));
    set(he,'String',num2str(ss.freqrange(1,1)));
else
    ss.freqrange(1,:)=[fa(ind) ind];
    set(he,'String',num2str(fa(ind)));
end

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];

%==========================================================================
function hf=uicbPinvHighFreqEdit(varargin)
% Make sure high frequency value is sane

global ss;

hf=varargin{1};
he=varargin{2};

% Get the frequency value from the edit box

freq=str2double(get(he,'String'));

fa=ss.fe(1).abscissa;
[freq,ind]=min(abs(fa-freq));

% Make sure this is sane

if fa(ind)<=ss.freqrange(1,1),
    uiwait(errordlg('High frequency value must be higher than low frequency value', ...
                    'Frequency Range Error'));
    set(he,'String',num2str(ss.freqrange(2,1)));
else
    ss.freqrange(2,:)=[fa(ind) ind];
    set(he,'String',num2str(fa(ind)));
end

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvSelectFullFreqRange(varargin)
% Select full frequency range

global ss;
ss.freqrange=[ss.fe(1).abscissa([1 end]) [1;ss.fe(1).numberelements]];

hf=varargin{1};

% Get handles

hr(1)=getappdata(hf,'HandlePartFreqRange');
hr(2)=getappdata(hf,'HandleFullFreqRange');
hp=getappdata(hf,'HandlePartialList');

% Set radio buttons

set(hr(1:2),'Value',0);
set(hr(2),'Value',1);

% Disable partial controls

set(hp,'Enable','off');

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvSelectReal(varargin)
% Select partial frequency range

global ss;
ss.realcomplex=1;

hf=varargin{1};

% Get handles

hr(1)=getappdata(hf,'HandleReal');
hr(2)=getappdata(hf,'HandleComplex');

% Set radio buttons

set(hr,'Value',0);
set(hr(1),'Value',1);

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvSelectComplex(varargin)
% Select partial frequency range

global ss;
ss.realcomplex=2;

hf=varargin{1};

% Get handles

hr(1)=getappdata(hf,'HandleReal');
hr(2)=getappdata(hf,'HandleComplex');

% Set radio buttons

set(hr,'Value',0);
set(hr(2),'Value',1);

% Tag that we are not done here

ss.done.pinv=false;
ss.pinv=[];


%==========================================================================
function hf=uicbPinvPlotFRF(varargin)
% Plot experimental FRF

global ss;

smac_view_frf(ss.fe,'freqrange');


%==========================================================================
function hf=uicbPinvSave(varargin)
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
function hf=uicbPinvQuit(varargin)
% Quit

hf=varargin{1};

rslt=questdlg('Exit SMAC?');
if strmatch(rslt,'Yes'),
    delete(hf);
end
