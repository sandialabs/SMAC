function out=smac_GUI_residuals(action,varargin)
% SMAC_GUI_RESIDUALS  Residual information setup
%
% OUT=SMAC_GUI_RESIDUALS
%
% SMAC_GUI_RESIDUALS allows the user to set up the various residual options
% to use when synthesizing.
%
% OUT provides the results of the GUI when it is closed.  If the user hits
% Okay, OUT is a residuals structure.  If the user hits Cancel, it is
% empty.

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
%  07-Jun-2005 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Blank out residual values when toggling on the form
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add frequency range to lowmode and highmode
%    o Don't allow 0 Hz frequency
%
%  27-Jul-2005 / ATA Engineering / Dan Hensley
%    o Make sure we check for lowmode and highmode frequency range
%
%  16-Sep-2005 / ATA Engineering / Dan Hensley
%    o Make sure lowmode and highmode frequency range selections are within
%      the frequency bounds of the fit
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
  action='';
end

% Register valid callbacks

bname='uicbResiduals';

funcs.setline='Setline';
funcs.validstartfrequency='ValidStartFrequency';
funcs.validendfrequency='ValidEndFrequency';
funcs.validfrequency='ValidFrequency';
funcs.validdamping='ValidDamping';
funcs.pick_one='PickOne';
funcs.pick_two='PickTwo';

funcs.okay='Okay';
funcs.cancel='Cancel';

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

hf=uicbResidualsCreatePanel(varargin);
notvalid=true;

while notvalid
    uiwait(hf);

    res = getappdata(hf,'ResidualData');

    % Check the data over for validity

    if ~isempty(res)
        out = checkdata(res);
        if ~isempty(out)
            notvalid = false;
        end
    else
        out      = [];
        notvalid = false;
    end
end

% Delete the GUI and return;

if ishandle(hf), delete(hf); end


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbResidualsCreatePanel(varargin)
% Create the Residuals panel

global ss;

hloc=varargin{1}{1};

ui=uiWidgets;
uip=ui.GetParams();

hf = figure(...
    'Color',uip.Form.Color, ...
    'NumberTitle','off', ...
    'MenuBar','none', ...
    'ToolBar','none', ...
    'Resize','on', ...
    'Name','SMAC Residuals', ...
    'BackingStore','on', ...
    'DoubleBuffer','on', ...
    'WindowStyle','modal', ...
    'Tag','SMAC_ResidualsForm', ...
    'CloseRequestFcn','smac_GUI_residuals(''cancel'',gcbf)', ...
    'Visible','off');
hf=handle(hf);

% Determine location on screen for figure

hp=handle(hloc.Parent);
hlocunits=hloc.Units;
hloc.Units='pixels';
hpunits=hp.Units;
hp.Units='pixels';

pos1=hp.Position;
pos2=hloc.Position;

hloc.Units=hlocunits;
hp.Units=hpunits;

hf.Position(1)=pos1(1)+pos2(1)+uip.Gap;
hf.Position(2)=pos1(2)+pos2(2)+uip.Gap;

%--------------------------------------------------------------------------
% Add Okay and Cancel buttons

hb(1)=ui.CreateButton(hf,'Okay',[170 40],uip.Gap*[1 1]);
hb(1).Callback='uiresume(gcbf)';
hb(1).TooltipString='Use current residual parameters';
ui.FitText(hb(1),[0 inf]);

hb(2)=ui.CreateButton(hf,'Cancel',[170 40],'right+gapright',hb(1));
hb(2).Callback='smac_GUI_residuals(''cancel'',gcbf)';
hb(2).TooltipString='Cancel current action';
ui.FitText(hb(2),[0 inf]);

%--------------------------------------------------------------------------
% Create frame that holds the rest of the controls

hfrs=ui.CreatePanel(hf,' Residuals',[500 200],'top+gaptop',hb(1));

%--------------------------------------------------------------------------
% Checkboxes

hc(1)=ui.CreateCheckbox(hfrs,'General Compliance',[170 uip.Editbox.Size(2)],[uip.Gap uip.Gap]);
hc(1).Callback='smac_GUI_residuals(''setline'',gcbo)';
hc(1).TooltipString='Use the specified frequency range to calculate compliance';
hc(1).Value=ss.residuals.compliance.active;
hc(1).Tag='compliance';
ui.FitText(hc(1),[0 inf]);
maxwidth=hc(1).Position(3);

% General Inertance

hc(end+1)=ui.CreateCheckbox(hfrs,'General Inertance',[170 uip.Editbox.Size(2)],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_residuals(''setline'',gcbo)';
hc(end).TooltipString='Use the specified frequency range to calculate inertance';
hc(end).Value=ss.residuals.inertance.active;
hc(end).Tag='inertance';
ui.FitText(hc(end),[0 inf]);
maxwidth=max(maxwidth,hc(end).Position(3));

% % Mode Compliance
% 
% hc(end+1)=ui.CreateCheckbox(hfrs,'Mode Compliance',[170 uip.Editbox.Size(2)],'top+gaptop',hc(end));
% hc(end).Callback='smac_GUI_residuals(''setline'',gcbo)';
% hc(end).TooltipString='Calculate compliance around each root';
% hc(end).Value=ss.residuals.modecompliance.active;
% hc(end).Tag='modeinertance';
% ui.FitText(hc(end),[0 inf]);
% maxwidth=max(maxwidth,hc(end).Position(3));
% 
% % Mode Inertance
% 
% hc(end+1)=ui.CreateCheckbox(hfrs,'Mode Inertance',[170 uip.Editbox.Size(2)],'top+gaptop',hc(end));
% hc(end).Callback='smac_GUI_residuals(''setline'',gcbo)';
% hc(end).TooltipString='Calculate inertance around each root';
% hc(end).Value=ss.residuals.modeinertance.active;
% hc(end).Tag='modecompliance';
% ui.FitText(hc(end),[0 inf]);
% maxwidth=max(maxwidth,hc(end).Position(3));

% High Mode

hc(end+1)=ui.CreateCheckbox(hfrs,'High Mode',[170 uip.Editbox.Size(2)],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_residuals(''setline'',gcbo)';
hc(end).TooltipString='Add an artificial mode at a high frequency';
hc(end).Value=ss.residuals.highmode.active;
hc(end).Tag='highmode';
ui.FitText(hc(end),[0 inf]);
maxwidth=max(maxwidth,hc(end).Position(3));

% Low Mode

hc(end+1)=ui.CreateCheckbox(hfrs,'Low Mode',[170 uip.Editbox.Size(2)],'top+gaptop',hc(end));
hc(end).Callback='smac_GUI_residuals(''setline'',gcbo)';
hc(end).TooltipString='Add an artificial mode at a low frequency';
hc(end).Value=ss.residuals.lowmode.active;
hc(end).Tag='lowmode';
ui.FitText(hc(end),[0 inf]);
maxwidth=max(maxwidth,hc(end).Position(3));

% Resize the checkboxes based on the maximum width

for k=1:length(hc);
    hc(k).Position(3)=maxwidth;
end

%--------------------------------------------------------------------------
% Add the options to each

% Compliance

clear ht he hb
ht(1)=ui.CreateText(hfrs,' -- Frequency (Hz)',[],'right+gapright',hc(1));
ui.FitText(ht(1),[0 inf]);

he(1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(1));
he(1).Callback='smac_GUI_residuals(''validstartfrequency'',gcbo,''Compliance'')';
he(1).String='';
he(1).TooltipString='Start of frequency range';
if ~isempty(ss.residuals.compliance.freq)
    he(1).String=num2str(ss.residuals.compliance.freq(1));
end

ht(2)=ui.CreateText(hfrs,' to ',[],'right',he(1));
ui.FitText(ht(2),[0 inf]);

he(2)=ui.CreateEdit(hfrs,'',[],'right',ht(2));
he(2).Callback='smac_GUI_residuals(''validendfrequency'',gcbo,''Compliance'')';
he(2).String='';
he(2).TooltipString='End of frequency range';
if ~isempty(ss.residuals.compliance.freq)
    he(2).String=num2str(ss.residuals.compliance.freq(2));
end

hb=ui.CreateButton(hfrs,'Pick Range',[],'right+gapright',he(2));
hb.Callback='smac_GUI_residuals(''pick_two'',gcbf,''Compliance'')';
hb.TooltipString='Pick frequency range from a plot';
ui.FitText(hb,[0 inf]);

setappdata(hc(1),'HandleRest',[he ht hb]);
setappdata(hf,'HandleComplianceData',he);
uicbResidualsSetline(hc(1));

%--------------------------------------------------------------------------
% Inertance

clear ht he hb
ht(1)=ui.CreateText(hfrs,' -- Frequency (Hz)',[],'right+gapright',hc(2));
ui.FitText(ht(1),[0 inf]);

he(1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(1));
he(1).Callback='smac_GUI_residuals(''validstartfrequency'',gcbo,''Inertance'')';
he(1).String='';
he(1).TooltipString='Start of frequency range';
if ~isempty(ss.residuals.inertance.freq)
    he(1).String=num2str(ss.residuals.inertance.freq(1));
end

ht(2)=ui.CreateText(hfrs,' to ',[],'right',he(1));
ui.FitText(ht(2),[0 inf]);

he(2)=ui.CreateEdit(hfrs,'',[],'right',ht(2));
he(2).Callback='smac_GUI_residuals(''validendfrequency'',gcbo,''Inertance'')';
he(2).String='';
he(2).TooltipString='End of frequency range';
if ~isempty(ss.residuals.inertance.freq)
    he(2).String=num2str(ss.residuals.inertance.freq(2));
end

hb=ui.CreateButton(hfrs,'Pick Range',[],'right+gapright',he(2));
hb.Callback='smac_GUI_residuals(''pick_two'',gcbf,''Inertance'')';
hb.TooltipString='Pick frequency range from a plot';
ui.FitText(hb,[0 inf]);

setappdata(hc(2),'HandleRest',[he ht hb]);
setappdata(hf,'HandleInertanceData',he);
uicbResidualsSetline(hc(2));

%--------------------------------------------------------------------------
% High Mode

clear ht he hb
ht(1)=ui.CreateText(hfrs,' -- Frequency (Hz)',[],'right+gapright',hc(end-1));
ui.FitText(ht(1),[0 inf]);

he(1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(1));
he(1).Callback='smac_GUI_residuals(''validfrequency'',gcbo,''Highmode'')';
he(1).TooltipString='Frequency for the high mode';
if isempty(ss.residuals.highmode.freq)
    ss.residuals.highmode.freq=ss.fe(1).abscissa(end);
end
he(1).String=num2str(ss.residuals.highmode.freq);

ht(2)=ui.CreateText(hfrs,', Damping (%)',[],'right',he(1));
ui.FitText(ht(2),[0 inf]);

he(2)=ui.CreateEdit(hfrs,'',[],'right',ht(2));
he(2).Callback='smac_GUI_residuals(''validdamping'',gcbo,''Highmode'')';
he(2).TooltipString='Damping for this mode';
if isempty(ss.residuals.highmode.damp)
    ss.residuals.highmode.damp=1e-8;
end
he(2).String=num2str(ss.residuals.highmode.damp*100);

hb=ui.CreateButton(hfrs,'Pick Freq',[],'right+gapright',he(2));
hb.Callback='smac_GUI_residuals(''pick_one'',gcbf,''Highmode'')';
hb.TooltipString='Pick frequency from a plot';
ui.FitText(hb,[0 inf]);

%

ht(end+1)=ui.CreateText(hfrs,' -- Frequency Range (Hz)',[],'right+gapright',hb);
ui.FitText(ht(end),[0 inf]);

he(end+1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(end));
he(end).Callback='smac_GUI_residuals(''validstartfrequency'',gcbo,''Highmode'')';
he(end).String='';
he(end).TooltipString='Start of frequency range';
if ~isempty(ss.residuals.highmode.frange)
    he(end).String=num2str(ss.residuals.highmode.frange(1));
end

ht(end+1)=ui.CreateText(hfrs,' to ',[],'right',he(end));
ui.FitText(ht(end),[0 inf]);

he(end+1)=ui.CreateEdit(hfrs,'',[],'right',ht(end));
he(end).Callback='smac_GUI_residuals(''validendfrequency'',gcbo,''Highmode'')';
he(end).String='';
he(end).TooltipString='End of frequency range';
if ~isempty(ss.residuals.highmode.frange)
    he(end).String=num2str(ss.residuals.highmode.frange(2));
end

hb(end+1)=ui.CreateButton(hfrs,'Pick Range',[],'right+gapright',he(end));
hb(end).Callback='smac_GUI_residuals(''pick_two'',gcbf,''Highmode'')';
hb(end).TooltipString='Pick frequency range from a plot';
ui.FitText(hb(end),[0 inf]);

%

setappdata(hc(end-1),'HandleRest',[he ht hb]);
setappdata(hf,'HandleHighmodeData',he(1:2));
setappdata(hf,'HandleHighmodeRangeData',he(3:4));
uicbResidualsSetline(hc(end-1));

%--------------------------------------------------------------------------
% Low Mode

clear ht he hb
ht(1)=ui.CreateText(hfrs,' -- Frequency (Hz)',[],'right+gapright',hc(end));
ui.FitText(ht(1),[0 inf]);

he(1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(1));
he(1).Callback='smac_GUI_residuals(''validfrequency'',gcbo,''Lowmode'')';
he(1).TooltipString='Frequency for the low mode';
if isempty(ss.residuals.lowmode.freq)
    ss.residuals.lowmode.freq=ss.fe(1).abscissa(1);
    if ss.residuals.lowmode.freq==0
        ss.residuals.lowmode.freq=ss.fe(1).abscissa(2);
    end
end
he(1).String=num2str(ss.residuals.lowmode.freq);

ht(2)=ui.CreateText(hfrs,', Damping (%)',[],'right',he(1));
ui.FitText(ht(2),[0 inf]);

he(2)=ui.CreateEdit(hfrs,'',[],'right',ht(2));
he(2).Callback='smac_GUI_residuals(''validdamping'',gcbo,''Lowmode'')';
he(2).TooltipString='Damping for this mode';
if isempty(ss.residuals.lowmode.damp)
    ss.residuals.lowmode.damp=1e-8;
end
he(2).String=num2str(ss.residuals.lowmode.damp*100);

hb=ui.CreateButton(hfrs,'Pick Freq',[],'right+gapright',he(2));
hb.Callback='smac_GUI_residuals(''pick_one'',gcbf,''Lowmode'')';
hb.TooltipString='Pick frequency from a plot';
ui.FitText(hb,[0 inf]);

%

ht(end+1)=ui.CreateText(hfrs,' -- Frequency Range (Hz)',[],'right+gapright',hb);
ui.FitText(ht(end),[0 inf]);

he(end+1)=ui.CreateEdit(hfrs,'',[],'right+gapright',ht(end));
he(end).Callback='smac_GUI_residuals(''validstartfrequency'',gcbo,''Lowmode'')';
he(end).String='';
he(end).TooltipString='Start of frequency range';
if ~isempty(ss.residuals.lowmode.frange)
    he(end).String=num2str(ss.residuals.lowmode.frange(1));
end

ht(end+1)=ui.CreateText(hfrs,' to ',[],'right',he(end));
ui.FitText(ht(end),[0 inf]);

he(end+1)=ui.CreateEdit(hfrs,'',[],'right',ht(end));
he(end).Callback='smac_GUI_residuals(''validendfrequency'',gcbo,''Lowmode'')';
he(end).String='';
he(end).TooltipString='End of frequency range';
if ~isempty(ss.residuals.lowmode.frange)
    he(end).String=num2str(ss.residuals.lowmode.frange(2));
end

hb(end+1)=ui.CreateButton(hfrs,'Pick Range',[],'right+gapright',he(end));
hb(end).Callback='smac_GUI_residuals(''pick_two'',gcbf,''Lowmode'')';
hb(end).TooltipString='Pick frequency range from a plot';
ui.FitText(hb(end),[0 inf]);

%

setappdata(hc(end),'HandleRest',[he ht hb]);
setappdata(hf,'HandleLowmodeData',he(1:2));
setappdata(hf,'HandleLowmodeRangeData',he(3:4));
uicbResidualsSetline(hc(end));

%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeControl(hfrs);
ui.ResizeFigure(hf);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';
setappdata(hf,'ResidualData',ss.residuals);

return;


%==========================================================================
%==========================================================================
function hf=uicbResidualsValidFrequency(varargin)
% Make sure value is ok

hc=varargin{1};
type=varargin{2};

hf = hc;
while ~strcmpi(get(hf,'Type'),'figure'), hf = get(hf,'Parent'); end
res=getappdata(hf,'ResidualData');

hd=getappdata(hf,sprintf('Handle%sData',type));

val=str2num(get(hc,'String'));
if isempty(val),
    if ~isempty(get(hc,'String'))
        uiwait(errordlg('Invalid number'));
    end
    set(hc,'String','');
    return;
end

% Make sure the value is valid

if  length(val)~=1 || val<=0
    uiwait(errordlg('Frequency must a be positive scalar'));
    set(hc,'String','');
    return;
end
   
% Store the frequency if everything is ok to this point    
    
res.(lower(type)).freq=val;
set(hc,'String',num2str(val));

% Store the data

setappdata(hf,'ResidualData',res);


%==========================================================================
function hf = uicbResidualsValidDamping(varargin)
% Make sure value is ok

hc   = varargin{1};
type = varargin{2};

hf = hc;
while ~strcmpi(get(hf,'Type'),'figure'), hf = get(hf,'Parent'); end
res = getappdata(hf,'ResidualData');

hd = getappdata(hf,sprintf('Handle%sData',type));

val = str2num(get(hc,'String'));
if isempty(val)
    if ~isempty(get(hc,'String'))
        uiwait(errordlg('Invalid number'));
    end
    set(hc,'String','');
    return;
end

% Make sure the value is valid

if numel(val)~=1 || val<=0 || val>50
    uiwait(errordlg('Damping must a be positive scalar between 0 and 50'));
    set(hc,'String','');
    return;
end
   
% Store the frequency if everything is ok to this point    
    
res.(lower(type)).damp = val/100;
set(hc,'String',num2str(val));

% Store the data

setappdata(hf,'ResidualData',res);


%==========================================================================
function hf=uicbResidualsValidStartFrequency(varargin)
% Make sure value is ok

global ss

hc   = varargin{1};
type = varargin{2};

% Get frequency range field name based on type

switch type
    case {'Lowmode','Highmode'}
        freqfield = 'frange';
        modif     = 'Range';
    otherwise
        freqfield = 'freq';
        modif     = '';
end

%

hf = hc;
while ~strcmpi(get(hf,'Type'),'figure'), hf = get(hf,'Parent'); end
res=getappdata(hf,'ResidualData');

hd = getappdata(hf,sprintf('Handle%sData',type));

val = str2num(get(hc,'String'));
if isempty(val)
    if ~isempty(get(hc,'String'))
        uiwait(errordlg('Invalid number'));
    end
    set(hc,'String','');
    return;
end

if numel(val)~=1
    uiwait(errordlg('Frequency must a be scalar'));
    set(hc,'String','');
    return;
end

% Make sure the value is valid

fa = ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2));
[~,ind]=min(abs(fa-val));

val=fa(ind);
if val==0
    val=fa(ind+1);
end
res.(lower(type)).(freqfield)(1)=val;
set(hc,'String',num2str(val));

if ~isempty(res.(lower(type)).(freqfield))
    if length(res.(lower(type)).(freqfield))>1 && val>=res.(lower(type)).(freqfield)(2)
        uiwait(errordlg('Start frequency must be less than end frequency'));
        set(hc,'String','');
        res.(lower(type)).(freqfield)(1)=inf;
    end
end

% Store the data

setappdata(hf,'ResidualData',res);


%==========================================================================
function hf=uicbResidualsValidEndFrequency(varargin)
% Make sure value is ok

global ss

hc=varargin{1};
type=varargin{2};

% Get frequency range field name based on type

switch type
    case {'Lowmode','Highmode'}
        freqfield='frange';
        modif='Range';
    otherwise
        freqfield='freq';
        modif='';
end

%

hf = hc;
while ~strcmpi(get(hf,'Type'),'figure'), hf = get(hf,'Parent'); end
res = getappdata(hf,'ResidualData');

hd  = getappdata(hf,sprintf('Handle%s%sData',type,modif));

val = str2num(get(hc,'String'));
if isempty(val)
    if ~isempty(get(hc,'String'))
        uiwait(errordlg('Invalid number'));
    end
    set(hc,'String','');
    return;
end

if numel(val)~=1
    uiwait(errordlg('Frequency must a be scalar'));
    set(hc,'String','');
    return;
end

% Make sure the value is valid

fa = ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2));
[~,ind]=min(abs(fa-val));

val=fa(ind);
if val==0
    val=fa(ind+1);
end
if isempty(res.(lower(type)).(freqfield))
    res.(lower(type)).(freqfield)=[inf val];
else
    res.(lower(type)).(freqfield)(2)=val;
end
set(hc,'String',num2str(val));

if ~isempty(res.(lower(type)).(freqfield))
    if val<=res.(lower(type)).(freqfield)(1)
        uiwait(errordlg('End frequency must be greater than start frequency'));
        set(hc,'String','');
        res.(lower(type)).(freqfield)(2)=inf;
    end
end

% Store the data

setappdata(hf,'ResidualData',res);


%==========================================================================
function hf=uicbResidualsPickOne(varargin)
% Pick a single frequency from a plot

hf   = varargin{1};
type = varargin{2};
res  = getappdata(hf,'ResidualData');

freq = guipick(1);

res.(lower(type)).freq = freq;

% Store the data

setappdata(hf,'ResidualData',res);

% Set the fields

hd = getappdata(hf,sprintf('Handle%sData',type));
set(hd(1),'String',num2str(freq(1)));


%==========================================================================
function hf=uicbResidualsPickTwo(varargin)
% Pick frequency range from a plot

global ss

hf   = varargin{1};
type = varargin{2};
res = getappdata(hf,'ResidualData');

% Get frequency range field name based on type

switch type
    case {'Lowmode','Highmode'}
        freqfield = 'frange';
        modif     = 'Range';
    otherwise
        freqfield = 'freq';
        modif     = '';
end

% Pick the frequencies

freq = guipick(2,ss.freqrange(:,1));

% Store the data

res.(lower(type)).(freqfield) = freq(:).';
setappdata(hf,'ResidualData',res);

% Set the fields

hd=getappdata(hf,sprintf('Handle%s%sData',type,modif));
set(hd(1),'String',num2str(freq(1)));
set(hd(2),'String',num2str(freq(2)));


%==========================================================================
%==========================================================================
function hf = uicbResidualsSetline(varargin)
% Enable or disable the current line

hc  = varargin{1};
val = double(get(hc,'Value'))+1;

enabled = {'off','on'};
setline(hc,'Enable',enabled{val});

% Set the residual data as well

hf = hc;
while ~strcmpi(get(hf,'Type'),'figure'), hf = get(hf,'Parent'); end

res = getappdata(hf,'ResidualData');
res.(get(hc,'Tag')).active   = get(hc,'Value');
res.(get(hc,'Tag')).residual = [];
setappdata(hf,'ResidualData',res);


%==========================================================================
function hf=uicbResidualsCancel(varargin)
% User hit Cancel

% Continue

hf = varargin{1};

setappdata(hf,'ResidualData',[]);
uiresume(hf);


%==========================================================================
%==========================================================================
function setline(hc,type,val)
% Set parameters related to the supplied handle

hr = getappdata(hc,'HandleRest');
if isempty(hr), return; end

set(hr,type,val);

%==========================================================================
function freq=guipick(num,frange)
% Select frequencies from a plot

global ss

% Plot the drive point FRF

refcoord=abs(allref(ss.fe(:)));
rescoord=abs(allres(ss.fe(:)));

ind=false(size(refcoord));
for k=1:length(ss.ref_coords)
    ind=ind | (refcoord == abs(ss.ref_coords(k)) & rescoord == abs(ss.ref_coords(k)));
end
if ~any(ind), ind(1)=true; end

hp=plot(ss.fe(ind),'xwin',ss.freqrange(:,1));

[x,~]=ginput(num);

freq=sort(x);

fa=ss.fe(1).abscissa;
if exist('frange','var')
    fa=fa(fa>=frange(1) & fa<=frange(2));
end

for k=1:length(freq)
    [~,ind]=min(abs(fa-freq(k)));
    freq(k)=fa(ind);
    if freq(k)==0
        freq(k)=fa(ind+1);
    end
end

% Delete the figure

delete(hp.hf(ishandle(hp.hf)));

%==========================================================================
function out=checkdata(res)
% Check over residual structure

out=[];

% Check for frequency range on general inertance and compliance

check={'General Inertance',  'inertance',  'freq'; ...
	   'General Compliance', 'compliance', 'freq'; ...
	   'Low Mode',           'lowmode',    'frange'; ...
	   'High Mode',          'highmode',   'frange'};

for k=1:size(check,1)
    if res.(check{k,2}).active
        if length(res.(check{k,2}).(check{k,3}))~=2 || any(isinf(res.(check{k,2}).(check{k,3})))
            uiwait(errordlg(sprintf('Must enter both frequencies for %s',check{k,1})));
            return;
        end
    end
end


out=res;
