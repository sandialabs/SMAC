function [out,cormax,refmax]=smac_GUI_addroot(action,varargin)
% SMAC_GUI_ADDROOT  Manually add a root
%
% [OUT,CORMAX,REFMAX]=SMAC_GUI_ADDROOT
%
% SMAC_GUI_ADDROOT provides a way for the user to manually add a root by
% performing frequency and damping fits around user-specified frequency and
% damping ranges.  It is an iterative process where the user narrows in on
% frequency and damping ranges until a parabola curve-fit of the resulting
% correlation coefficient tightly overlays.
%
% OUT is a 1x2 vector containing the frequency, damping fraction, and
% average correlation coefficient (between the frequency and damping fits) 
% for the root fit finally selected.  If the user cancels the form, OUT
% will be empty.  CORMAX is a vector of maximum correlation coefficients per
% reference.  REFMAX is the reference giving the maximum correlation
% coefficient.

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
%  07-Jun-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  08-Jun-2004 / ATA Engineering / Dan Hensley
%    o Add Help button
%    o Pass out correlation coefficient as well
%    o Don't allow invalid frequency and damping entries in the edit boxes
%
%  15-Jun-2004 / ATA Engineering / Dan Hensley
%    o Fix Best Fit placement to better work on different resolutions
%    o Show initial best fit frequency and damping as 0
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Only update one fit at a time, not both
%    o Add Calculate buttons to refit frequency or damping
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Move damping fit buttons to just above damping plot
%    o Suppress polynomial warnings
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for changed autofit algorithms (loop over references)
%    o Add output argument cormax
%
%  07-Jun-2005 / ATA Engineering / Dan Hensley
%    o Use new GUI functions in uiWidgets
%    o Fix bug when canceling Add root (output not assigned)
%
%  10-Sep-2005 / ATA Engineering / Dan Hensley
%    o Change maximum damping value to 25%
%
%  22-Nov-2013 / ATA Engineering / Bill Fladung
%    o Commented out one line that wasn't being used, but was causing an error.
%
%==========================================================================

% Check input arguments

if ~exist('action','var'),
  action='';
end

% Register valid callbacks

bname='uicbAddRoot';

funcs.addroot='AddRoot';

funcs.dampingfit='DampingFit';
funcs.frequencyfit='FrequencyFit';

funcs.getdamprange='GetRange';
funcs.getfreqrange='GetRange';

funcs.calculate='Calculate';

funcs.zoomoutdamp='ZoomOutDamp';
funcs.zoomoutfreq='ZoomOutFreq';

funcs.lowdampedit='LowDampEdit';
funcs.highdampedit='HighDampEdit';
funcs.lowfreqedit='LowFreqEdit';
funcs.highfreqedit='HighFreqEdit';

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

hf=uicbAddRootCreateFigure(varargin);

% Wait for the user to finish

uiwait(hf);

% Check for backup as well

if ~ishandle(hf),
  out=[];
  cormax=[];
  refmax=[];
else

  % Extract the best-fit frequency and damping and pass those out

  hfreq=getappdata(hf,'HandleBestFreq');
  hdamp=getappdata(hf,'HandleBestDamp');
  cc=[getappdata(hdamp,'CorrCoef') getappdata(hfreq,'CorrCoef')];

  out=[str2num(get(hfreq,'String')) str2num(get(hdamp,'String'))/100 mean(cc)];
  cormax=getappdata(hf,'CorMax');
  refmax=getappdata(hf,'RefMax');

  % Delete the figure and exit

  delete(hf);
end


%==========================================================================
% CALLBACKS
%==========================================================================
function hf=uicbAddRootCreateFigure(varargin)
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
	 'Name','SMAC Manual Add Root', ...
	 'BackingStore','on', ...
	 'DoubleBuffer','on', ...
	 'Tag','SMAC_AddRootForm', ...
	 'Visible','off');
hf=handle(hf);
setappdata(hf,'Backup',false);

%--------------------------------------------------------------------------
% Add Add Root button

hb=ui.CreateButton(hf,'Add Root',[170 40],uip.Gap*[1 1]+[50 0]);
hb.Callback='uiresume(gcbf)';
hb.TooltipString='Add the current root to the list';
%uiFitText(hb);
hb1=hb;

% Add damping axis

scr=get(0,'ScreenSize');
scr=scr(3:4);
asize=[800 300].*scr./[1600 1200];	% Resize based on 1600x1200 resolution

ha(1)=ui.CreateAxis(hf,'Damping',asize,'top+gaptop',hb);
ha(1).Position(2)=ha.Position(2)+70;
setappdata(hf,'HandleAxisDamping',ha(1));

% Add Cancel button

hb=ui.CreateButton(hf,'Cancel',[170 40],'right+gapright',hb);
hb.Callback='delete(gcbf)';
hb.TooltipString='Cancel without adding the current root';

% Add Damping Fit text

hr(1)=ui.CreateText(hf,'Damping Fit',[],'top+gaptop+gaptop',ha(1));
hr(1).TooltipString='Perform a fit of the damping value over the specified damping range';

% Add frequency axis

ha(2)=ui.CreateAxis(hf,'Frequency',asize,'top+gaptop',hr(1));
ha(2).Position(2)=ha(2).Position(2)+50;
setappdata(hf,'HandleAxisFrequency',ha(2));

%--------------------------------------------------------------------------
% Add Damping Fit and Frequency Fit rows

hr(2)=ui.CreateText(hf,'Frequency Fit',[],'top+gaptop+gaptop',ha(2));
hr(2).TooltipString='Perform a fit of the frequency value over the specified frequency range';

% Resize the width of both to the maximum

maxwidth=max(hr(1).Position(3),hr(2).Position(3));
hr(1).Position(3)=maxwidth;
hr(2).Position(3)=maxwidth;

%-------------------------------------------
% Add the rest of the damping controls

hb=ui.CreateButton(hf,'Select',[],'right+gapright',hr(1));
hb.Callback='smac_GUI_addroot(''getdamprange'',gcbf,''Damp'')';
hb.TooltipString='Select damping range from the plot';

he=ui.CreateEdit(hf,num2str(ss.fit.zeta(1)*100),[],'right+gapright',hb);
he.Callback='smac_GUI_addroot(''lowdampedit'',gcbf,gcbo)';
he.TooltipString='Enter starting damping percentage';
setappdata(hf,'HandleEditDampLow',he);

ht=ui.CreateText(hf,'-',[],'right+gapright',he);

he=ui.CreateEdit(hf,num2str(ss.fit.zeta(2)*100),[],'right+gapright',ht);
he.Callback='smac_GUI_addroot(''highdampedit'',gcbf,gcbo)';
he.TooltipString='Enter ending damping percentage';
setappdata(hf,'HandleEditDampHigh',he);

ht=ui.CreateText(hf,'%',[],'right+gapright',he);

hbz(1)=ui.CreateButton(hf,'Zoom Out',[],'right+gapright',ht);
hbz(1).Callback='smac_GUI_addroot(''zoomoutdamp'',gcbf)';
hbz(1).TooltipString='Double the current damping range';
ui.FitText(hbz(1),[0 inf]);

hbc(1)=ui.CreateButton(hf,'Calculate',[],'right+gapright',hbz(1));
hbc(1).Callback='smac_GUI_addroot(''calculate'',''Damp'',gcbf)';
hbc(1).TooltipString='Recalculate the damping using the current range';
hbc(1).Enable='off';
ui.FitText(hbc(1),[0 inf]);

% Remember the last entered values

setappdata(hf,'DampingOld',100*ss.fit.zeta(:).');

%-------------------------------------------
% Add the rest of the frequency controls

hb=ui.CreateButton(hf,'Select',[],'right+gapright',hr(2));
hb.Callback='smac_GUI_addroot(''getfreqrange'',gcbf,''Freq'')';
hb.TooltipString='Select frequency range from the plot';

he=ui.CreateEdit(hf,num2str(ss.freqrangecc(1,1)),[],'right+gapright',hb);
he.Callback='smac_GUI_addroot(''lowfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter starting frequency';
setappdata(hf,'HandleEditFreqLow',he);

ht=ui.CreateText(hf,'-',[],'right+gapright',he);

he=ui.CreateEdit(hf,num2str(ss.freqrangecc(2,1)),[],'right+gapright',ht);
he.Callback='smac_GUI_addroot(''highfreqedit'',gcbf,gcbo)';
he.TooltipString='Enter ending frequency';
setappdata(hf,'HandleEditFreqHigh',he);

ht=ui.CreateText(hf,'Hz',[],'right+gapright',he);

hbz(2)=ui.CreateButton(hf,'Zoom Out',[],'right+gapright',ht);
hbz(2).Callback='smac_GUI_addroot(''zoomoutfreq'',gcbf)';
hbz(2).TooltipString='Double the current frequency range';
ui.FitText(hbz(2),[0 inf]);

hbc(2)=ui.CreateButton(hf,'Calculate',[],'right+gapright',hbz(2));
hbc(2).Callback='smac_GUI_addroot(''calculate'',''Freq'',gcbf)';
hbc(2).TooltipString='Recalculate the frequency using the current range';
hbc(2).Enable='off';
ui.FitText(hbc(2),[0 inf]);
setappdata(hf,'HandleCalculate',hbc);

% Remember the last entered values

setappdata(hf,'FrequencyOld',ss.freqrangecc(:,1).');

% Move the two zoom buttons so they are aligned

posmax=max(hbz(1).Position(1),hbz(2).Position(1));
hbz(1).Position(1)=posmax;
hbz(2).Position(1)=posmax;

% Update the Calculate button positions

hbc(1).Position(1)=posmax+hbz(1).Position(3)+uip.Gap;
hbc(2).Position(1)=posmax+hbz(1).Position(3)+uip.Gap;

%--------------------------------------------------------------------------
% Add Best Fit panel

hfr=ui.CreatePanel(hf,' Best Fit',[400 100],'right+gapright+gapright+gapright+gapright',hbc(2));

he=ui.CreateEdit(hfr,'0',[],[uip.Gap uip.Gap]);
he.Enable='inactive';
he.BackgroundColor='yellow';
setappdata(hf,'HandleBestDamp',he);
setappdata(he,'CorrCoef',0);

ht=ui.CreateText(hfr,'%',[],'right+gapright',he);
% hl(end+1)=ht;  WAF131122: This is not used and was causing Randy Mayes an error.

he=ui.CreateEdit(hfr,'0',[],'top+gaptop',he);
he.Enable='inactive';
he.BackgroundColor='yellow';
setappdata(hf,'HandleBestFreq',he);
setappdata(he,'CorrCoef',0);

ht=ui.CreateText(hfr,'Hz',[],'right+gapright',he);

ui.ResizeControl(hfr);

% Shrink the axes to the Best Fit width

ha(1).Position(3)=hfr.Position(1)+hfr.Position(3)-ha(1).Position(1);
ha(2).Position(3)=hfr.Position(1)+hfr.Position(3)-ha(2).Position(1);

%--------------------------------------------------------------------------
% Store the FRF data with the figure (for performance reasons)

wr=2*pi*ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2)).';
frf=ss.fe.ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:);
frf=reshape(frf,[size(frf,1) size(ss.fe)]);

setappdata(hf,'Abscissa',wr);
setappdata(hf,'Ordinate',frf);

% Fill in the fit data with correlation coefficient information

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
setappdata(hf,'FreqFitData',[xf(:) ss.corr.corr(:,1) zeros(size(ss.corr.corr(:,1)))]);
uicbAddRootUpdatePlot(hf,'Freq');

setappdata(hf,'DampFitData',[ss.fit.zeta(:) zeros(length(ss.fit.zeta(:)),2)]);
uicbAddRootUpdatePlot(hf,'Damp');

%--------------------------------------------------------------------------
% Make sure the figure encompasses all of the controls

ui.ResizeFigure(hf);

%--------------------------------------------------------------------------
% Add a Help button

hb=ui.CreateButton(hf,'HELP',[],'right',hb1);
hb.Callback='smac_helper(''mfit_hlp'')';
hb.TooltipString='Get help on the Correlation Coefficient calculations';
hb.Position(1)=ha(1).Position(1)+ha(1).Position(3)-hb.Position(3);

% Make all uicontrols normalized

ui.SetUnits(hf,'normalized');

hf.Visible='on';


%==========================================================================
function hf=uicbAddRootGetRange(varargin)
% Get the damping range and perform a damping fit

global ss;

hf=varargin{1};
type=varargin{2};

switch type
  case 'Freq'
    t{1}='Frequency';
    t{2}='Damping';

  case 'Damp'
    t{1}='Damping';
    t{2}='Frequency';

end

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');

% Make sure only the plot we want can be selected

ha=getappdata(hf,['HandleAxis' t{1}]);
set(ha,'HitTest','on','HandleVisibility','on');
set(hf,'CurrentAxes',ha);

ha=getappdata(hf,['HandleAxis' t{2}]);
set(ha,'HitTest','off','HandleVisibility','off');

% Prompt for the range to select

[x,~]=ginput(2);
x=sort(x);

% Set the edit boxes appropriately

he(1)=getappdata(hf,['HandleEdit' type 'Low']);
he(2)=getappdata(hf,['HandleEdit' type 'High']);
set(he(1),'String',num2str(x(1)));
set(he(2),'String',num2str(x(2)));

% Perform the fits

switch type
  case 'Freq',
    uicbAddRootFreqFit(hf);

  case 'Damp',
    uicbAddRootDampFit(hf);

end


%==========================================================================
function hf=uicbAddRootZoomOutDamp(varargin)
% Zoom out on the damping axis

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditDampLow');
he(2)=getappdata(hf,'HandleEditDampHigh');

% Bump the damping range by half the difference

zeta=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];
dz=diff(zeta)/2;
zeta=zeta+dz*[-1 1];

% Make sure we have sane values

zeta(zeta<=0)=.001;
zeta(zeta>25)=25;

% Put these values back into the form

set(he(1),'String',num2str(zeta(1)));
set(he(2),'String',num2str(zeta(2)));

% Perform the fits

%uicbAddRootFreqFit(hf);
uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootZoomOutFreq(varargin)
% Zoom out on the damping axis

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditFreqLow');
he(2)=getappdata(hf,'HandleEditFreqHigh');

% Bump the frequency range by half the difference

freq=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];
df=diff(freq)/2;
freq=freq+df*[-1 1];

% Find the closest spectral line to the selected range
% FIXME:  No, don't do it this way

freq(freq<ss.freqrangecc(1,1))=ss.freqrangecc(1,1);
freq(freq>ss.freqrangecc(2,1))=ss.freqrangecc(2,1);

% Put these values back into the form

set(he(1),'String',num2str(freq(1)));
set(he(2),'String',num2str(freq(2)));

% Perform the fits

uicbAddRootFreqFit(hf);
%uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootLowDampEdit(varargin)
% Edit lower bound for damping

% FIXME:  Get valid damping range from the structure

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditDampLow');
he(2)=getappdata(hf,'HandleEditDampHigh');
dampold=getappdata(hf,'DampingOld');

% Make sure the value is a valid number

if isempty(str2num(get(he(1),'String'))),
  uiwait(errordlg('Invalid entry'));
  set(he(1),'String',num2str(dampold(1)));
  return;
end

% Get the damping values from the edit box

damp=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];

% Make sure the damping value makes sense

if damp(1)<=0 || damp(1)>25,
  uiwait(errordlg('Unreasonable damping','Damping Value Error'));
  damp(1)=dampold(1);
elseif damp(1)>damp(2),
  uiwait(errordlg('Low value for damping must be lower than high value','Damping Value Error'));
  damp(1)=dampold(1);
end
damp(damp<0)=0.001;
damp(damp>25)=25;

% Update the form

set(he(1),'String',num2str(damp(1)));
setappdata(hf,'DampingOld',damp);

% Perform the fits

%uicbAddRootFreqFit(hf);
uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootHighDampEdit(varargin)
% Edit higher bound for damping

% FIXME:  Get valid damping range from the structure

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditDampLow');
he(2)=getappdata(hf,'HandleEditDampHigh');
dampold=getappdata(hf,'DampingOld');

% Make sure the value is a valid number

if isempty(str2num(get(he(2),'String'))),
  uiwait(errordlg('Invalid entry'));
  set(he(2),'String',num2str(dampold(2)));
  return;
end

% Get the damping values from the edit box

damp=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];

% Make sure the damping value makes sense

if damp(2)<=0 || damp(2)>25,
  uiwait(errordlg('Unreasonable damping','Damping Value Error'));
  damp(2)=dampold(2);
elseif damp(2)<=damp(1),
  uiwait(errordlg('High value for damping must be higher than low value','Damping Value Error'));
  damp(2)=dampold(2);
end
damp(damp<0)=0.001;
damp(damp>25)=25;

% Update the form

set(he(2),'String',num2str(damp(2)));
setappdata(hf,'DampingOld',damp);

% Perform the fits

%uicbAddRootFreqFit(hf);
uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootLowFreqEdit(varargin)
% Edit lower bound for frequency

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditFreqLow');
he(2)=getappdata(hf,'HandleEditFreqHigh');
freqold=getappdata(hf,'FrequencyOld');

% Make sure the value is a valid number

if isempty(str2num(get(he(1),'String'))),
  uiwait(errordlg('Invalid entry'));
  set(he(1),'String',num2str(freqold(1)));
  return;
end

% Get the damping values from the edit box

freq=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];

% Make sure the frequency value makes sense

if freq(1)<ss.freqrangecc(1,1),
  uiwait(errordlg('Invalid frequency','Frequency Value Error'));
  freq(1)=freqold(1);
elseif freq(1)>freq(2),
  uiwait(errordlg('Low value for frequency must be lower than high value','Frequency Value Error'));
  freq(1)=freqold(1);
end
freq(freq<ss.freqrangecc(1,1))=ss.freqrangecc(1,1);
freq(freq>ss.freqrangecc(2,1))=ss.freqrangecc(2,1);

% Update the form

set(he(1),'String',num2str(freq(1)));
setappdata(hf,'FrequencyOld',freq);

% Perform the fits

uicbAddRootFreqFit(hf);
%uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootHighFreqEdit(varargin)
% Edit upper bound for frequency

global ss;

hf=varargin{1};
he(1)=getappdata(hf,'HandleEditFreqLow');
he(2)=getappdata(hf,'HandleEditFreqHigh');
freqold=getappdata(hf,'FrequencyOld');

% Make sure the value is a valid number

if isempty(str2num(get(he(2),'String'))),
  uiwait(errordlg('Invalid entry'));
  set(he(2),'String',num2str(freqold(2)));
  return;
end

% Get the damping values from the edit box

freq=[str2num(get(he(1),'String')) str2num(get(he(2),'String'))];

% Make sure the frequency value makes sense

if freq(2)>ss.freqrangecc(2,1),
  uiwait(errordlg('Invalid frequency','Frequency Value Error'));
  freq(2)=freqold(2);
elseif freq(2)<=freq(1),
  uiwait(errordlg('High value for frequency must be higher than low value','Frequency Value Error'));
  freq(2)=freqold(2);
end
freq(freq<ss.freqrangecc(1,1))=ss.freqrangecc(1,1);
freq(freq>ss.freqrangecc(2,1))=ss.freqrangecc(2,1);

% Update the form

set(he(2),'String',num2str(freq(2)));
setappdata(hf,'FrequencyOld',freq);

% Perform the fits

uicbAddRootFreqFit(hf);
%uicbAddRootDampFit(hf);

% Make the Calculate buttons active

hbc=getappdata(hf,'HandleCalculate');
set(hbc,'Enable','on');


%==========================================================================
function hf=uicbAddRootCalculate(varargin)
% Recalculate the best fit given the current settings

global ss;

type=varargin{1};
hf=varargin{2};

% Perform the fit

switch type
  case 'Freq',
    uicbAddRootFreqFit(hf);

  case 'Damp',
    uicbAddRootDampFit(hf);

end


%==========================================================================
function hf=uicbAddRootFreqFit(varargin)
% Perform a frequency fit on the given range

global ss;

hf=varargin{1};
hef(1)=getappdata(hf,'HandleEditFreqLow');
hef(2)=getappdata(hf,'HandleEditFreqHigh');
hed(1)=getappdata(hf,'HandleEditDampLow');
hed(2)=getappdata(hf,'HandleEditDampHigh');

% Get the frequency range for the fit

fmm=[str2num(get(hef(1),'String')) str2num(get(hef(2),'String'))];

% Get the damping value for the fit

hbd=getappdata(hf,'HandleBestDamp');
dmm=str2num(get(hbd,'String'))/100;
if isempty(dmm),
  dmm=mean([str2num(get(hed(1),'String')) str2num(get(hed(2),'String'))])/100;
end
if dmm==0,      % Should only hit this the first time through
  dmm=ss.corr.zeta;
end

numpts=ss.fit.freqpts;
Mcc=0;          % Work around a debugger issue with mcc
nref=length(ss.ref_coords);

wr=getappdata(hf,'Abscissa');
frf=getappdata(hf,'Ordinate');

% Perform the auto fit

smac_ffit_cauto;

% Store the maximum correlation coefficients

setappdata(hf,'CorMax',cormax(1,:));
setappdata(hf,'RefMax',max_ref);

% Do a polynomial curve-fit of the data and evaluate the points

ws=warning;
warning off MATLAB:polyfit:RepeatedPointsOrRescale;
p=polyfit(fvec(:),Mcc(:),2);
z=polyval(p,fvec(:));
warning(ws);

% Get the best-fit root by finding where the slope is 0

bf=roots(polyder(p));

% Calculate the correlation coefficient at this value

cc=polyval(p,bf);

% Update the best fit frequency edit box

hbf=getappdata(hf,'HandleBestFreq');
set(hbf,'String',num2str(bf));
setappdata(hbf,'CorrCoef',cc);

% Store the curve-fit data

setappdata(hf,'FreqFitData',[fvec(:) Mcc(:) z(:)]);

% Update the frequency plot

uicbAddRootUpdatePlot(hf,'Freq');


%==========================================================================
function hf=uicbAddRootDampFit(varargin)
% Perform a damping fit on the given range

global ss;

hf=varargin{1};
hef(1)=getappdata(hf,'HandleEditFreqLow');
hef(2)=getappdata(hf,'HandleEditFreqHigh');
hed(1)=getappdata(hf,'HandleEditDampLow');
hed(2)=getappdata(hf,'HandleEditDampHigh');

% Get the frequency value for the fit

hbf=getappdata(hf,'HandleBestFreq');
fmm=str2num(get(hbf,'String'));
if isempty(fmm),
  fmm=mean([str2num(get(hef(1),'String')) str2num(get(hef(2),'String'))]);
end

% Get the damping range for the fit

dmm=[str2num(get(hed(1),'String')) str2num(get(hed(2),'String'))]/100;

numpts=ss.fit.zetapts;
Mcc=0;          % Work around a debugger issue with mcc
nref=length(ss.ref_coords);

wr=getappdata(hf,'Abscissa');
frf=getappdata(hf,'Ordinate');

% Perform the auto fit

smac_dfit_cauto;

% Store the maximum correlation coefficients

setappdata(hf,'CorMax',cormax(1,:));
setappdata(hf,'RefMax',max_ref);

% Do a polynomial curve-fit of the data and evaluate the points

ws=warning;
warning off MATLAB:polyfit:RepeatedPointsOrRescale;
p=polyfit(dvec(:),Mcc(:),2);
z=polyval(p,dvec(:));
warning(ws);

% Get the best-fit root by finding where the slope is 0

bd=roots(polyder(p));

% Calculate the correlation coefficient at this value

cc=polyval(p,bd);

% Update the best fit damping edit box

hbd=getappdata(hf,'HandleBestDamp');
set(hbd,'String',num2str(bd*100));
setappdata(hbd,'CorrCoef',cc);

% Store the curve-fit data

setappdata(hf,'DampFitData',[dvec(:) Mcc(:) z(:)]);

% Update the frequency plot

uicbAddRootUpdatePlot(hf,'Damp');


%==========================================================================
function hf=uicbAddRootUpdatePlot(varargin)
% Update the plot

global ss;

hf=varargin{1};
type=varargin{2};

switch type,
  case 'Freq',
    ha=getappdata(hf,'HandleAxisFrequency');
    data=getappdata(hf,'FreqFitData');
    xlab='Frequency (Hz)';
    x=data(:,1);

  case 'Damp',
    ha=getappdata(hf,'HandleAxisDamping');
    data=getappdata(hf,'DampFitData');
    xlab='Damping (%)';
    x=data(:,1)*100;

  otherwise,
    error('Internal error:  Unknown plot type');
end

% Pull the data from the figure and update the plot

set(hf,'CurrentAxes',ha);

if ~isempty(data),
  hp=get(ha,'Children');

  if ~isempty(hp),      % Figure already exists
    hp=get(ha,'Children');
    set(hp(1),'XData',x,'YData',data(:,2),'LineStyle','-','Color',[0 0 1]);
    set(hp(2),'XData',x,'YData',data(:,3),'LineStyle','--','Color',[1 0 0]);
    set(ha,'XLim',x([1 end]));
  else                  % Draw new plots
    set(hf,'CurrentAxes',ha);
    hold on;
    plot(x,data(:,2),'b-');
    plot(x,data(:,3),'r--');
    xlim(x([1 end]));
    hold off;
    xlabel(xlab);
    ylabel('Correlation Coefficient');
  end
end
