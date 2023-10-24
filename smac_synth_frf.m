function h=smac_synth_frf(h,frfsyn,frfexp,rootlist,freqrange,mode,haveresiduals)
% SMAC_SYNTH_FRF  Synthesize FRF and plot against experimental data
%
% H=SMAC_SYNTH_FRF(H,FRFSYN,FRFEXP,ROOTLIST,FREQRANGE,XL,MODE,HAVERESIDUALS)
%
% SMAC_SYNTH_FRF will overlay the FRF from the supplied analytical curve-
% fit FRF in the imat_fn FRFSYN and the experimental data in the imat_fn
% FRFEXP.  H is a figure handle for the previous FRF plot.  If the handle
% exists and is valid, SMAC_SYNTH_FRF will delete the figure, create the
% new one, and set the figure position to the old one. Otherwise it will
% simply create a new figure.  ROOTLIST is an Nx2 vector.  Column 1 contains
% the root numbers, and column 2 contains the root frequencies of all of the
% roots used to synthesize the FRF.  XL is a 1x2 vector containing the X axis
% range over which to display the FRF.  MODE specifies whether the fit is
% real or complex (1=real, 2=complex).  HAVERESIDUALS is a logical that
% specifies whether residual terms are included in the analytical FRF.  If
% MODE specifies real, the imaginary portion of the FRF will be the default
% view on the plot.  If residual terms are included, both the real and
% imaginary portions of the FRF will be displayed.  If complex, the modulus
% and phase are displayed.
%
% H is the imat_fnplot handle of the created figure.

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
%  13-May-2004 / ATA Engineering / Dan Hensley
%    o Initial release
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Put filename in figure title
%
%  02-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make sure title appears above uppermost graph
%
%  17-Feb-2005 / ATA Engineering / Dan Hensley
%    o Fix bug where experimental and analytical plots got swapped after
%      hitting Prev or Next
%
%  15-Jun-2005 / ATA Engineering / Dan Hensley
%    o Only display the roots on the plot if they were supplied
%
%  27-Jul-2005 / ATA Engineering / Dan Hensley
%    o When uiselecting FRF, preselect the currently displayed FRF
%
%  10-Aug-2005 / ATA Engineering / Dan Hensley
%    o When fitting real modes with residual terms, plot both real and
%      imaginary of FRF
%
%  11-Sep-2006 / ATA Engineering / Dan Hensley
%    o Update to use IMAT v2.0 new plot handling capabilities
%
%==========================================================================

% Argument checking

narginchk(7,7);

if ~isa(frfsyn,'imat_fn');
  error('Synthesized FRF input must be an imat_fn');
end
if ~isa(frfexp,'imat_fn');
  error('Experimental FRF input must be an imat_fn');
end
if ~isempty(freqrange) && (~isnumeric(freqrange) || numel(freqrange)~=2),
  error('Frequency range must be a 1x2 numeric vector');
end

freqrange=sort(freqrange);

%--------------------------------------------------------------------------
% Build the plot

% If the figure exists, remember the old position and delete the figure

pos=[];
if ~isempty(h) && ishandle(h.hf),
  set(h.hf,'Units','normalized');
  pos=get(h.hf,'Position');
  delete(h.hf);
end

% Build the plot

ind=1;

% Build the root list plot

fa=frfsyn(ind).abscissa;
indr=zeros(size(rootlist,1),1);
tt=cell(size(indr));
for k=1:length(indr),
  [~,indr(k)]=min(abs(fa-rootlist(k,2)));
  tt{k}=sprintf('%d',rootlist(k,1));
end

% Set the frequency range if it's empty

if isempty(freqrange)
    freqrange=fa([1 end]);
end

if ~isempty(rootlist)
    froot=frfsyn(ind);
    froot.responsecoord='1root';
    froot.abscissa=rootlist(:,2);
    froot.ordinate=frfsyn(ind).ordinate(indr);
else
    froot=imat_fn([]);
end

% Get the default display type depending on type of fit

if mode==1,
    ctype='i';	% Imaginary for real normal modes with no residual terms
    if haveresiduals
        ctype='ri';		% Real and imaginary if residual terms are included
    end
    yscale={'yscale' 'lin'};
else
	ctype='m';	% Magnitude for complex modes
    yscale={};
end

h=plot(frfexp(ind),frfsyn(ind),froot,'xwin',freqrange,'complex',ctype,yscale{:});

% Put FRF data, current plot into appdata

setappdata(h.hf,'FRF_Synthesized',frfsyn);
setappdata(h.hf,'FRF_Experimental',frfexp);
setappdata(h.hf,'FRF_Index',ind);
setappdata(h.hf,'FRF_Roots',froot);
setappdata(h.hf,'FRF_Roots_Ind',indr);
setappdata(h.hf,'FRF_Roots_Text',tt);

% Update the root plot format

smac_frf_fig_format_root_plot(h);

% Put buttons and uimenus on

setappdata(h.hf,'Plot',h);
% setappdata(h.hf,'FigPrevHandle',@smac_frf_fig_prev);
% setappdata(h.hf,'FigNextHandle',@smac_frf_fig_next);
% setappdata(h.hf,'FigSelectHandle',@smac_frf_fig_select);
setappdata(h.hf,'FigPlotHandle',@smac_frf_fig_update_plot);
setappdata(h.hf,'FigUpdateHandle',@smac_frf_fig_update_labels);
setappdata(h.hf,'FigFormatRootsHandle',@smac_frf_fig_format_root_plot);

% Add Next and Previous uimenus

uimenu(h.hf, ...
       'Label','<<', ...
       'Callback',@smac_frf_fig_prev);
uimenu(h.hf, ...
       'Label','?', ...
       'Callback',@smac_frf_fig_select);
uimenu(h.hf, ...
       'Label','>>', ...
       'Callback',@smac_frf_fig_next);

% Add Next and Previous buttons

un=get(h.hf,'Units');
uicontrol(double(h.hf), ...
         'Units','normalized', ...
         'Position',[.01 .01 .05 .05], ...
         'String','<<', ...
         'TooltipString','Overlay the previous FRF', ...
         'Callback',@smac_frf_fig_prev);
uicontrol(double(h.hf), ...
         'Units','normalized', ...
         'Position',[.07 .01 .05 .05], ...
         'String','?', ...
         'TooltipString','Select the FRF to display from the list', ...
         'Callback',@smac_frf_fig_select);
uicontrol(double(h.hf), ...
         'Units','normalized', ...
         'Position',[.13 .01 .05 .05], ...
         'String','>>', ...
         'TooltipString','Overlay the next FRF', ...
         'Callback',@smac_frf_fig_next);
set(h.hf,'Units',un);

% Put on display information

smac_frf_fig_update_labels(h);

% Remember the position of the old figure if necessary

if ~isempty(pos),
  set(h.hf,'Units','normalized','Position',pos);
else
  set(h.hf,'Units','Normalized', ...
        'Position',[0.01 0.52 0.45 0.43]);
end

%==========================================================================
%% _____ CALLBACKS ________________________________________________________
%==========================================================================

%% smac_frf_fig_prev
function smac_frf_fig_prev(varargin)

h = getappdata(gcbf,'Plot');
hf = h.hf;

frf=getappdata(hf,'FRF_Synthesized');
ind=getappdata(hf,'FRF_Index');

% Decrement the index by 1

ind=ind-1;
if ind<1,
  ind=prod(size(frf));
end
setappdata(hf,'FRF_Index',ind);

% Update the plot and labels

feval(getappdata(hf,'FigPlotHandle'),h);
feval(getappdata(hf,'FigUpdateHandle'),h);

%--------------------------------------------------------------------------
%% smac_frf_fig_next
function smac_frf_fig_next(varargin)

h = getappdata(gcbf,'Plot');
hf = h.hf;

frf=getappdata(hf,'FRF_Synthesized');
ind=getappdata(hf,'FRF_Index');

% Increment the index by 1

ind=ind+1;
if ind>prod(size(frf)),
  ind=1;
end
setappdata(hf,'FRF_Index',ind);

% Update the plot and labels

feval(getappdata(hf,'FigPlotHandle'),h);
feval(getappdata(hf,'FigUpdateHandle'),h);

%--------------------------------------------------------------------------
%% smac_frf_fig_select
function smac_frf_fig_select(varargin)

h = getappdata(gcbf,'Plot');
hf = h.hf;

frf=getappdata(hf,'FRF_Synthesized');

% Prompt user

ind=getappdata(hf,'FRF_Index');
notdone=true;

while notdone,
  [f,ind]=uiselect(frf,ind);
  if isempty(ind), return; end

  notdone=false;
  if length(ind)>1,
    uiwait(errordlg('You may only select one from the list'));
	notdone=true;
	ind=ind(1);
  end
end

% Store the selected index

setappdata(hf,'FRF_Index',ind);

% Update the plot and labels

feval(getappdata(hf,'FigPlotHandle'),h);
feval(getappdata(hf,'FigUpdateHandle'),h);

%--------------------------------------------------------------------------
%% smac_frf_fig_update_labels
function smac_frf_fig_update_labels(h)
% Update the title and legend

frf=getappdata(h.hf,'FRF_Synthesized');
ind=getappdata(h.hf,'FRF_Index');

refc = frf(ind).referencecoord;
if ~iscell(refc), refc = {refc}; end
resc = frf(ind).responsecoord;
if ~iscell(resc), resc = {resc}; end
id1 = frf(ind).idline1;
if ~iscell(id1), id1 = {id1}; end

h.title(sprintf('%d: (%s,%s) "%s"',ind,refc{:},resc{:},id1{:}),'interpreter','none');

h.legend('Experimental Data','Analytical Fit',0)
h.xlabel('Frequency (Hz)');

%--------------------------------------------------------------------------
%% smac_frf_fig_update_plot
function smac_frf_fig_update_plot(h)
% Update the plot

frfsyn=getappdata(h.hf,'FRF_Synthesized');
frfexp=getappdata(h.hf,'FRF_Experimental');
froot=getappdata(h.hf,'FRF_Roots');
indr=getappdata(h.hf,'FRF_Roots_Ind');
ind=getappdata(h.hf,'FRF_Index');

% Update the root plot information

if ~isempty(froot)
    froot.ordinate=frfsyn(ind).ordinate(indr);
end
setappdata(h.hf,'FRF_Roots',froot);

% Replace the functions on the plot

h.replace([frfexp(ind); frfsyn(ind); froot]);

feval(getappdata(h.hf,'FigFormatRootsHandle'),h);

%--------------------------------------------------------------------------
%% smac_frf_fig_format_root_plot
function smac_frf_fig_format_root_plot(h)
% Update the plot

global ss;

if length(h.f) > 2
    h.style('Marker','*','LineStyle','none',3)
end

fname='';
if ~isempty(ss.filename), fname=[' (' ss.filename ')']; end
set(h.hf,'Name',['FRF' fname]);

% Draw the root number text--figure out the Y values from the plot

tt=getappdata(h.hf,'FRF_Roots_Text');
froot=getappdata(h.hf,'FRF_Roots');
ht=getappdata(h.hf,'HandleRootText');

if ~isempty(froot)
%     ha = get(h.hf,'CurrentAxes');
%     my=max(abs(diff(get(ha,'YLim'))));
%     if strcmpi(get(ha,'YScale'),'log'), my=log10(my); end

    fo=h.f(end).ordinate;
    fo=fo(1:length(tt));
    fo=fo+0.2*abs(fo);

    delete(ht(ishandle(ht)));
    ht = text(froot.abscissa,fo,tt);
    setappdata(h.hf,'HandleRootText',ht);
end
