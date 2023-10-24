function h=smac_synth_mmif(h,frfsyn,frfexp,rootlist,xl)
% SMAC_SYNTH_MMIF  Synthesize MMIF and plot against experimental data
%
% H=SMAC_SYNTH_MMIF(H,FRFSYN,FRFEXP,ROOTLIST[,XL])
%
% SMAC_SYNTH_MMIF will compute and plot the multivariate mode indicator
% functions from the supplied analytical curve-fit FRF in the imat_fn
% FRFSYN and the experimental data in the imat_fn FRFEXP.  H is a figure
% handle for the previous MMIF plot.  If the handle exists and is valid,
% SMAC_SYNTH_MMIF will delete the figure, create the new one, and set the
% figure position to the old one. Otherwise it will simply create a new
% figure.  ROOTLIST is an Nx2 vector.  Column 1 contains the root numbers,
% and column 2 contains the root frequencies of all of the roots used to
% synthesize the FRF.  XL is an optional 1x2 vector containing the X axis
% range over which to calculate and display the MMIF.
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
%  17-Sep-2004 / ATA Engineering / Dan Hensley
%    o Initial creation (adapt from smac_synth_nmif)
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Make experimental vs analytical traces the same color, and differ
%      only in linestyle
%    o Put filename in figure title
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for new .mif fields
%    o Use previously generated analytical MIF if it exists
%
%  15-Jun-2005 / ATA Engineering / Dan Hensley
%    o Only display the roots on the plot if they were supplied
%
%==========================================================================

global ss;

% Argument checking

narginchk(4,5);
if ~isa(frfsyn,'imat_fn'),
  error('Synthesized FRF must be an imat_fn');
end
if ~isa(frfexp,'imat_fn'),
  error('Experimental FRF must be an imat_fn');
end
if ~exist('xl','var'),
  xl=[];
end
if ~isempty(xl) && (~isnumeric(xl) || numel(xl)~=2),
  error('Input frequency range must be a 1x2 numeric vector');
end

% Get indices into functions for range to process

fa=frfexp(1).abscissa;
if isempty(xl),
  ind=[1 frfexp(1).numberelements];
  xl=frfexp(1).abscissa(ind);
else
  [~,ind(1)]=min(abs(fa-xl(1)));
  [~,ind(2)]=min(abs(fa-xl(2)));
end

%--------------------------------------------------------------------------
% Compute MMIF on experimental data

if ~isempty(ss.mif.exp.mmif),
  mmifexp=ss.mif.exp.mmif;
else
  mmifexp=mmif(frfexp,'nofp','silent');
  ss.mif.exp.mmif=mmifexp;
end

%--------------------------------------------------------------------------
% Compute NMIF on analytical fit

if ~isempty(ss.mif.ana.mmif),
  mmifsyn=ss.mif.ana.mmif;
else
  mmifsyn=mmif(frfsyn,'nofp','silent');
  ss.mif.ana.mmif=mmifsyn;
end

% If the figure exists, remember the old position and delete the figure

pos=[];
if ~isempty(h) && ishandle(h.hf),
  set(h.hf,'Units','normalized');
  pos=get(h.hf,'Position');
  delete(h.hf);
end

% Plot the MMIF

h=plot([mmifexp(:),mmifsyn(:)],'xwin',xl);
hl=sort(get(gca,'children'));
nfunc=length(hl);
for k=1:nfunc/2,
      h.style('Linewidth',2,k);
      col=get(hl(k),'Color');
      h.style('Linestyle','--','Color',col,k+nfunc/2);
end

h.ylabel('MMIF Magnitude');
h.xlabel('Frequency (Hz)');

fname='';
if ~isempty(ss.filename), fname=[' (' ss.filename ')']; end
titlestr=['MMIF' fname];
titlestr=strrep(titlestr,'\','\\');
h.title(titlestr,'interpreter','none');
set(h.hf,'Name',titlestr);
%title('MMIF Comparing Synthesis to Data');

%legend('Experimental Data','Analytical Fit',0);
legstr={};
for k=1:length(mmifexp),
  legstr{end+1}=sprintf('Experimental Data (%d)',k);
end
for k=1:length(mmifsyn),
  legstr{end+1}=sprintf('Analytical Fit (%d)',k);
end
h.legend(legstr,0);

% Add points showing the roots 

if ~isempty(rootlist)
    fa=mmifsyn(1).abscissa;
    indr=zeros(size(rootlist,1),1);
    tt=cell(size(indr));
    for k=1:length(indr),
        [tmp,indr(k)]=min(abs(fa-rootlist(k,2)));
        tt{k}=sprintf('%d',rootlist(k,1));
    end

    set(0,'CurrentFigure',h.hf);
    set(h.hf,'CurrentAxes',h.ha(1));
    hold on
    plot(rootlist(:,2),mmifsyn(1).ordinate(indr),'r* ');
    text(rootlist(:,2),mmifsyn(1).ordinate(indr)*0.97,tt)
    hold off
end

% Remember the position of the old figure if necessary

if ~isempty(pos),
  % Make sure the figure appears on the screen
  if any(pos(1:2)>1) || any(pos(1:2)<0),
    pos(1:2)=[.5 .1];
  end
  set(h.hf,'Units','normalized','Position',pos);
else
  set(h.hf,'Units','normalized','Position',[0.53  0.52  0.45  0.43]);
end
