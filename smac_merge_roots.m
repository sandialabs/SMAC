function merge_roots(fname)
% SMAC_MERGE_ROOTS  Merge roots from an ASH into current SMAC session
%
% SMAC_MERGE_ROOTS([FNAME])
%
% SMAC_MERGE_ROOTS will read shapes from an ASH file and merge them into
% the SMAC global structure ss.  If the optional filename FNAME is not
% supplied, the user will be prompted for the filename with a graphical
% dialog.
%
% The SMAC synthesize form must be open for this function to continue.  It
% will insert the roots from the ASH file into the SMAC data structure.  It
% ignores any residual modes in the ASH file (marked with damping of
% -100%).  If the imported roots lie outside the frequency range used to
% fit the current set, the frequency range will be expanded to 10 Hz around
% the frequency extremes.

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
%
%  31-Aug-2005 / ATA Engineering / Dan Hensley
%    o Expand the frequency range to the full range when reading in, rather
%      than 10 Hz around the root frequency extremes
%
%==========================================================================

global ss

% Check input arguments

if ~exist('fname','var')
	fname='*.ash';
end

% Find the form

hf=findobj(0,'Tag','SMAC_SynthForm');
if isempty(hf)
	uiwait(errordlg(sprintf(['Could not find the SMAC synthesis form\n' ...
		                     'Please make sure this form is open before using this function'])));
	return;
end

% Get the filename

if ~isempty(findstr(fname,'*')) || ~isempty(findstr(fname,'?'))
	fname=imat_getfile('*.ash','Select ASH file to merge');
	if isnumeric(fname)
		return;
	end
end

% Try reading in the file

try
	shp=readadf(fname);
catch
	uiwait(errordlg(lasterr));
	return;
end

%--------------------------------------------------------------------------
% Fill in the ss structure

shp=shp(shp.damping~=-1);	% Get rid of residual shapes

nshp=length(shp);
if nshp<1, return;, end

% Insert roots into the 
ss.fit.rootlistorig(end+nshp,end)=0;
ss.fit.rootlist=[ss.fit.rootlist; [shp.frequency shp.damping zeros(nshp,2)]];
ss.fit.corr_ref(end+nshp,end)=0;

% Sort frequencies

[tmp,ind]=sort(ss.fit.rootlist(:,1));
ss.fit.rootlistorig=ss.fit.rootlistorig(ind,:);
ss.fit.rootlist=ss.fit.rootlist(ind,:);
ss.fit.corr_ref=ss.fit.corr_ref(ind,:);

% Now expand the frequency ranges as necessary

freqbuffer=10;
if ss.freqrange(1,1)>ss.fit.rootlist(1,1)
	%ss.freqrange(1,:)=getnewfreqrange(ss.fe(1),ss.fit.rootlist(1,1)-freqbuffer);
	ss.freqrange(1,:)=getnewfreqrange(ss.fe(1),ss.fe(1).abscissa(1));
end
if ss.freqrange(2,1)<ss.fit.rootlist(end,1)
	%ss.freqrange(2,:)=getnewfreqrange(ss.fe(1),ss.fit.rootlist(end,1)+freqbuffer);
	ss.freqrange(2,:)=getnewfreqrange(ss.fe(1),ss.fe(1).abscissa(end));
end

if ss.freqrangecc(1,1)>ss.fit.rootlist(1,1)
	%ss.freqrangecc(1,:)=getnewfreqrange(ss.fe(1),ss.fit.rootlist(1,1)-freqbuffer);
	ss.freqrangecc(1,:)=getnewfreqrange(ss.fe(1),ss.fe(1).abscissa(1));
end
if ss.freqrangecc(2,1)<ss.fit.rootlist(end,1)
	%ss.freqrangecc(2,:)=getnewfreqrange(ss.fe(1),ss.fit.rootlist(end,1)+freqbuffer);
	ss.freqrangecc(2,:)=getnewfreqrange(ss.fe(1),ss.fe(1).abscissa(end));
end

ss.shape.rootsel=[];

% Clear out any previously synthesized data

smac_GUI_synth('clearanalysisdata',hf);

% Delete the form and reopen it

close(hf);
smac_GUI_synth;

%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================
function out=getnewfreqrange(f,newfreq,type)

fa=f(1).abscissa;
[tmp,ind]=min(abs(fa-newfreq));
if fa(ind)==0,
	ind=ind+1;
end
out=[fa(ind) ind];


