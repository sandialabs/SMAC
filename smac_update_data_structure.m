function out=smac_update_data_structure(ss)
% SMAC_UPDATE_DATA_STRUCTURE  Update SMAC data structure to current version
%
% SS=SMAC_UPDATE_DATA_STRUCTURE(SS_OLD)
%
% SMAC_UPDATE_DATA_STRUCTURE attempts to update the supplied SMAC data
% structure to the current version.  If it cannot it returns [].

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
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add fields for low and high mode frequency range
%
%  06-Mar-2006 / ATA Engineering / Dan Hensley
%    o Add support for version 9
%
%==========================================================================

out=[];

narginchk(1,1);
if ~isstruct(ss)
	error('Input must be a SMAC structure');
end

sstmp=smac_create_data_structure;

%--------------------------------------------------------------------------

if isfield(ss,'version')

	str=sprintf('*** Updating SMAC structure from old version %%d to current version %d',sstmp.version);
	tmp=ss;
	oldver=[];

	if ~ismember(ss.version,[1 2 3 4 5 6 7 8 9])
		return;
	end

	% Already current version

	if ss.version==sstmp.version
		str='';
	end

%% Upgrade to version 2

	if ss.version==1
		% FIXME: Can't remember what changed here
		return;
	end

%% Upgrade to version 3

	if ss.version==2
		% FIXME: Can't remember what changed here
		return;
	end

%% Upgrade to version 4

	if ss.version==3
		oldver=get_old_version(oldver,ss.version);

		tmp.residuals=sstmp.residuals;
		if ss.residuals
			tmp.residuals.modeinertance.active=true;
			tmp.residuals.modecompliance.active=true;
		end

		tmp.version=4;
	end

%% Upgrade to version 5

	if ss.version==4
		oldver=get_old_version(oldver,ss.version);

		tmp.shape.rootsel=[];

		tmp.residuals.use=false;
		tmp.residuals.lowmode.residual=[];
		tmp.residuals.highmode.residual=[];
		tmp.residuals.modeinertance.residual=[];
		tmp.residuals.modecompliance.residual=[];
		tmp.residuals.inertance.residual=[];
		tmp.residuals.compliance.residual=[];

		tmp.version=5;
	end

%% Upgrade to version 6

	if ss.version==5
		oldver=get_old_version(oldver,ss.version);

		tmp.residuals=rmfield(tmp.residuals,'modeinertance');
		tmp.residuals=rmfield(tmp.residuals,'modecompliance');

		tmp.version=6;
	end

%% Upgrade to version 7

	if ss.version==6
		oldver=get_old_version(oldver,ss.version);

		field={'lowmode','highmode'};
		for k=1:length(field)
			res=tmp.residuals.(field{k});
			residual=res.residual;
			res=rmfield(res,'residual');
			res.frange=[];
			if res.active
				res.frange=ss.freqrange(:,2).';
			end
			res.residual=residual;
			tmp.residuals.(field{k})=res;
		end

		tmp.version=7;
	end

%% Upgrade to version 8

	if ss.version==7
		oldver=get_old_version(oldver,ss.version);

		if size(tmp.corr.corr,2)==1
			tmp.corr.corr(end,2)=0;
		end

		if size(tmp.fit.rootlist,2)==3
			tmp.fit.rootlist(end,4)=0;
		end

		tmpshp=tmp.shape;
		fname=fieldnames(tmp.shape);
		tmp.shape=sstmp.shape;
		
		for k=1:length(fname)
			tmp.shape.(fname{k})=ss.shape.(fname{k});
		end

		tmp.version=8;
	end

%% Upgrade to version 9

	if ss.version==8
		oldver=get_old_version(oldver,ss.version);

        tmp.shape.psi=[];
        tmp.version=9;
	end

else
	return;
end

% Return the new updated structure

if ~isempty(str)
	fprintf(str,oldver);
	fprintf('\n');
end

out=tmp;

%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================
function oldver=get_old_version(oldver,in)

if isempty(oldver),
    oldver=in;
end
