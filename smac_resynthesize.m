function frfsyn=smac_resynthesize(fe,shp,residuals,recalc_residuals)
% SMAC_RESYNTHESIZE  Resynthesize analytical FRF from shapes and residuals
%
%  FRFSYN=SMAC_RESYNTHESIZE(FE,SHP,RESIDUALS[,RECALC])
%
%  SMAC_RESYNTHESIZE performs an FRF synthesis based on the 
%  and returns the result as FRF.
%
%  FE is an imat_fn containing the experimental functions.  Analytical
%  functions for each of the functions in FE will be generated.  SHP is an
%  imat_shp containing all of the mode shapes to include in the synthesis.
%  RESIDUALS is a SMAC residual structure specifying which residual terms
%  should be included in the synthesis, as well as the residuals
%  themselves.  RECALC is an optional logical specifying whether to
%  recalculate the residual terms, or to use the ones in the RESIDUALS
%  structure.  The default is FALSE.
%
%  FRFSYN is an imat_fn containing the synthesized FRF.

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
%  10-Aug-2005 / ATA Engineering / Dan Hensley
%    o Catch the case where a coordinate was not found in the synthesis
%
%  16-Sep-2005 / ATA Engineering / Dan Hensley
%    o Add code to recalculate inertance and compliance residuals, but keep
%      it turned off
%
%  07-Mar-2006 / ATA Engineering / Dan Hensley
%    o Finish support for recalculating residuals, and add the toggle as an
%      optional input argument
%
%==========================================================================

%% Argument checking

narginchk(3,4);

if nargin<4
    recalc_residuals=false;
end
if ~islogical(recalc_residuals) || numel(recalc_residuals)~=1
    error('RECALC must be a logical scalar');
end

fa=fe;
frfsyn=imat_fn([]);

% Extract some initial data and parameters

nfunc=length(fe);
ref_coords=unique(allref(fe));
nrefs=length(ref_coords);

% Get FRF abscissa and ordinate values

foall=fe.ordinate;
frall=fe(1).abscissa';
wrall=2*pi*frall;

% Get a list of roots

wroot=2*pi*shp(:).frequency;
zroot=shp(:).damping;

numroots=length(wroot);
lambda = zeros(numroots,1);

% Determine whether we are doing a real or a complex fit

realcomplex=1;
if any(strcmpi(shp.shapetype,'complex'))
    realcomplex=2;
end

%==========================================================================
%% FRF Synthesis

% Calculate the kernel
% NOTE:  Assume modal mass = 1

j=sqrt(-1);

if realcomplex==1       % Real normal mode kernel
    kernel=realkernel(wrall,wroot,zroot);
else                    % Complex mode kernel
    kernel=complexkernel(wrall,wroot,zroot);
end

% Get the list of shape reference and response coordinates

shpref=imat_ctrace(shp.referencecoord);
shpres=imat_ctrace(shp.responsecoord);

%--------------------------------------------------------------------------
%% Loop over all FRF

for k=1:prod(size(fa))

    % Extract the reference coefficients

	shp1=zeros(1,length(shp));
	refcoord=imat_ctrace(fa(k).referencecoord);
	for m=1:length(shp)
		shptmp=shp(m);
		try
			if shpref(m)==refcoord
				refcoord=shpres(m);
				if sign(shpref(m))*sign(shpres(m))<0
					refcoord=-refcoord;
				end
			end
			shp1(m)=shptmp{refcoord};
		catch
			uiwait(errordlg(sprintf('Could not extract reference coordinate:\n\n%s\n\nSee Matlab window for more details.',lasterr)));
			return;
		end
	end
	
	% Extract the response coefficients
	
	try
		shp2=shp{fa(k).responsecoord};
	catch
		uiwait(errordlg(sprintf('Could not extract response coordinate:\n\n%s\n\nSee Matlab window for more details.',lasterr)));
		return;
	end

%-------------------------------------------
%% Complex mode FRF synthesis

    if realcomplex==2,

        top=(shp1.*shp2);
        top=[top conj(top)];
        top=top(ones(length(wrall),1),:);
        
        frf=sum(top.*kernel,2);
        
        % Add in low mode residual effects
        
        if residuals.lowmode.active
            wr=2*pi*residuals.lowmode.freq;
            zr=residuals.lowmode.damp;
            
            kern=complexkernel(wrall,wr,zr);
            res=residuals.lowmode.residual(k);
            res=[res conj(res)];
            frf=frf+sum(res(ones(length(wrall),1),:).*kern,2);
        end
        
        % Add in high mode residual effects
        
        if residuals.highmode.active
            wr=2*pi*residuals.highmode.freq;
            zr=residuals.highmode.damp;
            
            kern=complexkernel(wrall,wr,zr);
            res=residuals.highmode.residual(k);
            res=[res conj(res)].';
            
            frf=frf+kern*res;
        end

%===========================================
%% Real normal mode FRF synthesis
    else

        top=(shp1.*shp2).';
        frf=kernel*top;
        
        frf_to_fit=[];
        kern_to_fit=[];
        kern=[];
        
        % Add in low mode residual effects
        
        if residuals.lowmode.active
            wr=2*pi*residuals.lowmode.freq;
            zr=residuals.lowmode.damp;
            
            kernl=realkernel(wrall,wr,zr);
            kern=kernl;

            if ~recalc_residuals
                frf=frf+residuals.lowmode.residual(k)*kernl;
            else
                ind=find(frall>=residuals.lowmode.frange(1) & frall<=residuals.lowmode.frange(2));
                frf_to_fit=[frf_to_fit; real(fe(k).ordinate(ind)-frf(ind))];
                kern_to_fit(end+1:end+length(ind),end+1)=real(kernl(ind));
            end
        end
        
        % Add in high mode residual effects
        
        if residuals.highmode.active
            wr=2*pi*residuals.highmode.freq;
            zr=residuals.highmode.damp;

            kernh=realkernel(wrall,wr,zr);
            kern=[kern kernh];

            if ~recalc_residuals
                frf=frf+residuals.highmode.residual(k)*kernh;
            else
                ind=find(frall>=residuals.highmode.frange(1) & frall<=residuals.highmode.frange(2));
                frf_to_fit=[frf_to_fit; real(fe(k).ordinate(ind)-frf(ind))];
                kern_to_fit(end+1:end+length(ind),end+1)=real(kernh(ind));
            end
        end
        
        % Calculate the low and high mode terms if need be
        
        if recalc_residuals
            resvals=kern_to_fit\frf_to_fit;
            if ~isempty(resvals)
                frf=frf+kern*resvals;
            end
        end
    end

%----------------------------------------------------------------------
%% Add in residual effects

    % : Inertance :

    inert=0;
    if residuals.inertance.active
        if ~recalc_residuals    % Use existing value
            inert=inert+residuals.inertance.residual(k);
        else                    % Recalculate
            ind=find(frall>=residuals.inertance.freq(1) & frall<=residuals.inertance.freq(2));
            frf_to_fit=real(fe(k).ordinate(ind)-frf(ind));
            res=ones(length(ind),1);
            inert=inert + res\frf_to_fit;
        end
    end

    %     if residuals.modeinertance.active
    %         inert=inert+residuals.modeinertance.residual(k);
    %     end

    if inert~=0
        frf=frf+inert;
    end

    % : Compliance :

    compl=0;
    if residuals.compliance.active
        if ~recalc_residuals    % Use existing value
            compl=compl+residuals.compliance.residual(k);
        else                    % Recalculate
            ind=find(frall>=residuals.compliance.freq(1) & frall<=residuals.compliance.freq(2));
            frf_to_fit=real(fe(k).ordinate(ind)-frf(ind));
            res=-wrall(ind(:))'.^2;

            compl=compl + res\frf_to_fit;
        end
    end

    %     if residuals.modecompliance.active
    %         compl=compl+residuals.modecompliance.residual(k);
    %     end

    if compl~=0
        frf=frf+(compl.*(-wrall(:).^2));
    end

    % Store the result

    fa(k).ordinate=frf;

end

%--------------------------------------------------------------------------
%% Return the synthesized FRF

frfsyn=fa;

%==========================================================================
%% AUXILLIARY FUNCTIONS
%==========================================================================

%% realkernel
function kernel=realkernel(wrall,wroot,zroot)
% Compute the real kernel

j=sqrt(-1);

num=-wrall'.^2*ones(1,length(wroot));
denom=ones(length(wrall),1)*(wroot.^2)'+ num;
denom=denom+2*j*wrall'*(wroot'.*zroot');

kernel=num./denom;

%--------------------------------------------------------------------------
%% complexkernel
function kernel=complexkernel(wrall,wroot,zroot)
% Compute the complex kernel

j=sqrt(-1);

lambda = -wroot(:).*zroot(:) + j.*(wroot(:).*sqrt(1 - zroot(:).^2));
lambda = [ lambda; conj(lambda) ].';

num=-wrall'.^2*ones(1,2*length(wroot));
denom=j*wrall'*ones(1,2*length(wroot));
denom=denom - ones(length(wrall),1)*lambda;

kernel=num./denom;

%--------------------------------------------------------------------------
