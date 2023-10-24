function [frfsyn,residues,residuals]=smac_synthesize(rootlist,residuals)
% SMAC_SYNTHESIZE  Synthesize analytical fit
%
%  [FRFSYN,RESIDUES[,RESIDUALS]]=SMAC_SYNTHESIZE(ROOTLIST,RESIDUALS)
%
%  SMAC_SYNTHESIZE performs a synthesis of an analytical fit of the data
%  and returns the result as FRF.
%
%  ROOTLIST is an nx3 vector containing the roots to use in the synthesis.
%  The first column contains the root numbers, the second contains the
%  frequencies, and the third contains the damping values.  RESIDUALS is a
%  structure specifying what residual terms should be included.  RESIDUES

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
%  03-May-2004 / ATA Engineering / Dan Hensley
%    o Turn into a function
%    o Add comments and help to the function
%
%  27-May-2004 / ATA Engineering / Dan Hensley
%    o Clean up code and modify to use the new data structure
%
%  02-Jun-2004 / ATA Engineering / Dan Hensley
%    o Reorganize the code again to move plotting out
%
%  10-Sep-2004 / ATA Engineering / Dan Hensley
%    o Update for multiple references
%
%  03-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add support for low mode, high mode, and mode inertance and
%      compliance residue types
%
%  08-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add support for mode and general inertance and compliance in real
%      mode synthesis
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Return residuals structure as optional third output
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%    o Fix inertance and compliance residual calculations to include the
%      kernel instead of 0's
%    o Use auxilliary functions realkernel and complexkernel to simplify code
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Change lowmode and highmode residual fits to use the real portion of
%      the FRF (real modes), and the given frequency range rather than the
%      closest points to the resonant peak
%
%  10-Aug-2005 / ATA Engineering / Dan Hensley
%    o Rearrange real mode residual fitting to occur after the resonance
%      residue calculations and to occur inertance first, then compliance
%
%==========================================================================

global ss;

% Argument checking

narginchk(2,2);

% Extract some initial data and parameters

nfunc=length(ss.fe);
nrefs=length(ss.ref_coords);

% Get FRF abscissa and ordinate values

foall=ss.fe(:).ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:);
frall=ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2))';
wrall=2*pi*frall;

% Get frequency range to synthesize and plot

frdisp=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2))';
freqrange=[min(frdisp) max(frdisp)];
wrdisp=2*pi*frdisp;

wroot=2*pi*rootlist(:,2);
zroot=rootlist(:,3);

numroots=size(rootlist,1);
lambda = zeros(numroots,1);

% Find parts of frfs to fit from root table and put the
% omegas to fit in a vector and frfs to fit in matrix

%==========================================================================
% Complex mode FRF synthesis

if ss.realcomplex==2
    nl=ss.corr.nl;

    % Preallocate memory

    wfit=zeros(1,nl*numroots);
    frf_to_fit=zeros(nl*numroots,nfunc*nrefs);

    for k=1:numroots

        % Get the closest spectral line to the root frequency

        [~,sind]=min(abs(frall-rootlist(k,2)));

        % Make sure this spectral line isn't too close to the edge

        if ( (sind - (round(nl/2) - 1)) <= 0 )
            sind = round(nl/2);
        end
        if ( (sind + (nl-round(nl/2))) >= max(size(frall)) )
            sind = max(size(frall)) - (nl-round(nl/2));
        end

        % Calculate

        wfit(1,(k-1)*nl+(1:nl))=wrall(sind-(round(nl/2)-(1:nl)));
        frf_to_fit(((k-1)*nl+(1:nl)),:)=foall(sind-(round(nl/2)-(1:nl)),:);
    end

    % Calculate the synthesized FRF kernel

    kernel=complexkernel(wfit,wroot,zroot);

    % Calculate residues from frfs including negative frequency values

    kernel1=complexkernel(-wfit,wroot,zroot);

    %  Calculate kernel for display of frf synthesis

    kernel_disp=complexkernel(wrdisp,wroot,zroot);

    %--------------------------------------------------------
    %  Add in residual effects as requested

    res=[];
    res1=[];
    res_disp=[];

    % Add low mode

    if residuals.lowmode.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.lowmode.frange(1) & frall<=residuals.lowmode.frange(2));
        is=size(frf_to_fit,1)+1;
        ie=is+length(ind)-1;

        % Expand the FRF and kernel matrices to include these lines

        frf_to_fit=[frf_to_fit; foall(ind,:)];
        kernel(is:ie,1:end)=complexkernel(wrall(ind),wroot,zroot);
        kernel1(is:ie,1:end)=complexkernel(-wrall(ind),wroot,zroot);

        % Calculate residual matrix terms

        res(is:ie,end+1:end+2)=complexkernel(wrall(ind),2*pi*residuals.lowmode.freq,residuals.lowmode.damp);
        res1(is:ie,end+1:end+2)=complexkernel(-wrall(ind),2*pi*residuals.lowmode.freq,residuals.lowmode.damp);
        res_disp=[res_disp complexkernel(wrdisp,2*pi*residuals.lowmode.freq,residuals.lowmode.damp)];
    end

    % Add high mode

    if residuals.highmode.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.highmode.frange(1) & frall<=residuals.highmode.frange(2));
        is=size(frf_to_fit,1)+1;
        ie=is+length(ind)-1;

        % Expand the FRF and kernel matrices to include these lines

        frf_to_fit=[frf_to_fit; foall(ind,:)];
        kernel(is:ie,1:end)=complexkernel(wrall(ind),wroot,zroot);
        kernel1(is:ie,1:end)=complexkernel(-wrall(ind),wroot,zroot);

        % Calculate residual matrix terms

        res(is:ie,end+1:end+2)=complexkernel(wrall(ind),2*pi*residuals.highmode.freq,residuals.highmode.damp);
        res1(is:ie,end+1:end+2)=complexkernel(-wrall(ind),2*pi*residuals.highmode.freq,residuals.highmode.damp);
        res_disp=[res_disp complexkernel(wrdisp,2*pi*residuals.highmode.freq,residuals.highmode.damp)];
    end

    %     % Add mode inertance (don't include low mode and high mode here)
    %
    %     if residuals.modeinertance.active
    %         res=[res ones(length(wfit),1)];
    %         res_disp=[res_disp ones(length(wrdisp),1)];
    %         res(:,end)=res(:,end).*mult_fit;
    %     end

    %     % Add mode compliance (don't include low mode and high mode here)
    %
    %     if residuals.modecompliance.active
    %         res=[res -wfit(:).^2];
    %         res_disp=[res_disp -wrdisp(:).^2];
    %         res(:,end)=res(:,end).*mult_fit;
    %     end

    % Add general inertance

    if residuals.inertance.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.inertance.freq(1) & frall<=residuals.inertance.freq(2));
        is=size(frf_to_fit,1)+1;
        ie=is+length(ind)-1;

        % Expand the FRF and kernel matrices to include these lines

        frf_to_fit=[frf_to_fit; foall(ind,:)];
        kernel(is:ie,1:end)=complexkernel(wrall(ind),wroot,zroot);
        kernel1(is:ie,1:end)=complexkernel(-wrall(ind),wroot,zroot);

        % Calculate residual matrix terms

        res(is:ie,end+1)=ones(length(ind),1);
        res1(is:ie,end+1)=res(is:ie,end);
        res_disp=[res_disp ones(length(wrdisp),1)];
    end

    % Add general compliance

    if residuals.compliance.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.compliance.freq(1) & frall<=residuals.compliance.freq(2));
        is=size(frf_to_fit,1)+1;
        ie=is+length(ind)-1;

        % Expand the FRF and kernel matrices to include these lines

        frf_to_fit=[frf_to_fit; foall(ind,:)];

        kernel(is:ie,1:end)=complexkernel(wrall(ind),wroot,zroot);
        kernel1(is:ie,1:end)=complexkernel(-wrall(ind),wroot,zroot);

        % Calculate residual matrix terms

        res(is:ie,end+1)=-wrall(ind(:)).^2;
        res1(is:ie,end+1)=res(is:ie,end);
        res_disp=[res_disp -wrdisp(:).^2];
    end

    %--------------------------------------------------------
    %  Calculate residues

    kernel=[[kernel; kernel1] [res; res1]];
    kernel_disp=[kernel_disp res_disp];

    frf_to_fit = [frf_to_fit; conj(frf_to_fit)];

    residues=kernel\frf_to_fit;

%==========================================================================
%==========================================================================
% Real normal mode FRF synthesis

else

    % Preallocate memory

    nl=4;
    wfit=zeros(1,nl*numroots);
    ifrf_to_fit=zeros(nl*numroots,nfunc*nrefs);
    rfrf_to_fit=zeros(nl*numroots,nfunc*nrefs);

    % Loop over all of the roots

    for k=1:numroots

        % Get the closest spectral line to the root frequency

        [~,sind]=min(abs(frall-rootlist(k,2)));

        % Make sure this spectral line isn't too close to the edge

        if sind<2, sind=2; end
        if sind+2>length(frall), sind=length(frall)-2; end

        % Calculate

        wfit(4*k-3:4*k)=wrall(sind-1:sind+2);
        ifrf_to_fit(4*k-3:4*k,:)=imag(foall(sind-1:sind+2,:));
        rfrf_to_fit(4*k-3:4*k,:)=real(foall(sind-1:sind+2,:));
    end

    residues=zeros(numroots,nfunc);

    % Calculate the synthesized FRF kernel

    kernel=realkernel(wfit,wroot,zroot);

    %  Calculate kernel for display of frf synthesis

    kernel_disp=realkernel(wrdisp,wroot,zroot);

    %--------------------------------------------------------
    %  Add in residual effects as requested

    res=[];
    res_disp=[];
    usekernelreal=false;
    kerneladd=[];
	wrootadd=[];
	zrootadd=[];

    % Add low mode

    if residuals.lowmode.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.lowmode.frange(1) & frall<=residuals.lowmode.frange(2));
        is=size(ifrf_to_fit,1)+1;
        ie=is+length(ind)-1;

        ifrf_to_fit=[ifrf_to_fit; real(foall(ind,:))];
        kerneladd=[kerneladd; real(realkernel(wrall(ind),wroot,zroot))];

        % Calculate residual matrix terms

        res(is:ie,end+1)=real(realkernel(wrall(ind),2*pi*residuals.lowmode.freq,residuals.lowmode.damp));
        res_disp=[res_disp realkernel(wrdisp,2*pi*residuals.lowmode.freq,residuals.lowmode.damp)];
		
		% Append frequency and damping to additional storage terms
		
		wrootadd(end+1)=2*pi*residuals.lowmode.freq;
		zrootadd(end+1)=residuals.lowmode.damp;
    end

    % Add high mode

    if residuals.highmode.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.highmode.frange(1) & frall<=residuals.highmode.frange(2));
        is=size(ifrf_to_fit,1)+1;
        ie=is+length(ind)-1;

        ifrf_to_fit=[ifrf_to_fit; real(foall(ind,:))];
        kerneladd=[kerneladd; real(realkernel(wrall(ind),wroot,zroot))];

        % Calculate residual matrix terms

        res(is:ie,end+1)=real(realkernel(wrall(ind),2*pi*residuals.highmode.freq,residuals.highmode.damp));
        res_disp=[res_disp realkernel(wrdisp,2*pi*residuals.highmode.freq,residuals.highmode.damp)];
		
		% Append frequency and damping to additional storage terms
		
		wrootadd(end+1)=2*pi*residuals.highmode.freq;
		zrootadd(end+1)=residuals.highmode.damp;
    end

    %     % Add mode inertance (don't include low mode and high mode here)
    %
    %     if residuals.modeinertance.active
    %         is=length(wfit)+1;
    %         ie=is+length(wfit)-1;
    %         res(is:ie,end+1)=ones(length(wfit),1);
    %         res_disp=[res_disp ones(length(wrdisp),1)];
    %         res(is:ie,end)=res(is:ie,end).*mult_fit;
    %         ifrf_to_fit(is:ie,:)=rfrf_to_fit;
    %         usekernelreal=true;
    %     end

    %     % Add mode compliance (don't include low mode and high mode here)
    %
    %     if residuals.modecompliance.active
    %         is=length(wfit)+1;
    %         ie=is+length(wfit)-1;
    %         res(is:ie,end+1)=-wfit(:).^2;
    %         res_disp=[res_disp -wrdisp(:).^2];
    %         res(is:ie,end)=res(is:ie,end).*mult_fit;
    %         ifrf_to_fit(is:ie,:)=rfrf_to_fit;
    %         usekernelreal=true;
    %     end
	
    %--------------------------------------------------------
    %  Calculate residues for what we have so far from imaginary part of frfs only

    kerneli=imag(kernel);
    kernelr=real(kernel);

    if ~usekernelreal
        kernelr=[];
    end

    kernel=[[kerneli; kernelr; kerneladd] res];
    residues=kernel\ifrf_to_fit;

	res_add=[];		% Additional residues
	
    %--------------------------------------------------------
    % Add general inertance

    if residuals.inertance.active

        % Calculate where the new data fits in

        ind=find(frall>=residuals.inertance.freq(1) & frall<=residuals.inertance.freq(2));

        % Subtract the analytical FRF containing the roots from the experimental
		
		frf_rem=realkernel(wrall(ind),[wroot(:); wrootadd(:)],[zroot(:); zrootadd(:)])*residues(:,:);
		frf_to_fit=real(foall(ind,:)-frf_rem);

        % Calculate residual matrix terms
		
		res=ones(length(ind),1);
        res_disp=[res_disp ones(length(wrdisp),1)];
		
		res_add=res\frf_to_fit;
    end

    %--------------------------------------------------------
    % Add general compliance

    if residuals.compliance.active
		
		% Calculate where the new data fits in

        ind=find(frall>=residuals.compliance.freq(1) & frall<=residuals.compliance.freq(2));

		% Subtract the analytical FRF containing the roots from the experimental

		frf_rem=realkernel(wrall(ind),[wroot(:); wrootadd(:)],[zroot(:); zrootadd(:)])*residues(:,:);
		frf_to_fit=foall(ind,:)-frf_rem;
		if ~isempty(res_add)
			frf_to_fit=frf_to_fit - ones(length(ind),1)*res_add;
		end	
		frf_to_fit=real(frf_to_fit);
		
		% Calculate residual matrix terms

        res=-wrall(ind(:)).^2;
		res=res(:);
        res_disp=[res_disp -wrdisp(:).^2];

		res_add(end+1,:)=res\frf_to_fit;
	end

    %--------------------------------------------------------
    %  Now put together all of the terms
	
	residues=[residues; res_add];
    kernel_disp=[kernel_disp res_disp];

end

%--------------------------------------------------------------------------
% Store the residuals in the output structure

ind=ss.realcomplex*numroots+1;

if residuals.lowmode.active
    residuals.lowmode.residual=residues(ind,:);
    ind=ind+ss.realcomplex;
end

if residuals.highmode.active
    residuals.highmode.residual=residues(ind,:);
    ind=ind+ss.realcomplex;
end

% if residuals.modeinertance.active
%     residuals.modeinertance.residual=residues(ind,:);
%     ind=ind+1;
% end

% if residuals.modecompliance.active
%     residuals.modecompliance.residual=residues(ind,:);
%     ind=ind+1;
% end

if residuals.inertance.active
    residuals.inertance.residual=residues(ind,:);
    ind=ind+1;
end

if residuals.compliance.active
    residuals.compliance.residual=residues(ind,:);
    ind=ind+1;
end

%--------------------------------------------------------------------------
% Calculate the synthesized FRF

frf_fit_all=kernel_disp*residues(:,:);

% Put the synthesized FRF into an imat_fn

frfsyn=ss.fe;			% Take attributes from experimental FRF
frfsyn.ordinate=frf_fit_all;
frfsyn.abscissa=frdisp;

%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================

function kernel=realkernel(wrall,wroot,zroot)
% Compute the real kernel

j=sqrt(-1);

num=-wrall'.^2*ones(1,length(wroot));
denom=ones(length(wrall),1)*(wroot.^2)'+ num;
denom=denom+2*j*wrall'*(wroot'.*zroot');

kernel=num./denom;

%--------------------------------------------------------------------------
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
