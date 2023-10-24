% SMAC_FFIT_CAUTO (script)  Perform frequency curve-fit
%
% SMAC_FFIT_CAUTO is called from SMAC_AUTOFIT.  It requires that some
% variables are set up before calling it.  It iterates over a vector of
% frequencies and a fixed damping value to determine the frequency that
% maximizes the correlation coefficient.

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
%  26-May-2004 / ATA Engineering / Dan Hensley
%    o Initial creation; adapted from ffit_cauto.m
%    o Renamed variables, reorganized, added comments
%
%  10-Sep-2004 / ATA Engineering / Dan Hensley
%    o Update for multiple references
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Rearrange CC code to look like smac_corr (loop over each reference)
%
%==========================================================================

% Get the frequency vector over which to iterate (single damping)

fvec=linspace(fmm(1),fmm(2),numpts);
dvec=dmm(1);

wvec=2*pi*fvec';

% Calculate the transfer functions for this frequency vector

nlines=length(wr);
nfunc=length(wvec);
nref=length(ss.ref_coords);

H = -(ones(nfunc,1)*wr.^2 ./ ...
	(wvec.^2*ones(1,nlines) + j*2*dvec*wvec*wr - ones(nfunc,1)*wr.^2));

% Loop over each reference

cormax=zeros(2,nref);
for k=1:nref

	% Compute the reciprocal modal vector Psi and Hp (SDOF FRF)

	if ss.realcomplex==2,	% Complex
		psi = ss.pinv(:,:,k)*H.';
		Hp = frf(:,:,k)*psi;
	else			% Real normal
		Hs = [real(H) imag(H)];
		psi = ss.pinv(:,:,k)*Hs';
		Hp = frf(:,:,k)*psi;
	end

	% Initialize variables

	Mcc=zeros(nfunc,1);
	nl=ss.corr.nl;

	% Compute the correlation coefficients

	ws=warning;
	warning off MATLAB:divideByZero;

	for i=1:nfunc
		[Y,ind] = max(abs(H(i,:)));
		if ind+nl > nlines		% On the right side of the FRF
			R = corrcoef([abs(H(i,max(1,ind-nl):nlines)).' abs(Hp(max(1,ind-nl):nlines,i))]);
		elseif ind-nl <= 0		% On the left side of the FRF
			R = corrcoef([abs(H(i,1:min(ind+nl,nlines))).' abs(Hp(1:min(ind+nl,nlines),i))]);
		else				% In the middle
			R = corrcoef([abs(H(i,max(1,ind-nl):min(ind+nl,nlines))).' abs(Hp(max(1,ind-nl):min(ind+nl,nlines),i))]);
		end

		Mcc(i,1) = R(1,2);
	end
	warning(ws);

	% Store

	[cormax(1,k),itmp]=max(Mcc);
	cormax(2,k)=fvec(itmp);

end

% Find the frequency with the maximum correlation coefficient

[max_cc,max_ref]=max(cormax(1,:));
freq_opt=cormax(2,max_ref);
