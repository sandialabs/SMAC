function smac_corr
% SMAC_CORR  Calculate correlation coefficients
%
% Reference:  Mays, R. L., and Johansen, D. D., "A Modal Parameter Extraction
% Algorithm Using Best-Fit Reciprocal Vectors", IMAC 1998.

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
%    o Renamed variables, reorganized, added comments
%    o Cleaned up code
%    o Adapt to new data structure
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Make it work with multiple references
%
%  13-Aug-2004 / ATA Engineering / Dan Hensley
%    o Compute overall correlation coefficient using max instead of mean
%
%  24-Sep-2004 / ATA Engineering / Dan Hensley
%    o Turn off divide by zero warnings
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Fix subscripting for the case where we have just a few lines
%
%==========================================================================

global ss;

% Get the frequency range for the coefficient calculations

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));

% If the first frequency is 0, move up one spectral line

if xf(1)==0,
  xf=xf(2:end);
  ss.freqrangecc(1,:)=[xf(1) 2];
end

xw=2*pi*xf;
bs=length(xf);

% Get the experimental FRF abscissa

wr=2*pi*ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2));
wr=wr(:).';
nel=length(wr);

% Calculate analytical SDOF FRF, assuming A=1, zeta
%    NOTE:  Each row is an FRF with resonant frequency at the spectral line
%           frequency through range of fit

%           This is Hp in equation (5) of the paper

H = -(ones(bs,1)*wr.^2 ./ ...
     (xw.^2*ones(1,nel) + j*2*ss.corr.zeta*xw*wr - ones(bs,1)*wr.^2));

% Loop over each reference

nref=length(ss.ref_coords);
ss.corr.corr_ref=zeros(bs,nref);

for k=1:nref,

  % Compute the reciprocal modal vector Psi and Hp (SDOF FRF)

  if ss.realcomplex==2,	% Complex
    psi = ss.pinv(:,:,k)*H.';
    Hp = ss.fe(:,k).ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:)*psi;
  else			% Real normal
    Hs = [real(H) imag(H)];
    psi = ss.pinv(:,:,k)*Hs';
    Hp = ss.fe(:,k).ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:)*psi;
  end

  % Initialize variables

  Mcc=zeros(bs,1);
  nl=ss.corr.nl;

  % Calculate correlation coefficients for each spectral line in the fit
  % Only look at nl lines around each peak

  ws=warning;
  warning off MATLAB:divideByZero;

  for i=1:bs
    [Y,ind] = max(abs(H(i,:)));
    if ind+nl > nel		% On the right side of the FRF
      R = corrcoef([abs(H(i,max(1,ind-nl):nel)).' abs(Hp(max(1,ind-nl):nel,i))]);
    elseif ind-nl <= 0		% On the left side of the FRF
      R = corrcoef([abs(H(i,1:min(ind+nl,nel))).' abs(Hp(1:min(ind+nl,nel),i))]);
    else				% In the middle
      R = corrcoef([abs(H(i,max(1,ind-nl):min(ind+nl,nel))).' abs(Hp(max(1,ind-nl):min(ind+nl,nel),i))]);
    end
  
    Mcc(i,1) = R(1,2);
  end
  warning(ws);

  % Store

  ss.corr.corr_ref(:,k)=Mcc;

end

% Average the individual correlation coefficient vectors

%ss.corr.corr=mean(ss.corr.corr_ref,2);
[ss.corr.corr(:,1) ss.corr.corr(:,2)]=max(ss.corr.corr_ref,[],2);
