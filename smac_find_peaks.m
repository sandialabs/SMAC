function [pv,pind]=smac_find_peaks(fn,thresh)
% SMAC_FIND_PEAKS  Find peaks in the supplied function
%
% [PV,PIND]=SMAC_FIND_PEAKS(FN[,THRESH])
%
% SMAC_FIND_PEAKS searches for peaks in the supplied numeric vector FN.
% THRESH is an optional input specifying a value in FN below which to
% ignore peaks.  SMAC_FIND_PEAKS looks for points where the sign of the
% slope of the function changes from positive to negative.
%
% PV is a vector containing the peak values found.  PIND is a vector the
% same size as PV containing the indices into FN where the peak values
% were found.

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
%  22-May-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%==========================================================================

% Input argument checking

narginchk(1,2)

% Find where the slope changes from positive to negative

pind=find(diff(sign(diff(fn)))==-2)+1;

% Keep the peaks above the threshhold

if exist('thresh','var'),
  pind=pind(fn(pind)>thresh);
end

pv=fn(pind);
