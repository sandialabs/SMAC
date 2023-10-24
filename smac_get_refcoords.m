function out=smac_get_refcoords(in,flag)
% SMAC_GET_REFCOORDS  Select reference coordinates from the supplied imat_ctrace
%
% OUT=SMAC_GET_REFCOORDS(IN[,FLAG])
%
% SMAC_GET_REFCOORDS prompts the user for which reference coordinate to use
% from the supplied list of coordinates in the imat_ctrace IN.  FLAG is
% a logical specifying whether only one reference can be selected.  If
% TRUE, SMAC_GET_REFCOORDS will produce an error if more than one reference
% is selected and ask again.  If not specified, FLAG is TRUE.
%
% OUT is an imat_ctrace containing the selected reference.  If the user
% cancels the form, OUT will be numeric.
%
% The user may only select a single reference coordinate at this time.

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
%    o When reading an AFU or UNV, if there are multiple reference prompt
%      the user for which one to use
%    o Start code cleanup and additional comments
%
%==========================================================================

% Argument checking

narginchk(1,2);
if ~isa(in,'imat_ctrace'),
    error('Input IN must be an imat_ctrace of at least one');
end

if ~exist('flag','var'),
    flag = true;
end
flag=logical(flag);
if numel(flag)~=1,
    error('FLAG must be a logical scalar');
end

% Set up some variables prior to asking user

out=-1;
notdone=true;
plural='(s)';
if flag, plural=''; end

% Prompt user for reference(s) to use
while notdone,

    % Prompt the user for the reference to use

    out=uiselect(in,sprintf('Select Reference%s to Use',plural),1:length(in));
    if isnumeric(out), return; end
    
    if flag && length(out)~=1,
        uiwait(errordlg('You can only select one reference'));
    elseif isempty(out),
        uiwait(errordlg('You must select at least one reference'));
    else
        notdone=false;
    end
end
