function smac_list_ref_info
% SMAC_LIST_REF_INFO  List CC and shape selection reference information
%
% SMAC_LIST_REF_INFO lists a table of information regarding which reference
% generated the highest correlation coefficient and which reference was
% ultimately used to select the shape.
%
% The first column shows the root number.  This matches the roon number
% shown in the Synthesis form.  The second and third columns are the
% frequency and damping for each root, respectively.  The fourth column
% shows the final correlation coefficient value.  The fifth columns lists
% the reference coordinate that produced that correlation coefficient.  If
% a shape was generated for that root, the sixth column lists the reference
% that generated this shape.

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

global ss

% List the table data

colwidth=[5 9 7 6 10 10];

% Display the table title

fmtttl=sprintf('%%-%ds  %%-%ds  %%-%ds  %%-%ds    %%-%ds  %%-%ds\n',colwidth);

fprintf('\n');
fprintf(fmtttl,'Root','Frequency','Damping','','Max CC','Shape');
fprintf(fmtttl,'Num','(Hz)','(%cr)','Max CC','Reference','Reference');
dashes=cell(length(colwidth),1);
for k=1:length(dashes)
	dashes{k}=repmat('-',1,colwidth(k));
end
fprintf(fmtttl,dashes{:});

% Display the table data

fmtdata=sprintf('%%%dd  %%%d.3f  %%%d.3f  %%%d.3f    %%-%ds  %%-%ds\n',colwidth);

for k=1:size(ss.fit.rootlist,1)
	data=[k ss.fit.rootlist(k,1:3).*[1 100 1]];
	
	ccind=ss.fit.rootlist(k,4);
	ccref='';
	if ismember(ccind,1:length(ss.ref_coords))
		ccref=char(ss.ref_coords(ccind));
	end
	
	shpref='';
	ind=find(k==ss.shape.rootsel);
	if ~isempty(ind) && ~isempty(ss.shape.refpick)
		shpref=char(ss.ref_coords(ss.shape.refpick(ind)));
	end
	fprintf(fmtdata,data,ccref,shpref);
end

fprintf('\n');
