function smac_setprogress(prog)
% SMAC_SETPROGRESS  Set progress indicator and clear fields in SMAC structure
%
% SMAC_SETPROGRESS(PROG)
%
% SMAC_SETPROGRESS sets certain fields in the SMAC structure and clears
% others, depending on the string supplied in PROG.  It is intended to make
% the structure representative of the location in the curve-fitting process
% as if the user had started from the beginning, in case the user has
% backed up.  This ensures that the structure data is all consistent.

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
%  03-Jun-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for new storage
%
%  09-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for new storage (.nroots,.shape,.corr,.corr_ref)
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Several updates for new data structure members
%
%  31-Aug-2005 / ATA Engineering / Dan Hensley
%    o Update for .shape.refpick
%
%==========================================================================

global ss

% Clear out information depending on where we are in the process

switch prog,
  case 'pinv'

    ss.done.pinv=true;
    ss.done.corr=false;
    ss.done.selroots=false;
    ss.done.autofit=false;
    
    ss.corr.corr=[];
    ss.corr.corrind=[];
    ss.corr.nroots=[];
    ss.corr.corr_ref=[];
    ss.freqrangecc=[];
    
    ss.fit.rootlist=[];
    ss.fit.rootlistorig=[];
    ss.fit.corr_ref=[];

    ss.fa=imat_fn;
    ss.mif.ana.nmif=imat_fn([]);
    ss.mif.ana.cmif=imat_fn([]);
    ss.mif.ana.mmif=imat_fn([]);

    ss.shape.residues=[];
    ss.shape.shape=imat_shp([]);
    ss.shape.refpick=[];
    ss.shape.rootsel=[];

    sstmp=smac_create_data_structure;
    ss.residuals=sstmp.residuals;

  case 'corrcoef'

    ss.done.corr=true;
    ss.done.selroots=false;
    ss.done.autofit=false;
    
%     ss.corr.corrind=[];
%     ss.corr.nroots=[];
    
    ss.fit.rootlist=[];
    ss.fit.rootlistorig=[];
    ss.fit.corr_ref=[];

    ss.fa=imat_fn([]);
    ss.mif.ana.nmif=imat_fn([]);
    ss.mif.ana.cmif=imat_fn([]);
    ss.mif.ana.mmif=imat_fn([]);
    
    ss.shape.residues=[];
    ss.shape.shape=imat_shp([]);
    ss.shape.refpick=[];
    ss.shape.rootsel=[];

    sstmp=smac_create_data_structure;
    ss.residuals=sstmp.residuals;
    
  case 'plotcorr'

    ss.done.corr=true;
    ss.done.selroots=false;
    ss.done.autofit=false;
    
    ss.fit.rootlist=[];
    ss.fit.rootlistorig=[];
    ss.fit.corr_ref=[];

    ss.fa=imat_fn([]);
    ss.mif.ana.nmif=imat_fn([]);
    ss.mif.ana.cmif=imat_fn([]);
    ss.mif.ana.mmif=imat_fn([]);

    ss.shape.residues=[];
    ss.shape.shape=imat_shp([]);
    ss.shape.refpick=[];
    ss.shape.rootsel=[];

  case 'selroots'

    ss.done.selroots=true;
    ss.done.autofit=false;

    ss.fit.rootlist=[];
    ss.fit.rootlistorig=[];
    ss.fit.corr_ref=[];
    
    ss.fa=imat_fn([]);
    ss.mif.ana.nmif=imat_fn([]);
    ss.mif.ana.cmif=imat_fn([]);
    ss.mif.ana.mmif=imat_fn([]);

    ss.shape.residues=[];
    ss.shape.shape=imat_shp([]);
    ss.shape.refpick=[];

  case 'autofit'

    ss.done.autofit=true;
    
    ss.fa=imat_fn([]);
    ss.mif.ana.nmif=imat_fn([]);
    ss.mif.ana.cmif=imat_fn([]);
    ss.mif.ana.mmif=imat_fn([]);

    ss.shape.residues=[];
    ss.shape.shape=imat_shp([]);
    ss.shape.refpick=[];

  otherwise,
    uiwait(errordlg(sprintf('INTERNAL ERROR:  Unknown process step ''%s''',prog)))

end
