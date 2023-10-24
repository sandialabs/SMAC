function ss=smac_resynth
% SMAC_RESYNTH  Front end for resynthesizing FRF from mode shapes
%
% SMAC_RESYNTH will first prompt the user for a shape file (ASH) from which
% to work.  It will then prompt the user for an experimental FRF file (AFU)
% containing the functions to which to compare.  SMAC_RESYNTH will parse
% the FRF down to the common set specified in the shape file.

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
%  17-Jun-2005 / ATA Engineering / Dan Hensley
%    o Initial release
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Update for new low and high mode information in IDLine3
%    o Fix bug that caused a crash if some of the residual shapes were not
%      there (because they weren't used)
%
%==========================================================================

global ss;
global SMAC_SAVE_V6

% Put IMAT patch into the directory tree only if we need to

try
    [v,n]=imat_ver;
    if strcmp(n,'1.31')
        pp=path;
        ps=mfilename('fullpath');
        ps=ps(1:end-length(mfilename)-1);
        np=fullfile(ps,'imat_v131_patch');
        if isempty(findstr(np,pp))
            addpath(np);
        end
    end
    if strcmp(n,'1.41') || strcmp(n,'1.42')
        pp=path;
        ps=mfilename('fullpath');
        ps=ps(1:end-length(mfilename)-1);
        np=fullfile(ps,'imat_v142_patch');
        if isempty(findstr(np,pp))
            addpath(np);
        end
    end
catch
    error('You need IMAT v1.31 or higher to run this version of SMAC');
end

% Initialize some variables

g=imat_fn([]);
ssin=false;

vv=ver('Matlab');
if vv.Version(1)=='7',
  SMAC_SAVE_V6=true;
else
  SMAC_SAVE_V6=false;
end

% Initialize SMAC structure

ss=smac_create_data_structure;

% Display the splash screen

hsplash=smac_splash;


%--------------------------------------------------------------------------
% Read the shape file and store the shapes

fname = imat_getfile('*.ash','Load I-deas Shape Data file');
if isnumeric(fname),
    if ishandle(hsplash), delete(hsplash); end
    return;
end
shp = readadf(fname);
ss.shape.shape = shp(shp.damping~=-1);

% Read the function file

fname = imat_getfile('*.afu','Load I-deas FRF Data file');
if isnumeric(fname),
    if ishandle(hsplash), delete(hsplash);, end
    return;
end
frf = readadf(fname);

% Filter to keep only FRF 
filt=imat_filt({'FunctionType','==','Frequency Response Function'; ...
                 'ResponseNode','~=',9999                           });
frf=frf{filt};
if isempty(frf)
    error('No FRF to process!');
end

% Keep only the FRF that have the DOF in the shape

doffrf = abs(allres(frf(:)));
dofshp = abs(alldof(ss.shape.shape(:)));

ind = ismember(id(doffrf),id(dofshp),'rows');
frf = frf(ind);
if isempty(frf)
    error('FRF and shapes have no DOF in common!');
end

% Rearrange the functions by reference and common responses

[frf,refs] = process_functions(frf);

% Look for residual information in the shape file

res=shp(shp.damping==-1);

if ~isempty(res)
    
    res = res(strcmpi(res.idline5,'smac residual information'));
    if isempty(res)
        fprintf('*** Shape file has shapes with damping = -100% but they are not residuals\n');
    end
    
    % Extract the residual information from the shapes
    
    [refc,resc,residuals,frf] = process_shape_residuals(res,frf);
    ss.residuals = residuals;
else
    refc=unique(allref(frf));
end

% Store the data into the SMAC data structure

ss.ref_coords = refc;
ss.fe = frf;

if any(strcmpi(ss.shape.shape.shapetype,'complex'))
    ss.realcomplex = 2;
else
    ss.realcomplex = 1;
end

% Now start the resynthesis

if ishandle(hsplash), delete(hsplash);, end
smac_GUI_resynth;


%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================

function [f,refs]=process_functions(f)
% Determine which FRF to use

refs=unique(allref(f));
refsc=cellstr(refs);

% Get the common set of references and responses

if length(refs)>1,
    g2=f;

    ind=false(length(f),length(refs));
    for k=1:length(refs),
        ind(:,k)=in_filter(imat_filt('ReferenceCoord','=',refs(k)),f);
    end

    % Get the responses for each reference

    res=allres(f(ind(:,1)));
    for k=2:length(refs),
        ar=allres(f(ind(:,k)));
        [tmp,rind]=intersect(res.id,ar.id,'rows');
        res=res(rind);
    end

    % Get the common set of response coordinates

    if ~isempty(res),

        ftmp=f;
        f=imat_fn(length(res),length(refs));
        for k=1:length(refs),
            ff=ftmp(ind(:,k));
            f(:,k)=ff{res};
        end

        if prod(size(f))~=prod(size(g2)),
            fprintf('** Warning:  Some FRF removed to have a consistent set of responses for each reference\n');
        end
        notdone=false;
    else
        uiwait(errordlg('There must be at least one common response from each reference'));

    end
end


%==========================================================================
function [refc,resc,residuals,frf]=process_shape_residuals(shp,frf)
% Extract the residual information from the shapes

tmp=smac_create_data_structure;
residuals=tmp.residuals;
use_res=false;

% Look for the index vector

shpvec=shp(shp.frequency==1000);
if length(shpvec)>1
    error('Too many DOF indexing shapes');
end
idline1=shpvec.idline1;
if ~iscell(idline1), idline1={idline1}; end
if isempty(shpvec) | isempty(findstr(idline1{:},'DOF indexing shape'))
    error('Could not find DOF indexing shape');
end

% Get the response coordinates in order

resc=alldof(shpvec);
ii=shpvec.shape;
resc(ii<0)=-resc(ii<0);
resc=resc(ii~=0);
[tmp,ind]=sort(abs(ii(ii~=0)));
resc=resc(ind);

% Now get the reference coordinates from IDLine2

str=shpvec.idline2;
if ~iscell(str), str={str}; end
tmp=strread(str{:},'%s','delimiter',' ');
if length(tmp)<2
    error('Invalid IDLine2 in DOF indexing shape (expecting reference list)');
end
refc=imat_ctrace(tmp(2:end));

%--------------------------------------------------------------------------
% Throw out the FRF whose response coordinates don't match the residuals

[ind,loc]=ismember(id(allres(frf(:,1))),id(resc),'rows');
if isempty(frf)
    error('No FRF response coordinates match the residual coordinates!');
end
if ~all(ind)
    fprintf('** Throwing out %d FRF whose response coordinates do not match the residual coordinates\n', ...
            length(find(~ind))*size(frf,2));
end
frf=frf(ind,:);

% Now sort the responses to match the residual order
% NOTE:  indres says which residues to keep

[indres,loc]=ismember(id(resc),id(allres(frf(:,1))),'rows');
frf=frf(loc(loc~=0),:);

if ~all(indres)
    fprintf('** Throwing out %d residual terms whose DOF are not present in the FRF\n',length(find(~indres))*length(refc))
end

%-------------------------------
% Rearrange FRF according to matching reference coordinates

[ind,loc]=ismember(id(allref(frf(1,:)')),id(refc),'rows');
frf=frf(:,loc(loc~=0));
if isempty(frf)
    error('No FRF reference coordinates match the residual coordinates!');
end
if ~all(ind)
    fprintf('** Throwing out %d FRF whose reference coordinates do not match the residual reference(s)\n', ...
            length(find(~ind))*size(frf,1));
end

% Throw out residual shapes whose references don't match what are in the FRF

[indref,locref]=ismember(id(refc),id(allref(frf(1,:)')),'rows');
refc=refc(locref(locref~=0));

% Build residual keep index

nref=length(refc);
nres=size(frf,1);
indres=find(indres(:));
indres=repmat(indres,1,nref)+repmat(0:nres:(nref*(nres-1)),length(indres),1);
indres=indres(:).';

%--------------------------------------------------------------------------
% Process each of the supplied residual shapes

% : Low and High Mode :

shplh={'lowmode',  100, 200; ...
       'highmode', 200, 300 };
   
for k=1:size(shplh,1)
    s2=shp(shp.frequency>shplh{k,2} & shp.frequency<shplh{k,3});
    if ~isempty(s2)
        s2=s2(locref(locref~=0));
        if length(s2)~=length(refc)
            error('Inconsistency between # of references and # of %s residual shapes',shplh{k,1});
        end
        
        % Sort the shapes by increasing frequency

        [tmp,ind]=sort(s2.frequency);
        s2=s2(ind);

        idline4=s2(1).idline4;
        if ~iscell(idline4), idline4={idline4}; end
        val=str2num(idline4{:});
        if length(val)~=2
            error('%s residual shape IDLine4 does not contain frequency and damping',shplh{k,1});
        end

        % Extract frequency range
        
        idline3=s2(1).idline3;
        if ~iscell(idline3), idline3={idline3}; end
        frange=str2num(idline3{:});
        if length(frange)~=2
            error('%s residual shape IDLine3 does not contain the frequency range',shplh{k,1});
        end


        % Extract the frequency and damping and residual values

        residuals.(shplh{k,1}).active=true;
        residuals.(shplh{k,1}).freq=val(1);
        residuals.(shplh{k,1}).damp=val(2);
        residuals.(shplh{k,1}).frange=frange;
        residuals.(shplh{k,1}).residual=reshape(s2{resc},1,[]);
        residuals.(shplh{k,1}).residual=residuals.(shplh{k,1}).residual(indres);       
        use_res=true;
    end
end

% % : Mode Inertance and Compliance :
% 
% shpmic={'modeinertance',  300, 400; ...
%         'modecompliance', 400, 500 };
%    
% for k=1:size(shpmic,1)
%     s2=shp(shp.frequency>shpmic{k,2} & shp.frequency<shpmic{k,3});
%     if ~isempty(s2)
%         s2=s2(locref(locref~=0));
%         if length(s2)~=length(refc)
%             error('Inconsistency between # of references and # of %s residual shapes',shpmic{k,1});
%         end
%         
%         % Sort the shapes by increasing frequency
%         
%         [tmp,ind]=sort(s2.frequency);
%         s2=s2(ind);
%         
%         % Extract the residual values
%         
%         residuals.(shpmic{k,1}).active=true;
%         residuals.(shpmic{k,1}).residual=reshape(s2{resc},1,[]);
%         residuals.(shpmic{k,1}).residual=residuals.(shpmic{k,1}).residual(indres);       
%         use_res=true;
%     end
% end

% : General Inertance and Compliance :

shpic={'inertance',  500, 600; ...
       'compliance', 600, 700 };
   
for k=1:size(shpic,1)
    s2=shp(shp.frequency>shpic{k,2} & shp.frequency<shpic{k,3});
    if ~isempty(s2)
        s2=s2(locref(locref~=0));
        if length(s2)~=length(refc)
            error('Inconsistency between # of references and # of %s residual shapes',shpic{k,1});
        end
        
        % Sort the shapes by increasing frequency
        
        [tmp,ind]=sort(s2.frequency);
        s2=s2(ind);
        
        % Extract the frequency range
       
        idline4=s2(1).idline4;
        if ~iscell(idline4), idline4={idline4}; end
        val=str2num(idline4{:});
        if length(val)~=2
            error('%s residual shape IDLine4 does not contain frequency range',shpic{k,1});
        end
        residuals.(shpic{k,1}).freq=val(:).';
        
        % Extract the residual values
        
        residuals.(shpic{k,1}).active=true;
        residuals.(shpic{k,1}).residual=reshape(s2{resc},1,[]);
        residuals.(shpic{k,1}).residual=residuals.(shpic{k,1}).residual(indres);       
        use_res=true;
    end
end

% Set residual use flag

residuals.use=use_res;

%==========================================================================
