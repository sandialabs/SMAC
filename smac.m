function ss=smac(input)
%              SMAC Modal Extraction Package
%
%  This code allows the user to load the frf matrix and initiate
%  the Graphical User Interface for SMAC

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
%  03-Jun-2004 / ATA Engineering / Dan Hensley
%    o Make sure the code handles partially finished processes (pick up
%      where you left off)
%    o Clean out remaining old code
%
%  07-Jun-2004 / ATA Engineering / Dan Hensley
%    o Display splash screen
%
%  28-Jun-2004 / ATA Engineering / Dan Hensley
%    o Clear out remaining code using the old structure
%
%  09-Jul-2004 / ATA Engineering / Dan Hensley
%    o Make sure at least one FRF is selected
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Code changes to handle multiple references
%    o Turn into a function so we can better organize things
%    o Allow input argument to SMAC
%
%  09-Aug-2004 / ATA Engineering / Dan Hensley
%    o Work around IMAT v1.31 not having @ideas_ctrace/intersect
%
%  13-Aug-2004 / ATA Engineering / Dan Hensley
%    o Fix minor bug that caused it to fail for the single-reference case
%
%  02-Nov-2004 / ATA Engineering / Dan Hensley
%    o Make sure we pick up where we left off with the new GUI system
%
%  30-Nov-2004 / ATA Engineering / Dan Hensley
%    o Minor fix to get around case where there are multiple functions with the
%    same reference/response pair
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Add check for Matlab version 7 (so we can save .mat with -v6)
%
%  24-Feb-2005 / ATA Engineering / Dan Hensley
%    o Preselect all FRF when presenting to user
%
%  02-Jun-2005 / ATA Engineering / Dan Hensley
%    o Only apply IMAT v131 patch to version 1.31 (add path)
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o If loading a .mat file, make sure the structure is the latest version
%
%==========================================================================

imat_fn;    % Initialize IMAT

%setunits('IN');
global ss;
global SMAC_SAVE_V6

% Make sure the IMAT pre-requisites are here

if isempty(which('imat_fn'))
    error('SMAC requires IMAT to run');
end
if isempty(which('nmif')) || isempty(which('mmif')) || isempty(which('cmif')) || isempty(which('ortho'))
    error('You need the IMAT ''examples'' directory in your MATLAB path (for cmif, mmif, nmif, and ortho)');
end

% Initialize some variables

g=imat_fn([]);
ssin=false;

SMAC_SAVE_V6=false;     % Save version 6 MAT file?

% Post copyright notice

fprintf('SMAC  Copyright (C) 2008  Sandia Corporation\n');
fprintf('This program comes with ABSOLUTELY NO WARRANTY.\n');
fprintf('This is free software, and is licensed under the Mozilla Public License (MPL) version 1.1.\n');
fprintf('You are welcome to redistribute it under certain conditions.\n\n');

% Initialize SMAC structure

ss=smac_create_data_structure;

% Display the splash screen

hsplash=smac_splash;

% Start at the beginning of the curve-fitting process

smacproc=1;

% FIXME:  If variable exists in workspace, see if we want to start from
% there

%==========================================================================
% Determine which type of input file will be used in SMAC

if ~exist('input','var'),
    In_buttonName=questdlg('Enter the type of data file for your SMAC Analysis',...
        'Welcome to Sythesize Modes and Correlate (SMAC)', ...
        '*.afu','*.unv','*.mat','*.afu');
    if isempty(In_buttonName),
        if ishandle(hsplash), delete(hsplash); end
        return;
    end
else

    % Look at the input arguments to decide what to do

    if     isa(input,'imat_fn'),
        g=input;
        In_buttonName='*.afu';
    elseif isstruct(input),	% FIXME:  Make sure the structure is a valid ss structure
        ss=input;
        ssin=true;
        In_buttonName='*.mat';
    end
end

%==========================================================================
switch In_buttonName,

    %------------------------------------------------------------------------
    case '*.afu',

        % Use browser to load IDEAS FRF Data file into SMAC

        if isempty(g),
            [fname,pname] = uigetfile('*.afu',...
                'Load IDEAS FRF Data file');
            if isnumeric(fname),
                if ishandle(hsplash), delete(hsplash); end
                return;
            end
            ss.filename=fullfile(pname,fname);

            g = readadf(ss.filename);
        end

        % Process FRFs

        [f,refs]=process_functions(g);

        if isnumeric(f),
            if ishandle(hsplash), delete(hsplash); end
            return;
        end

        % Put the relevant information into the structure

        ss.fe=f;
        ss.freqrange=[ss.fe(1).abscissa([1 end]) [1;ss.fe(1).numberelements]];
        ss.ref_coords=refs;


        %------------------------------------------------------------------------
    case '*.unv',

        % Use browser to load Universal FRF Data file into SMAC

        [fname,pname] = uigetfile('*.unv',...
            'Load Universal FRF Data file');
        if isnumeric(fname),
            if ishandle(hsplash), delete(hsplash); end
            return;
        end
        ss.filename=fullfile(pname,fname);

        g = readunv(ss.filename,[58 1858]);

        % Process FRFs

        [f,refs]=process_functions(g);

        if isnumeric(f),
            if ishandle(hsplash), delete(hsplash); end
            return;
        end

        % Put the relevant information into the structure

        ss.fe=f;
        ss.freqrange=[ss.fe(1).abscissa([1 end]) [1;ss.fe(1).numberelements]];
        ss.ref_coords=refs;


        %------------------------------------------------------------------------
    case '*.mat',

        % Use browser to load MATLAB FRF Data file into SMAC

        if ~ssin,
            [fname,pname] = uigetfile('*.mat','Select MATLAB .mat file containing SMAC structure');
            if isnumeric(fname),
                if ishandle(hsplash), delete(hsplash); end
                return;
            end
            fname=fullfile(pname,fname);

            clear ss;
            load(fname);

            if ~exist('ss','var'),
                uiwait(errordlg(sprintf('File ''%s'' does not contain a SMAC structure',fname)));
                return
            elseif ~isstruct(ss),
                uiwait(errordlg('SMAC structure variable is not a structure'));
            end
            
            % Make sure the structure is the latest version
            
            ss=smac_update_data_structure(ss);
            if isempty(ss)
                error('SMAC structure is invalid or could not be updated to the latest version');
            end            

        end

        % Figure out where to start in the SMAC process

        if ss.done.pinv, smacproc=2; end
        if ss.done.corr, smacproc=3; end
        if ss.done.selroots, smacproc=4; end
        if ss.done.autofit, smacproc=5; end

end

%==========================================================================
% Now start the curve-fitting process

switch smacproc,

    case 1,	% Pseudo-inverse
        smac_GUI_pinv;

    case 2,	% Get correlation coefficients
        smac_GUI_corrcoef;

    case 3,	% Get initial frequency list
        smac_GUI_plotcorr;

    case 4,	% Perform initial auto-fit on the selected peaks
        smac_GUI_autofit;

    case 5,	% Synthesize FRF and final root maintenance
        smac_GUI_synth;

    otherwise,
        uiwait(errordlg('Unknown step in the process'));
        return;

end

% Display some final information

disp(' ')
fprintf('	Low frequency of experimental data:  %g Hz\n',ss.fe(1).abscissa(1));
fprintf('      High frequency of experimental data:  %g Hz\n',ss.fe(1).abscissa(end));
fprintf('   Experimental data frequency resolution:  %g Hz\n',ss.fe(1).abscissainc);
disp(' ')

if ishandle(hsplash), delete(hsplash); end


%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================
function [f,refs]=process_functions(g)
% Determine which FRF to use

f=-1;

% Prompt for which reference to use if multiple references found

filt=imat_filt({'FunctionType','==','Frequency Response Function'; ...
                 'ResponseNode','~=',9999                           });

g=g{filt};
refs=unique(allref(g));
if length(refs)>1,
    %refs=smac_get_refcoords(unique(allref(g)),logical(1));
    refs=smac_get_refcoords(unique(allref(g)),false);
    if isnumeric(refs), return; end
end
refsc=cellstr(refs);

filt=imat_filt('ReferenceCoord','=',refs);

% Display filter message

disp('smac.m: Setting Filters...');
disp('FunctionType == Frequency Response Function');
disp('ResponseNode ~= 9999');
fprintf('ReferenceCoord == %s\n',refsc{:});

notdone=true;
while notdone,

    f=imat_fn([]);
    while isempty(f),
        f=uiselect(g{filt},1:length(g{filt}));
        if isnumeric(f), return; end

        if length(unique(allref(f)))~=length(refs), f=imat_fn([]); end

        if isempty(f),
            uiwait(errordlg('You must select at least one FRF from each reference!'));
        end
    end

    % Make sure we have the same number of responses for each reference

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
                %f(:,k)=ff{'ResponseCoord','=',res};
                f(:,k)=ff{res};
            end

            if prod(size(f))~=prod(size(g2)),
                fprintf('** Warning:  Some FRF removed to have a consistent set of responses for each reference\n');
            end
            notdone=false;
        else
            uiwait(errordlg('There must be at least one common response from each reference'));

        end
    else		% Single reference
        notdone=false;

    end

end

%==========================================================================
