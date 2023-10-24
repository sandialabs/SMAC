function [shape,refpick,psi]=smac_create_shapes(rootlist,residuals,fittype,diagtoggle)
% SMAC_CREATE_SHAPES  Generate mode shapes for the selected roots
%
% [SHP,REFPICK,PSI]=SMAC_CREATE_SHAPES(ROOTLIST,RESIDUALS,FITTYPE[,DIAGTOGGLE])
%
% SMAC_CREATE_SHAPES generates mode shapes for the roots supplied in
% ROOTLIST.  ROOTLIST is an nx3 matrix, where column 2 contains the
% frequencies and column 3 contains the damping values.  RESIDUALS
% is a structure containing information for the various types of residuals
% to include in the calculations.  FITTYPE specifies the type of shape
% fitting to use for multi-reference curve fits.  In all cases the mode
% shapes are fit from all references.  The final modes are selected from
% one of three methods.
% Choices are
%    FITTYPE    SHAPE SELECTION ALGORITHM
%     1       - Select the shape whose drive point magnitude is highest
%     2       - Select the shape whose correlation coefficient closest
%               to the mode frequency is highest
%     3 (def) - Select the shape whose synthesis most closely matches the
%               analytical FRF synthesized from all modes selected
%
% DIAGTOGGLE is a logical specifying whether diagnostic information such as
% the fit tables and graphs are displayed when generating the shapes.  The
% default is FALSE.
%
% In all cases, SMAC_CREATE_SHAPES will select shapes using all three
% methods and will display summaries from each, but it will only return
% shapes from the method specified.
%
% SHP is an imat_shp containing the mode shapes for all of the roots
% provided.  REFPICK is a vector specifiying the reference that created
% each shape that was selected and returned in SHP.  PSI is the reciprocal
% modal vectors calculated from the shape parameters

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
%  02-Jun-2004 / ATA Engineering / Dan Hensley
%    o Initial creation; adapted from shaper.m and shaper1.m
%
%  28-Jun-2004 / ATA Engineering / Dan Hensley
%    o Just keep the first drive point index if more than one match is
%      found
%
%  27-Jul-2004 / ATA Engineering / Dan Hensley
%    o Fix bug that caused an error if a drive point match was not found
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for multiple references--use reference with highest amplitude
%      drive point for shape selection
%
%  02-Nov-2004 / ATA Engineering / Dan Hensley
%    o Select shape based on synthesized FRF that most closely matches
%      analytical drive point FRF (correlation coefficient)
%    o Make sure drive point selection works when reference and response
%      coordinates have different signs
%
%  08-Nov-2004 / ATA Engineering / Dan Hensley
%    o Add shape selection method by reference with max CC at that
%      frequency
%    o Allow user to select one of three shape fitting methods
%    o Update FRF synthesis plot legend
%    o Minor math updates to synthesis
%    o Don't take absolute value of drive point shape coefficient unless
%      it's a real normal mode.  If negative, print a warning message.
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Add mask for negative drive point coefficient mask and don't
%      consider mode shapes that have negative coefficients (real normal
%      modes only)
%
%  17-Feb-2005 / ATA Engineering / Dan Hensley
%    o When selecting shapes via max CC, use max value and not max absolute
%      value
%
%  22-Feb-2005 / ATA Engineering / Dan Hensley
%    o Fix bug when selecting too many DP FRF for a given reference
%    o Support toggle for displaying diagnostic information
%
%  24-Feb-2005 / ATA Engineering / Dan Hensley
%    o Add some drawnow's to force the plots to be displayed (Matlab 7)
%
%  07-Mar-2005 / ATA Engineering / Dan Hensley
%    o Multiply synthesized FRF by dpsign so phase is correct on Best Synth
%      diagnostic plots
%
%  11-Mar-2005 / ATA Engineering / Dan Hensley
%    o When doing a Best Synth fit, go out to the half power points for the
%      CC comparison
%
%  03-Jun-2005 / ATA Engineering / Dan Hensley
%    o Fix bug where drive point shape coefficient could be stored as
%      negative even when the residue was positive
%    o Add support for various residual types
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Support new call to smac_synthesize (add 3rd output)
%
%  14-Sep-2005 / ATA Engineering / Dan Hensley
%    o Greatly improve output tables, highlighting negative DP entries and
%      which entry was selected
%    o Fix bug in Max CC selection where it didn't index in to the roots
%      actually selected for processing
%    o Trap and warn user when DP residue is negative for all references
%      for a given mode
%
%  06-Mar-2006 / ATA Engineering / Dan Hensley
%    o Calculate psi vectors when done creating the shapes
%
%==========================================================================

global ss;
shape=imat_shp([]);
deffittype=3;
refpick=[];
psi=[];
dpnegwarn=false;

% Argument checking

narginchk(2,4);
if nargin==2
    fittype=deffittype;
end
if ~exist('diagtoggle','var')
    diagtoggle=false;
end

% See if we can determine the drive point FRF

dpind=[];
refcoord=ss.ref_coords;		% Get reference coordinate(s)
nres=size(ss.fe,1);
nref=length(refcoord);

ctref=cell(nref,1);
ctres=cell(nref,1);

for k=1:nref
    ctref{k}=imat_ctrace(ss.fe(:,k).referencecoord);
    ctres{k}=imat_ctrace(ss.fe(:,k).responsecoord);
    dbref(:,:,k)=double(ctref{k});
    dbres(:,:,k)=double(ctres{k});

    refcoorda=abs(ctref{k});
    rescoorda=abs(ctres{k});

    %ind=strmatch(refcoord{k},char(rescoorda));	% Single ref
    ind=find(rescoorda==abs(refcoord(k)));	% Multi ref
    if isempty(ind), ind=0; end
    dpind(k)=ind(1);
end

%% Get the drive point FRF index for each reference

for k=1:nref
    inddp=[];
    while length(inddp)~=1
        selind=dpind(k);
        if selind<1, selind=[]; end
        [tmp,inddp]=uiselect(ss.fe(:,k),sprintf('Select the drive point for reference %d (%s)', ...
            k,char(refcoord(k))),selind);
        if isnumeric(tmp), return; end
        if length(inddp)~=1
            uiwait(errordlg('You may only select a single FRF as the drive point'));
        else
            dpind(k)=inddp;
        end
    end
end
drawnow;

% Synthesize the FRF

if isempty(ss.fa)
    [frfsyn,residues_all,residuals]=smac_synthesize(rootlist,residuals);
    ss.fa=frfsyn;
    ss.shape.residues=residues_all;
    if ss.residuals.use
        ss.residuals=residuals;
    end
else
    frfsyn=ss.fa;
    residues_all=ss.shape.residues;
    residuals=ss.residuals;
end

%--------------------------------------------------------------------------
%% Loop over each reference

numroots=size(rootlist,1);
shape=imat_shp(numroots,nref);
dp=zeros(numroots,nref);
dpsign=zeros(nref,1);
negdpmask=false(numroots,nref);

for k=1:nref

    % Select driving point for this reference

    inddp=dpind(k);

    % Multiply residues by the sign of driving point

    irind=nres*(k-1);
    residues=residues_all(:,irind+1:irind+nres);
    dpsign(k)=sign(dbres(inddp,2,k)*dbref(inddp,2,k));	% Sign of drive point response
    res=residues;

    % If the drive point response coordinate is negative and the sign is
    % negative, flip the coefficients on the whole shape
    
    if dpsign(k)<0 && dbres(inddp,2,k) < 0
        ctres{k}(inddp)=abs(ctres{k}(inddp));
        res(:,inddp)=-res(:,inddp);
    end
    
    % When checking coefficients, a negative reference will produce valid
    % negative residues that should not be trapped below
    
    if dpsign(k)<0 && dbref(inddp,2,k)<0
        negfact=-1;
    else
        negfact=1;
    end

    res=res.';
    res=res(:,1:numroots);

    %  Drive point shape coefficients

    if ss.realcomplex==1
        %ii=find(negfact*dpsign(k)*res(inddp,:)<0); % 3-31-11 error we think
        ii=find(negfact*res(inddp,:)<0);
        negdpmask(ii,k)=true;
        if ~isempty(ii)
            fprintf('*** Drive point residue for mode %d reference %d is negative\n',[ii(:) k*ones(size(ii(:)))]');
        end
        phii=sqrt(abs(res(inddp,:)));
    else
%         ii=find(negfact*angle(res(inddp,:))>0);
%         negdpmask(ii,k)=true;
%         if ~isempty(ii)
%             fprintf('*** Drive point residue for mode %d reference %d is negative\n',[ii(:) k*ones(size(ii(:)))]');
%         end
        phii=sqrt(res(inddp,:));
    end

    % Make sure the drive point shape coefficient is non-zero

    if any(phii==0)
        ind=find(phii==0);
        rootnum=ss.shape.rootsel(ind);
        str=sprintf( [...
            'Drive point residue for reference %d (%s) is 0.  This indicates a severe\n' ...
            'problem with the synthesis.  You may have too many roots for the number\n' ...
            'of spectral lines available.\n\n'],k,char(ctref{k}(1)));
        plural='';
        if length(rootnum)>1, plural='s'; end
        str=[str sprintf('This affects root%s %d',plural,rootnum(1))];
        if length(rootnum)>1
            str=[str sprintf(', %d',rootnum(2:end-1)) sprintf(' and %d.',rootnum(end))];
        else
            str=[str '.'];
        end
        str=[str sprintf('\n\nPlease correct this before continuing.')];
        uiwait(errordlg(str));
        return;
    end

    % Compute the mode shapes by scaling residues with shape coefficients

    phij=res./repmat(phii,length(ctres{k}),1);
    dp(:,k)=phij(inddp,:).';

    % Store all of this in an imat_shp

    shape(:,k) = build_shape(ctres{k},phij);

    % Set whether this is a real or complex shape

    if ss.realcomplex== 2
        shape.shapetype = 'Complex';
    else
        shape.shapetype = 'Real';
    end

    % Set some attributes

    shape(:,k).frequency=rootlist(:,2);
    shape(:,k).damping=rootlist(:,3);
    shape(:,k).modalmassreal=1;
    shape(:,k).referencecoord=char(refcoord(k));
    shape(:,k).responsecoord=char(ctres{k}(inddp));
    shape(:,k).orddendatatype=ss.fe(inddp).orddendatatype;
    shape(:,k).ordnumdatatype=ss.fe(inddp).ordnumdatatype;
    shape(:,k).idline1=sprintf('Generated from reference %s',char(refcoord(k)));

    if diagtoggle
        figure('NumberTitle','off', ...
            'Name',sprintf('Self-MAC Mode Shapes from reference %d (%s)',k,char(refcoord(k))));
        bar3o(ortho(shape(:,k)));
        drawnow;
    end

end

refpicktmp=ones(numroots,1);

% If we have one reference, we have the shapes and references we want

if nref==1
    shpkeep=shape;
    refpick=refpicktmp;
else

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Method 1: Max Drive Point coefficient

    % We pick the shape with the highest amplitude drive point shape coefficient

    shpkeeptmp=imat_shp(numroots);
    dptmp=dp;
    dptmp(negdpmask)=-inf;	% Don't consider shapes with negative drive point coefficients

    for k=1:numroots
        [~,ind]=max(dptmp(k,:));
        shpkeeptmp(k)=shape(k,ind);
        refpicktmp(k)=ind;
    end

    % Display output:

    if fittype==1 || diagtoggle
        fprintf('\n');
        fprintf('Shapes selected by drive point amplitude (DP listed per reference)\n');
        fprintf('   NOTE:  DP''s surrounded by # indicate negative drive point residue, and\n');
        fprintf('          DP''s with a trailing * incidate the one for which the shape was kept.\n');
        fprintf('\n');

        cref=cellstr(refcoord);
        under=repmat({'-------'},1,nref+3);
        refchar=strjust(char(refcoord),'right');
        fmtref=['%' sprintf('%ds',size(refchar,2)+4)];

        fprintf([repmat('%11s',1,nref+2) '%' sprintf('%ds',size(refchar,2)+4) '\n'], ...
            'Mode #','Freq (Hz)',cref{:},'Ref',under{:});
        for k=1:numroots
            str=sprintf('%11d%11.3f',rootlist(k,1),shpkeeptmp(k).frequency);
            for m=1:nref
                if ~isinf(dptmp(k,m))
                    str=[str sprintf('%10.3f',dptmp(k,m))];
                else
                    str=[str sprintf('%10s',sprintf('#%.3f#',dp(k,m)))];
                end
                if m==refpicktmp(k)
                    str=[str '*'];
                else
                    str=[str ' '];
                end
            end
            str=[str sprintf(fmtref,char(refcoord(refpicktmp(k))))];
            fprintf('%s\n',str);
        end
    end

    % Store shape if desired

    if fittype==1
        fprintf('** Selecting shapes based on maximum drive point amplitude\n');
        shpkeep=shpkeeptmp;
        refpick=refpicktmp;
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Method 2: Highest correlation coefficient

    % We pick the shape with the highest correlation coefficient

    shpkeeptmp=imat_shp(numroots);
    cc=ss.fit.corr_ref(rootlist(:,1),:);
    cctmp=cc;
    cctmp(negdpmask)=-inf;

    frall=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2))';
    for k=1:numroots
        %[tmp,ind2]=max(abs(cctmp(k,:)));
        [~,ind2]=max(cctmp(k,:));

        shpkeeptmp(k)=shape(k,ind2);
        refpicktmp(k)=ind2;
    end

    % Display output:

    if fittype==2 || diagtoggle
        fprintf('\n');
        fprintf('Shapes selected by max correlation coefficient (CC listed per reference)\n');
        fprintf('   NOTE:  CC''s surrounded by # indicate negative drive point residue, and\n');
        fprintf('          CC''s with a trailing * incidate the one for which the shape was kept.\n');
        fprintf('\n');

        cref=cellstr(refcoord);
        under=repmat({'-------'},1,nref+3);
        refchar=strjust(char(refcoord),'right');
        fmtref=['%' sprintf('%ds',size(refchar,2)+4)];

        fprintf([repmat('%11s',1,nref+2) '%' sprintf('%ds',size(refchar,2)+4) '\n'], ...
            'Mode #','Freq (Hz)',cref{:},'Ref',under{:});
        for k=1:numroots
            str=sprintf('%11d%11.3f',rootlist(k,1),shpkeeptmp(k).frequency);
            for m=1:nref
                if ~isinf(cctmp(k,m))
                    str=[str sprintf('%10.3f',cc(k,m))];
                else
                    str=[str sprintf('%10s',sprintf('#%.3f#',cc(k,m)))];
                end
                if m==refpicktmp(k)
                    str=[str '*'];
                else
                    str=[str ' '];
                end
            end
            str=[str sprintf(fmtref,char(refcoord(refpicktmp(k))))];
            fprintf('%s\n',str);
        end
    end

    % Store shape if desired

    if fittype==2
        fprintf('** Selecting shapes based on highest correlation coefficient\n');
        shpkeep=shpkeeptmp;
        refpick=refpicktmp;
    end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Method 3:  Best correlation between FRF

    % Synthesize FRF from the fit shapes from each reference

    frall=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2))';
    deltaF=diff(frall(1:2));
    wrall=2*pi*frall(:);
    dpfrf=imat_fn(nref);
    for k=1:nref
        dpfrf(k)=ss.fe(dpind(k),k);
        foall(:,k)=dpfrf(k).ordinate(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
        %   ff=ss.fe(dpind(k),k);
        %   ff.ordinate=ff.ordinate(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
        %   [tmp,indp{k},indps{k},indpe{k}]=findpeaks(ff,.3,[],'rel',diagtoggle);
    end

    j=complex(0,1);
    nfreq=length(frall);

    shpkeeptmp=imat_shp(numroots);
    cc=zeros(numroots,nref);
    cckeep=zeros(numroots,nref);

    % Loop over each root

    for k=1:numroots
        wn=2*pi*shape(k,1).frequency;
        zeta=shape(k,1).damping;

        % Figure out how many spectral lines to do

        %   maxpts=15;        % MAXIMUM 41--15 on each side of the peak
        minpts= 2;        % MINIMUM  5-- 2 on each side of the peak

        [~,ind]=min(abs(wrall-wn));

        %   sind=ind;
        %   eind=ind;

        % Determine the optimal range of points to use

        %   for m=1:nref
        %       [tmp,indd]=min(abs(indp{m}-ind));
        %       if tmp<=2,
        %           sind=min(sind,indps{m}(indd));
        %           eind=max(eind,indpe{m}(indd));
        %       end
        %   end

        % Calculate number of spectral lines to half power points (actually go
        % twice as far as the half-power points)

        df=shape(k,1).frequency*shape(k,1).damping;
        nspec=ceil(df/deltaF);
        sind=ind-nspec;
        eind=ind+nspec;

        % Make sure we use more than one point in the comparison

        sind=min(sind,ind-minpts);
        eind=max(eind,ind+minpts);

        % Make sure we don't use too many points

        %   sind=max(sind,ind-maxpts);
        %   eind=min(eind,ind+maxpts);

        % Make sure we are still within bounds

        if sind<1
            eind=eind+(1-sind);
            sind=1;
        end
        if eind>nfreq
            sind=sind-(nfreq-sind);
            eind=nfreq;
        end

        % Calculate the baseline FRF

        wna=wn(ones(nfreq,1));
        zetaa=zeta(ones(nfreq,1));

        j=sqrt(-1);
        if ss.realcomplex==1		% Real normal modes
            Hk=-wrall.^2 ./ (wna.^2 + 2*j*zetaa.*wna.*wrall - wrall.^2);

        else				        % Complex modes
            lambda=j.*wna.*sqrt(1-zetaa.^2);
            Hk(:,1)= -wrall.^2 ./ (j*wrall + zetaa.*wna - lambda);
            Hk(:,2)= -wrall.^2 ./ (j*wrall + zetaa.*wna + lambda);
        end

        % Loop over each reference to see which FRF best matches the synthesized FRF

        Rmax=[-inf 0];

        frf=dpfrf;
        frf.abscissa=frall;

        for m=1:nref
            if ss.realcomplex==1
                sfrf=Hk*dp(k,m)^2;
            else
                sfrf=Hk(:,1)*dp(k,m)^2 + Hk(:,2)*conj(dp(k,m))^2;
            end
            sfrf=dpsign(m)*sfrf;

            frf(m).ordinate=sfrf;

            % Calculate correlation coefficient between this function and the
            % analytical function

            R = corrcoef([abs(sfrf(sind:eind)) abs(foall(sind:eind,m))]);
            cc(k,m)=R(1,2);
            cckeep(k,m)=R(1,2);
            cc(negdpmask)=-inf;		% Throw away shapes with negative drive point coefficients
            if cc(k,m)>Rmax(1)
                Rmax=[cc(k,m) m];
            end

        end

        if Rmax(2)==0
            Rmax(2)=1;	% Just select the first one.
        end

        if diagtoggle

            % Figure out frequency window

            %xwin=frall(ind)+[min(-3,-0.1*frall(ind)) max(3,0.1*frall(ind))];
            xwin=[frall(max(1,sind-2*nspec)) frall(min(length(frall),eind+2*nspec))];

            % Plot the original analytical FRF and overlay fits

            fmax = frf(Rmax(2));
            fmax.ordinate = fmax.ordinate(sind:eind);
            fmax.abscissa = frall(sind:eind);

            h = plot(dpfrf,frf,fmax,'xwin',xwin,'complex','mp');
            set(h.hf,'Name',sprintf('Mode %d: %g Hz, %g%% damping',k,shape(k,1).frequency,shape(k,1).damping*100), ...
                'NumberTitle','off');

            % Make the line colors match

            hl=sort(get(gca,'children'));
            for m = 1:nref
                h.style('Linestyle','--','Color',get(hl(m),'Color'),m+nref);
            end
            h.style('Marker','o','LineStyle','none','Color',get(hl(Rmax(2)),'Color'),nref*2+1)
           
            % Add legend

            legstr=cell(nref,1);
            for m=1:nref
                if Rmax(2)==m
                    str2=sprintf(' (CC=%.2f) <<O>>',Rmax(1));
                else
                    if ~isinf(cc(k,m))
                        str2=sprintf(' (CC=%.2f)',cc(k,m));
                    else
                       str2=sprintf(' #%.3f#',cckeep(k,m));
                    end
                end
                legstr{m}=sprintf('%s%s',char(refcoord(m)),str2);
            end
            h.legend(legstr,0);

        end

        % Store this shape

        shpkeeptmp(k)=shape(k,Rmax(2));
        refpicktmp(k)=Rmax(2);

        if diagtoggle
            hold off;
        end

    end

    % Display output:

    if fittype==3 || diagtoggle
        fprintf('\n');
        fprintf('Shapes selected by correlation coefficient from synthesized FRF (CC listed per reference)\n');
        fprintf('   NOTE:  CC''s surrounded by # indicate negative drive point residue, and\n');
        fprintf('          CC''s with a trailing * incidate the one for which the shape was kept.\n');
        fprintf('\n');

        cref=cellstr(refcoord);
        under=repmat({'-------'},1,nref+3);
        refchar=strjust(char(refcoord),'right');
        fmtref=['%' sprintf('%ds',size(refchar,2)+4)];

        fprintf([repmat('%11s',1,nref+2) '%' sprintf('%ds',size(refchar,2)+4) '\n'], ...
            'Mode #','Freq (Hz)',cref{:},'Ref',under{:});
        for k=1:numroots
            str=sprintf('%11d%11.3f',rootlist(k,1),shpkeeptmp(k).frequency);
            for m=1:nref
                if ~isinf(cc(k,m))
                    str=[str sprintf('%10.3f',cc(k,m))];
                else
                    str=[str sprintf('%10s',sprintf('#%.3f#',cckeep(k,m)))];
                end
                if m==refpicktmp(k)
                    str=[str '*'];
                else
                    str=[str ' '];
                end
            end
            str=[str sprintf(fmtref,char(refcoord(refpicktmp(k))))];
            fprintf('%s\n',str);
        end
    end

    % Store shape if desired

    if fittype==3
        fprintf('** Selecting shapes based on best correlation between FRF\n');
        shpkeep=shpkeeptmp;
        refpick=refpicktmp;
    end
end

%--------------------------------------------------------------------------
%% Store the output shapes

ss.shape.shape_all=shape;	% FIXME:  Remove when done
shape=shpkeep;

% Display a warning dialog if all DP coefficients were negative for a given root

dpsum=sum(negdpmask,2);
ind=find(dpsum==nref);
if any(ind)
    str=[];
    str{1}='Warning:  The following modes had negative drive point coefficients';
    str{2}='for all references.  These shapes should not be considered valid:';
    str{3}='';
    for k=1:length(ind)
        str{end+1}=sprintf('Mode %3d - %g Hz',ind(k),shpkeep(ind(k)).frequency);
    end

    uiwait(warndlg(str));
end

%--------------------------------------------------------------------------
%% Calculate the psi vectors

% Get the frequency range for the coefficient calculations

% xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
% if xf(1)==0
%   xf=xf(2:end);
% end
% 
% xw=2*pi*xf;
% bs=length(xf);

% Get the experimental FRF abscissa

wr=2*pi*ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2));
wr=wr(:).';
nel=length(wr);

% Get mode frequencies and damping

xf=shape.frequency(:);
zeta=shape.damping(:);
bs=length(xf);
xw=2*pi*xf;

% Calculate analytical SDOF FRF, assuming A=1, zeta
%    NOTE:  Each row is an FRF with resonant frequency at the spectral line
%           frequency through range of fit

%           This is Hp in equation (5) of the paper

j=sqrt(-1);
H = -(ones(bs,1)*wr.^2 ./ ...
     (xw.^2*ones(1,nel) + j*2*(zeta.*xw)*wr - ones(bs,1)*wr.^2));
 
% Loop over each reference

nref=length(ss.ref_coords);
psi2=zeros(size(ss.fe,1),length(shape),nref);

for k=1:nref

    % Compute the reciprocal modal vector Psi

    if ss.realcomplex==2	% Complex
        psi2(:,:,k) = ss.pinv(:,:,k)*H.';
    else			% Real normal
        Hs = [real(H) imag(H)];
        psi2(:,:,k) = ss.pinv(:,:,k)*Hs';
    end
end

% Now just keep the psi vector for the appropriate reference

psi=zeros(size(ss.fe,1),length(shape));
for k=1:length(shape)
    psi(:,k)=psi2(:,k,refpick(k));
end
