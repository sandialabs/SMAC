function nroots=smac_getnumroots(ss,freq,miftype)
% SMAC_GETNUMROOTS  Get number of repeated roots at the supplied frequencies
%
% NROOTS=SMAC_GETNUMROOTS(SS,FREQ,MIFTYPE)
%
% SMAC_GETNUMROOTS will attempt to automatically pick the number of
% repeated roots at the specified frequency list.  SS is the SMAC data
% structure, FREQ is a vector of frequencies to check, and MIFTYPE is a
% string giving the type of MIF to use for the check ('nmif', 'mmif', or
% 'cmif').  It will return the number of roots in NROOTS, which is a vector
% the same length as FREQ.
%
% NMIF is generally not useful for multiple reference data.
% SMAC_GETNUMROOTS assumes that the corresponding MIF has already been
% calculated and is present in the SS structure.

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
%  03-Sep-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  01-Oct-2004 / ATA Engineering / Dan Hensley
%    o Change algorithm so that it counts peaks in each of the indicator
%      functions
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for new .mif fields
%
%==========================================================================

debug=false;

% Input error checking

narginchk(3,3);
if ~isstruct(ss),
  error('SS must be a SMAC structure');
end
if ~isnumeric(freq),
  error('Frequency vector must be numeric');
end
if ~ischar(miftype),
  error('MIFTYPE must be a string: ''nmif'', ''mmif'', or ''cmif''');
end

%--------------------------------------------------------------------------

freq=freq(:);

% fprintf('...Searching for multiple roots at frequencies: ');
% fprintf('%g, ',freq(1:end-1));
% fprintf('%g\n',freq(end));

% Find the closest index for each frequency into the data

xf=ss.fe(1).abscissa;

freqind=zeros(size(freq));
for k=1:length(freq);
  [~,freqind(k)]=min(abs(xf-freq(k)));
end

% Get the indices into the correlation coefficient frequencies

xfc=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));

fcind=zeros(size(freq));
for k=1:length(freq);
  [~,fcind(k)]=min(abs(xfc-freq(k)));
end

% Set the output equal to the current input

[tmp,ind]=ismember(fcind,ss.corr.corrind);
ind=ind(ind~=0);
nroots=ones(size(freq));
ss.corr.nroots(length(ss.corr.corrind)+1)=0;
nroots(tmp)=ss.corr.nroots(ind);

% Count roots based on MIF type

switch miftype

  case 'nmif',
    if length(ss.ref_coords)>1,
      uiwait(errordlg('Cannot look for repeated roots with NMIF'));
      return;
    end
    nroots=ones(size(freq));


  case 'mmif',

    %nspec=1;
    freqfrac=0.50/100;		% Frequency must be within this range to be considered near the peak
    nroots=zeros(size(freq));

    % Invert the MMIF

    ws=warning;
    warning('off');
      mmifi=1/ss.mif.exp.mmif;
    warning(ws);
    xf=mmifi(1).abscissa;
    ainc=diff(xf(1:2));

    % First find all the peaks in each of the indicator functions
    %    - store in a cell array

    % Find the peaks in each of the indicator functions

    nrefs=length(ss.ref_coords);
    peaks=cell(nrefs,1);
    for k=1:nrefs,
      [~,peaks{k}]=findpeaks(mmifi(k),0.03,[],'rel');
      if debug,
        plot(ss.mif.exp.mmif(k));
        hold on;
          plot(ss.mif.exp.mmif(k).abscissa(peaks{k}),ss.mif.exp.mmif(k).ordinate(peaks{k}),'r* ');
        hold off;
      end
    end

    % For each frequency passed in, see if we are within FREQFRAC
    % percentage from one of the listed peaks.  Cycle through each of
    % the MIF and count the same.

    for k=1:nrefs,
      for m=1:length(freq),
        [fdiff,ind]=min(abs(xf(peaks{k})-freq(m)));
        %if round(fdiff/ainc)<=nspec,
        if fdiff/freq(m)<=freqfrac,
          if k>1
            % To be repeated,
            %    The higher MIF must drop to at least 60% of the value
            %    The MMIF value must also be below 0.9
            %    The primary MMIF must have also shown a peak
            if (1-ss.mif.exp.mmif(k).ordinate(peaks{k}(ind)))/(1-ss.mif.exp.mmif(k-1).ordinate(peaks{k}(ind)))>=0.60 & ...
               ss.mif.exp.mmif(k).ordinate(peaks{k}(ind))<0.9 & ...
               nroots(m)>0,	% We must have first counted 
              nroots(m)=nroots(m)+1;
            end
          else
            nroots(m)=nroots(m)+1;
          end
        end
      end
    end

    % For diagnostics, plot the MMIF with the number of roots counted at
    % each

    if debug,
      ind=[];
      for k=1:length(freq),
        [tmp,ind(k)]=min(abs(xf-freq(k)));
      end
      plot(ss.mif.exp.mmif);
      hold on;
      plot(freq,ss.mif.exp.mmif(1).ordinate(ind),'r*');
      text(freq,ss.mif.exp.mmif(1).ordinate(ind)-.02,cellstr(num2str(nroots)));
    end

    nroots(nroots==0)=1;

%     threshval=0.5;
% 
%     mifvals=ss.mif.exp.mmif(1:length(ss.ref_coords)).ordinate(freqind,:);
%     % FIXME:  Ignore any initial ones above 0.9??
%     mifthresh=repmat(mifvals(:,1),1,length(ss.ref_coords)) + ...
%               repmat((1-mifvals(:,1)),1,length(ss.ref_coords)) .* ...
%               repmat([0 1-threshval.^(1:length(ss.ref_coords)-1)],length(freq),1);
%     nroots=sum(mifvals<=mifthresh,2);


  case 'cmif',

    %nspec=1;
    freqfrac=0.50/100;		% Frequency must be within this range to be considered near the peak
    nroots=zeros(size(freq));

    xf=ss.mif.exp.cmif(1).abscissa;
    ainc=diff(xf(1:2));

    % First find all the peaks in each of the indicator functions
    %    - store in a cell array

    % Find the peaks in each of the indicator functions

    nrefs=length(ss.ref_coords);
    peaks=cell(nrefs,1);
    for k=1:nrefs,
      [~,peaks{k}]=findpeaks(ss.mif.exp.cmif(k),0.03,[],'rel');
      if debug,
        plot(ss.mif.exp.cmif(k));
        hold on;
          plot(ss.mif.exp.cmif(k).abscissa(peaks{k}),ss.mif.exp.cmif(k).ordinate(peaks{k}),'r* ');
        hold off;
      end
    end

    % For each frequency passed in, see if we are within FREQFRAC
    % percentage from one of the listed peaks.  Cycle through each of
    % the MIF and count the same.

    for k=1:nrefs,
      for m=1:length(freq),
        [fdiff,ind]=min(abs(xf(peaks{k})-freq(m)));
        %if round(fdiff/ainc)<=nspec,
        if fdiff/freq(m)<=freqfrac,
          if k>1
            % To be repeated,
            %    The lower MIF must be within 10% of the higher value
            %    The primary MIF must have also shown a peak
            if ss.mif.exp.cmif(k).ordinate(peaks{k}(ind))/ss.mif.exp.cmif(k-1).ordinate(peaks{k}(ind))>0.2 & ...
               nroots(m)>0,	% We must have first counted 
              nroots(m)=nroots(m)+1;
            end
          else
            nroots(m)=nroots(m)+1;
          end
        end
      end
    end

    % For diagnostics, plot the CMIF with the number of roots counted at
    % each

    if debug,
      ind=[];
      for k=1:length(freq),
        [tmp,ind(k)]=min(abs(xf-freq(k)));
      end
      plot(ss.mif.exp.cmif);
      hold on;
      plot(freq,ss.mif.exp.cmif(1).ordinate(ind),'r*');
      text(freq,ss.mif.exp.cmif(1).ordinate(ind)-.02,cellstr(num2str(nroots)));
    end

    nroots(nroots==0)=1;

%     threshval=0.5;
% 
%     mifvals=ss.mif.exp.cmif(1:length(ss.ref_coords)).ordinate(freqind,:);
%     mifthresh=repmat(mifvals(:,1),1,length(ss.ref_coords)).* ...
%               repmat([1 threshval.^(1:length(ss.ref_coords)-1)],length(freq),1);
%     nroots=sum(mifvals>=mifthresh,2);


  otherwise
    error('Unknown MIF Type:  choose from NMIF, MMIF, or CMIF');

end
