function smac_autofit(showconvplot)
% SMAC_AUTO_FIT  Perform auto-SMAC curve fit by iterating over frequency and damping
%
% SMAC_AUTOFIT([DOGRAPH])
%
% SMAC_AUTO_FIT takes the initial root list selected from the correlation
% coefficient calculation and estimates each root by oscillating between a
% frequency and damping optimization until the correlation coefficient is
% maximized based on the damping and frequency convergence tolerances
% specificed in this function.
%
% DOGRAPH is an optional logical specifying whether to display the
% convergence graphs for each root while autofitting each root

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
%  25-May-2004 / ATA Engineering / Dan Hensley
%    o Initial creation; adapted from auto_smac_c.m
%    o Renamed variables, reorganized, added comments
%
%  30-May-2004 / ATA Engineering / Dan Hensley
%    o Sort non-zero roots at the end
%
%  10-Sep-2004 / ATA Engineering / Dan Hensley
%    o Add interior loop to loop over multiple roots at a given peak
%    FIXME:  More work needs to be done with convergence and root keeping
%    o Add max # of iterations
%
%  24-Sep-2004 / ATA Engineering / Dan Hensley
%    o When splitting a peak due to multiple roots, put the split
%      frequencies into the original root list
%    o Fix bug that occurred when expanding the original root list
%    o Add original damping to the information stored
%    o Change root splitting algorithm to what I discussed with Randy M
%      today (split into 3, +/- Df/2, and don't do an initial frequency
%      fit on the Df splits)
%    o Add graphical debugging measure showing frequency and damping
%      convergence
%
%  15-Oct-2004 / ATA Engineering / Dan Hensley
%    o Add toggle for convergence graphs as input argument
%    o Add correlation coefficient graph to the convergence graphs
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Update for changed autofit algorithms (loop over references)
%
%  07-Mar-2005 / ATA Engineering / Dan Hensley
%    o Calculate actual CC for split roots instead of using the initial
%      root's value
%
%==========================================================================

global ss;

% Check input arguments

narginchk(0,1);
if ~exist('showconvplot','var'),
  showconvplot=false;
end
if ~islogical(showconvplot) || numel(showconvplot)~=1,
  error('Convergence graph display toggle must be a logical scalar');
end

maxiter=200;

% Set convergence multipliers 

fconmult=1.5;
% dconmult=1.5;
dconmult=2.0;
% FIXME:  Change fconmult to 2.0 and dconmult to 3.0 based on conversation
% with Randy 09/26/05

%--------------------------------------------------------------------------
% Extract the information that is constant in these calculations
%  Used in smac_ffit_cauto and smac_dfit_cauto

wr=2*pi*ss.fe(1).abscissa(ss.freqrange(1,2):ss.freqrange(2,2)).';
frf=ss.fe.ordinate(ss.freqrange(1,2):ss.freqrange(2,2),:);
frf=reshape(frf,[size(frf,1) size(ss.fe)]);

% Extract information from the structure

xf=ss.fe(1).abscissa(ss.freqrangecc(1,2):ss.freqrangecc(2,2));
peaks=xf(ss.corr.corrind)';
initcorr=ss.corr.corr(ss.corr.corrind);

% Initialize variables

%rootlist=zeros(length(peaks),3);
rootlist=zeros(0,4);
corr_ref=zeros(0,length(ss.ref_coords));

%niter=ones(1,length(peaks));
niter=[];

% Store original information

rootlistorig=[peaks(:) ss.corr.zeta*ones(length(peaks),1) initcorr(:)];

%--------------------------------------------------------------------------
% Loop over each peaks to get optimal frequency and damping values

imat_progress(0,'Auto Fit Progress');
rootiter=0;
originc=1;

for kk=1:length(peaks)

  imat_progress(2,sprintf('Peak %d (%g Hz)',kk,peaks(kk)));
  icor=initcorr(kk);

  if ss.corr.nroots(kk)>1,  % Multiple roots:  split and insert into table
    fprintf(' >>>>>Multiple roots at peak %d:  splitting to generate multiple seeds and\n',kk);
    fprintf('                                  tightening convergence tolerance\n');
    imat_progress(2,sprintf('Peak %d (%g Hz) --multiple roots--',kk,peaks(kk)));

% ORIGINAL SPLIT ALGORITHM, DIDN'T WORK THAT WELL
%     mm=5;
%     freqp=ss.fit.freqp;
%     peakval=peaks(kk)+[0 peaks(kk)*freqp/100*[-1 -1 1 1]];
%     zeta=[ss.corr.zeta ss.fit.zeta(1)*[1 1] ss.fit.zeta(2)*[1 1]];	% FIXME:  Use different seed damping values??
%     diffz=diff(ss.fit.zeta)/2;
%     zetar=[ss.fit.zeta;	...	% Make a matrix
%            ss.fit.zeta(1)/2 ss.fit.zeta(1)+diffz; ...
%            ss.fit.zeta(2)-diffz ss.fit.zeta(2)*1.1; ...
%            ss.fit.zeta(1)/2 ss.fit.zeta(1)+diffz; ...
%            ss.fit.zeta(2)-diffz ss.fit.zeta(2)*1.1   ];

    mm=3;
    freqp=ss.fit.freqp;
    %peakval=peaks(kk)+[0 peaks(kk)*(freqp/2/100)*[-1 1]];
    %zeta=[ss.corr.zeta; ss.fit.zeta(2)*ones(mm-1,1)];
    %zeta=[ss.corr.zeta; mean(ss.fit.zeta)*ones(mm-1,1)];
    peakval=peaks(kk)+(peaks(kk)*(freqp/2/100)*[-1 0 1]);
    zeta=[ss.fit.zeta(1); ss.corr.zeta; ss.fit.zeta(2)];
    zetar=repmat(ss.fit.zeta,mm,1);
    corr=ss.corr.corr(ss.corr.corrind(kk))*ones(size(peakval));

    % Calculate initial correlation for repeated roots
    
    for m=2:length(peakval);
        fmm(1:2)=peakval(m);
        dmm(1)=zeta(m);
        numpts=1;
        smac_ffit_cauto;
        corr(m)=max_cc;
    end
    
    % Expand the original root list

    rootlistorig=[rootlistorig(1:originc-1,:); ...
                  repmat(rootlistorig(originc,:),mm,1); ...
                  rootlistorig(originc+1:end,:) ];
    rootlistorig(originc:originc+mm-1,:)=[peakval(:) zeta(:) corr(:)];
    originc=originc+mm;

    % Specify damping and frequency convergence values (use half of nominal)

%     fcon=ss.fit.freqconv*.5;
%     dcon=ss.fit.zetaconv*.5;
    fcon=ss.fit.freqconv;
    dcon=ss.fit.zetaconv;

  else                  % Single root
    mm=1;
    peakval=peaks(kk);
    freqp=ss.fit.freqp;
    zeta=ss.corr.zeta;
    zetar=ss.fit.zeta;
    corr=ss.corr.corr(ss.corr.corrind(kk));

    % Specify damping and frequency convergence values

    fcon=ss.fit.freqconv;
    dcon=ss.fit.zetaconv;
    originc=originc+1;

  end

  %------------------------------------------------------------------------
  % Loop over all of the potential roots

  for mm=1:length(peakval),

    % Set up convergence graph if necessary

    if showconvplot,
      hf=figure('color','white','backingstore','on','doublebuffer','on');
      set(hf,'Name',sprintf('Root %d:  Start at %.3f Hz, %.4f damping',size(rootlist,1)+1,peakval(mm),zeta(mm)));
      subplot(3,1,1);
      hlc=plot([-1 -1],corr(mm)*[1 1],'b-');
      set(hlc,'Linewidth',2);
      hold on;
      xlim([-1 0]);
      xlabel('Iteration');
      ylabel('Correlation Coefficient');

      subplot(3,1,2);
      hlf=plot([-1 -1],peakval(mm)*[1 1],'b-');
      set(hlf,'Linewidth',2);
      hold on;
      plot(repmat([-1 maxiter+1]',1,2),repmat(peakval(mm) + peakval(mm).*[-1 1]/2*(freqp/100),2,1),'g--');
      hold off;
      xlim([-1 0]);
      xlabel('Iteration');
      ylabel('Frequency (Hz)');

      subplot(3,1,3);
      hld=plot([-1 -1],zeta(mm)*[1 1],'r-');
      set(hld,'Linewidth',2);
      hold on;
      plot(repmat([-1 maxiter+1]',1,2),repmat(zetar(mm,:),2,1),'g--');
      hold off;
      xlim([-1 0]);
      xlabel('Iteration');
      ylabel('Damping (ratio)');
      ha=get(hf,'Children');
      
      drawnow;
    end

    % Set iteration counter

    rootiter=rootiter+1;
    niter(rootiter)=1;

    % Perform the initial frequency fit using the default damping value

    if mm==1,			% Only do this for the first root of any multi-root splits
      fmm=peakval(mm) + peakval(mm).*[-1 1]/2*(freqp/100);
      dmm=zeta(mm);		% Use the default damping fraction for the frequency iteration
      numpts=ss.fit.freqpts;
      
      smac_ffit_cauto
    else
      freq_opt=peakval(mm);
    end

    % Perform the damping fit using the optimal frequency from the frequency fit

    dmm=zetar(mm,:);		% Damping range
    fmm=freq_opt;		% Use the optimal frequency from smac_ffit_cauto
    numpts=ss.fit.zetapts;

    smac_dfit_cauto

    %------------------------------------------
    % Initialize variables to refine the search

    freq_old=inf;
    zeta_old=inf;
    rootcor=max_cc;   	%  changed from Mcc(3) to max(Mcc) aug 1 03 rlm

    numpts=3;   %  This reduces the refinement in ffit and dfit to the tolerances fcon and dcon

    % Add this data to the convergence plot

    if showconvplot,
      xx=get(hlc,'XData');
      yy=get(hlc,'YData');
      xx=[xx(1); max(xx)+1];
      yy=[yy(1); rootcor];
      set(hlc,'XData',xx,'YData',yy);

      xx=get(hlf,'XData');
      yy=get(hlf,'YData');
      xx=[xx(1); max(xx)+1];
      yy=[yy(1); freq_opt];
      set(hlf,'XData',xx,'YData',yy);
      
      yy=get(hld,'YData');
      yy=[yy(1); zeta_opt];
      set(hld,'XData',xx,'YData',yy);
      set(ha,'XLim',[-1 max(xx)+1]);
      
      drawnow;
    end

    % Perform increasingly refined fits until the damping and frequency do
    % not change

    while ( zeta_opt~=zeta_old || freq_opt~=freq_old )

      niter(rootiter)=niter(rootiter)+1;
      icor=rootcor;

      % Remember last frequency and damping value

      zeta_old=zeta_opt;
      freq_old=freq_opt;

      % Iterate over frequencies

      fmm=freq_opt*(1+fcon*[-1 1]);
      dmm=zeta_opt;

      smac_ffit_cauto

      % Iterate over damping values

      fmm=freq_opt;
      dmm=zeta_opt*(1+dcon*[-1 1]);

      smac_dfit_cauto

      % Remember the last optimal correlation coefficient

      rootcor=max_cc;
      %fprintf('......Iteration %d:  freq=%.3f, damping =%.4f\n',niter(rootiter),freq_opt,zeta_opt);

      % Add this data to the convergence plot

      if showconvplot,
	    xx=get(hlc,'XData');
	    yy=get(hlc,'YData');
	    xx=[xx(:); max(xx)+1];
	    yy=[yy(:); rootcor];
	    set(hlc,'XData',xx,'YData',yy);

	    xx=get(hlf,'XData');
	    yy=get(hlf,'YData');
	    xx=[xx(:); max(xx)+1];
	    yy=[yy(:); freq_opt];
	    set(hlf,'XData',xx,'YData',yy);
      
	    yy=get(hld,'YData');
	    yy=[yy(:); zeta_opt];
	    set(hld,'XData',xx,'YData',yy);
        set(ha,'XLim',[-1 max(xx)+1]);
        pause(0);
      
        drawnow;
      end

      % Bail out if we reach the maximum number of iterations

      if niter(rootiter)>maxiter,
        fprintf('***** Number of iterations exceeded %d:  breaking\n',maxiter);
        break;
      end

    end

    %  diagnostic plots below of comparison of H analytical and Hp after optimization
    %   Also you need to turn on the "figure" command and the "close" command
    %   at the beginning and end of this routine
    %     semilogy(freq_rad(ind-s.nl:ind+s.nl)/(2*pi),abs(H(k,ind-s.nl:ind+s.nl).'),'r',freq_rad(ind-s.nl:ind+s.nl)/(2*pi),abs(Hp(ind-s.nl:ind+s.nl,k)),'b')   %  What is frequency vector called
    %         legend('analytical','experimental')
    %         xlabel('Frequency - Hz')
    %         pause    

    % Display information about the optimal values for this peak

    icor=rootcor;       % Before this point icor is the next to last value
    fprintf('...Peak %d: Frequency=%g Hz, damping=%.3f%% corr coeff=%.3f (%d iterations)\n', ...
            kk,freq_opt,zeta_opt*100,icor,niter(rootiter));

    % Throw away peak if we diverge too much on frequency or the damping is
    % negative or the frequency is out of the fit range

    if abs((freq_old-peakval(mm))/peakval(mm))>.03 || ...
            zeta_old<=0 || ...
            zeta_opt<0 || ...
            freq_opt<ss.freqrangecc(1,1) || ...
            freq_opt>ss.freqrangecc(2,1)
      freq_old=0;
      zeta_old=0;
      
      freq_opt=0;   % Added 3/14/06 DPH
      zeta_opt=0;
      icor=0;
      fprintf('......Removing this root--diverged too much or frequency is out of bounds\n');
    end

    % Store the information for this peak

    %rootlist(kk,:)=[freq_opt zeta_opt icor];
    rootlist(end+1,:)=[freq_opt zeta_opt icor max_ref];
    corr_ref(end+1,:)=cormax(1,:);

    % Review this root compared to the others and remove it if it is
    % identical within our tolerances

    rr=rootlist(end,:);
    wv=warning;
	warning off MATLAB:divideByZero;
    ind=abs((rootlist(:,1)-rr(1))/rr(1))<fconmult*fcon & abs((rootlist(:,2)-rr(2))/rr(2))<dconmult*dcon;
    warning(wv);
    ind(end)=false;
    if any(ind),
      fprintf('......Removing this root--duplicate\n');
    end
    rootlist(ind,:)=0;

  end

  % Update the progress bar

  imat_progress(1,kk/length(peaks));

end
imat_progress(-1);

%--------------------------------------------------------------------------
% Remove identical roots from this root list, where
%     Frequency is less than 1.5*fcon AND damping is less than 1.5*dcon

% FIXME:  This loop removes the first root that matches the tolerance.  Do
% we want to do this or keep the first and remove the rest?  Maybe just
% keep the one with the highest tolerance?

% FIXME:  How to handle multiple roots??

% for kk=1:size(rootlist,1)
%   rr=rootlist(kk,:);
%   if rr(1)~=0,
%     ind=abs((rootlist(:,1)-rr(1))/rr(1))<fconmult*fcon & abs((rootlist(:,2)-rr(2))/rr(2))<dconmult*dcon;
%     ind(kk)=false;
%     if any(ind),
%       ind=kk;
%       fprintf('...Removing root %d: %g Hz, damping=%g\n',kk,rootlist(ind,1:2));
%     end
%     rootlist(ind,:)=0;
%   end
% end

%--------------------------------------------------------------------------

% Sort the root list (keep 0 frequency roots in place)

tmpind=find(rootlist(:,1)~=0);
[tmp,ind]=sortrows(rootlist(tmpind,:));
rootlist(tmpind,:)=tmp;

% Store the initial and final root list in the global structure

ss.fit.rootlist=rootlist;
ss.fit.rootlistorig=rootlistorig;
ss.fit.corr_ref=corr_ref;

return;
