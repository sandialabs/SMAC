function hh = bar3(varargin)
%BAR3   3-D bar graph.
%   BAR3(Y,Z) draws the columns of the M-by-N matrix Z as vertical 3-D
%   bars.  The vector Y must be monotonically increasing or
%   decreasing. 
%
%   BAR3(Z) uses the default value of Y=1:M.  For vector inputs,
%   BAR3(Y,Z) or BAR3(Z) draws LENGTH(Z) bars.  The colors are set by
%   the colormap.
%
%   BAR3(Y,Z,WIDTH) or BAR3(Z,WIDTH) specifies the width of the
%   bars. Values of WIDTH > 1, produce overlapped bars.  The default
%   value is WIDTH=0.8
%
%   BAR3(...,'detached') produces the default detached bar chart.
%   BAR3(...,'grouped') produces a grouped bar chart.
%   BAR3(...,'stacked') produces a stacked bar chart.
%   BAR3(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
%
%   H = BAR3(...) returns a vector of surface handles.
%
%   Example:
%       subplot(1,2,1), bar3(peaks(5))
%       subplot(1,2,2), bar3(rand(5),'stacked')
%
%   See also BAR, BARH, and BAR3H.

%   Mark W. Reichelt 8-24-93
%   Revised by CMT 10-19-94, WSun 8-9-95
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2004-11-01 14:01:56 -0700 (Mon, 01 Nov 2004) $
%   Adapted by Dan Hensley / ATA Engineering / 05-26-04

narginchk(1,4);

%%%% New DPH %%%%
mat=varargin{1};
cmap=[0.8,  1 0 0 ; ...	% Red
      0.6,  1 1 0 ; ...	% Yellow
      0.4,  0 1 0 ; ...	% Green
      0.2,  0 1 1 ; ...	% Cyan
      0.0,  0 0 1 ] ;	% Blue
if max(mat(:))-min(mat(:))>1,
  crange=[0 100];
else
  crange=[0 1];
end
%%%% End DPH %%%%

[msg,x,y,xx,yy,linetype,plottype,barwidth,zz] = makebars(varargin{:},'3');

if ~isempty(msg)
    error(msg);
end

[n,m] = size(y);
% Create plot
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;
edgec = get(gcf,'defaultaxesxcolor');
facec = 'flat';
h = []; 
cc = ones(size(yy,1),4);

if ~isempty(linetype)
    facec = linetype;
end

for i=1:size(yy,2)/4

    % Set the color based on the amplitude (all vertices) DPH
    co=mat(:,i).';
    co=co(ones(6,1),:);
    co=co(:);
    co=co(:,ones(1,4));
    % End new

    h = [h,surface('xdata',xx+x(i),...
            'ydata',yy(:,(i-1)*4+(1:4)), ...
            'zdata',zz(:,(i-1)*4+(1:4)),...
            'cdata',co, ...
...%             'cdata',i*cc, ...
            'FaceColor',facec,...
            'EdgeColor',edgec,...
            'tag','bar3',...
            'parent',cax)];
end

if length(h)==1
    set(cax,'clim',[1 2]);
end

if ~hold_state, 
  % Set ticks if less than 16 integers
  if all(all(floor(y)==y)) & (size(y,1)<16) 
      set(cax,'ytick',y(:,1));
  end
  
 xTickAmount = sort(unique(x(1,:))); 
 if length(xTickAmount)<2
     set(cax,'xtick',[]);
 elseif length(xTickAmount)<=16
      set(cax,'xtick',xTickAmount);
 end  %otherwise, will use xtickmode auto, which is fine
  
  hold off, view(3), grid on
  set(cax,...
      'NextPlot',next,...
      'ydir','reverse');
  if plottype==0,
    set(cax,'xlim',[1-barwidth/m/2 max(x)+barwidth/m/2])
  else
    set(cax,'xlim',[1-barwidth/2 max(x)+barwidth/2])
  end
  grid on
  dx = diff(get(cax,'xlim'));
  dy = size(y,1)+1;
  if plottype==2,
    set(cax,'PlotBoxAspectRatio',[dx dy (sqrt(5)-1)/2*dy])
  else
    set(cax,'PlotBoxAspectRatio',[dx dy (sqrt(5)-1)/2*dy])
  end
end

if nargout==1
    hh = h;
end
