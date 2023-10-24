function hf=smac_splash(time)
% SMAC_SPLASH  Display SMAC splash screen

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
%  06-Jul-2004 / ATA Engineering / Dan Hensley
%    o Automatically remove the splash screen after 20 seconds if it's
%      still there
%
%==========================================================================

% Figure out where the splash screen is located

pname=which(mfilename);
ind=findstr(pname,filesep);
pname=pname(1:ind(end));

fname=fullfile(pname,'snllineblk.png');
if ~exist(fname,'file'),
  fprintf('*** Error:  Could not find logo file ''%s''\n',fname);
  return;
end

% Load the file

ss=get(0,'ScreenSize');
im=imread(fname);
if isempty(im),
  return;
end

% Display the image

hf=figure('Color','white', ...
          'NumberTitle','off', ...
          'MenuBar','none', ...
          'Name',' ', ...
          'Visible','off');

% Put the SMAC text below the image

ui=uiWidgets;

ht=ui.CreateText(hf,'SMAC (Synthesize Modes And Correlate)',[],[1 1]);
ht.HorizontalAlignment='center';
ht.FontSize=20;
ht.FontName='Courier New';
ht.FontWeight='bold';
ht.BackgroundColor='white';

ui.FitText(ht);

% Display the image

image(double(im)/255);
imsize=[size(im,2) size(im,1)];

axis equal;
ha=handle(gca);
set(ha,'XTick',[],'YTick',[],'XColor','white','YColor','white','Box','off','Units','pixels')
ha.Position(3)=600*ss(3)/1600;
axis tight;
ha.Position(4)=ha.Position(3)*imsize(2)/imsize(1);
ha.Position(1:2)=[1 ht.Position(2)+ht.Position(4)+10];


% Resize the figure and move it to the upper right corner

ui.ResizeFigure(hf,0);
hf=handle(hf);
hf.Position(1:2)=[1 ss(4)-hf.Position(4)-25];

% Center the two controls

ht.Position(1)=(hf.Position(3)-ht.Position(3))/2;
ha.Position(1)=(hf.Position(3)-ha.Position(3))/2;

set(hf,'Visible','on');

% Store the figure handle for later

%delfigure(hf);
ht=timer('ExecutionMode','singleshot','StartDelay',20,'TimerFcn','1;','StopFcn',@delfigure,'UserData',hf);
start(ht);

return;

%==========================================================================
% AUXILLIARY FUNCTIONS
%==========================================================================
function delfigure(ht,event,tmp)

hf=get(ht,'UserData');
if ishandle(hf),
  delete(hf);
end
delete(ht);
