function fig = smac_helper(type)
% HELPER  Provide GUI containing help on the SMAC forms
%
% FIG=HELPER(TYPE)
%
% HELPER provides a listbox GUI containing the help for the form specified
% in TYPE.  TYPE is a string containing the desired help topic.
%
% FIG is the handle of the HELPER GUI.

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
%  14-May-2004 / ATA Engineering / Dan Hensley
%    o Remove dependence on helper.mat
%
%  08-Jun-2004 / ATA Engineering / Dan Hensley
%    o Add help for manual fit GUI
%    o Update help for the other GUIs
%
%  17-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add help for resynthesis
%
%  18-February-2008 / Sandia National Labs / Randy Mayes
%	o pseudo_inv - frequency range, data condensation
%	o nsmac_hlp - damping estimate, bandwidth of fit, frequency lines
%	o corre_plot 
%	o sf-hlp - frequency range, damping range
%	o sh_hlp - CMIF, MMIF
%
%  TO DO:
%    o Rename to smac_help
%    o Have listbox pulldown that allows you to select which help to view
%
%==========================================================================

ui=uiWidgets;
uip=ui.GetParams();

h0 = figure('Color',uip.Form.Color, ...
	'Units','pixels', ...
	...'Colormap',mat0, ...
	...'FileName','/export/home/seklenk/codes/smac_pac/helper.m', ...
	'Name','Help Window', ...
	...'PaperPosition',[18 180 576 432], ...
	...'PaperUnits','points', ...
	'Position',[10 300 700 755], ...
	'Tag','Helper', ...
	'NumberTitle','off', ...
	'MenuBar','none', ...
	'ToolBar','none');
help_hand = uicontrol('Parent',h0, ...
	'Units','pixels', ...
	'BackgroundColor',[1 1 1], ...
	...'ForegroundColor',[0 0 1], ...
	'FontName','Courier New', ...
	'FontSize',12, ...
	'HorizontalAlignment','left', ...
	...'Max',20, ...
	'Position',[5 50 690 700], ...
	'Style','listbox', ...
	'Tag','help_list');
h1 = uicontrol('Parent',h0, ...
	'Units','pixels', ...
	...'BackgroundColor',[0 1 1], ...
	'Callback','close', ...
	'FontSize',12, ...
	'FontWeight','bold', ...
	...'ForegroundColor',[1 0 0], ...
	'ListboxTop',0, ...
	...'Position',[0.416015625 0.07552083333333334 0.212890625 0.07031250000000001], ...
	'Position',[315 5 70 40], ...
	'String','OK', ...
	'Tag','out_button');
if nargout > 0, fig = h0; end

%--------------------------------------------------------------------------
switch type,
    %--------------
    case 'pseudo_inv',
	he.str = { ...
          'This interface is the initial SMAC window that allows the';
		  'the user to select an analysis method and to modify the';
		  'existing frequency response function (frf) data before';
 		  'beginning SMAC.'; 
		  '';
		  'SOLUTION METHOD';
		  'The user can select SMAC to fit the experimental data as';
		  'either REAL MODES or COMPLEX MODES. This selection does';
		  'change the form of the FRF matrix used in the pseudo-';
		  'inverse and ultimately in the SMAC analysis.  REAL modes';
		  'is always recommended.  If the results are not satisfactory,';
		  'then try complex modes'; 
		  '';
		  'FREQUENCY RANGE';
		  'These selection toggles allow the user to keep the full';
		  'experimental frequency range (FULL RANGE) in the SMAC ';
		  'analysis or to modify the range (PARTIALRANGE). Normally';
		  'full range is recommended unless there are more modes than';
		  'sensors, in which case the user might consider reducing the';
		  'range. The user can select the partial frequency range';
          'either by typing frequencies in, or by graphical selection';
          'using the Select button.  If REAL MODES are selected, the';
          'NMIF will be displayed.  If COMPLEX MODES are selected, the';
          'CMIF will be displayed. The pseudo-inverse is only calculated';
          'over the selected range.';
          '';
...		  'DATA CONDENSATION';
...		  'This pushbutton is not recommended because condensation';
...		  'of the FRF data generally degrades the modal parameter';
...		  'accuracy.  It may be used to reduce the size of the data';
...		  'set or focus on only well excited modes, but results can';
...		  'be unpredictable';
...		  'This pushbutton enables further condensation of the FRF';
...		  'data through a SINGULAR VALUE DECOMPOSITION (SVD). When the';
...		  'SINGULAR VALUE DECOMPOSITION button is selected, an edit';
...		  'box containing a normalized singular value cutoff can be';
...		  'modified before calculating the pseudo-inverse. A figure';
...		  'will be displayed showing a plot of the normalized singular';
...		  'values when the user executes the pseudo-inverse calculation.';
...		  'A question box is then displayed which allows the user to add';
...		  'or retain additional singular values in the analysis. Once';
...		  'the user selects a value (from 0 to [total singular value -';
...		  'singular value cutoff]), the pseudo-inverse is calculated.';
...		  '';
		  'PLOT FRFs';
		  'This pushbutton enables the user to view all the existing';
		  'FRF data that has been loaded into SMAC. Any number of data'; 
		  'channels can be overlayed to view response characteristics.';
		  '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.'; 
		  '';
		  'EXECUTE PSEUDO INVERSE';
		  'This pushbutton will execute the pseudo-inverse calculation';
		  'and will initiate a call to the SMAC Correlation Coefficient';
		  'interface.'}; 

    %--------------
    case 'nsmac_hlp',
	he.str = { ...
          'This interface is used to set initial calculation values and';
		  'the frequency of fit for the SMAC Correlation Calculations. ';
		  '';
		  'DAMPING ESTIMATE';
		  'An initial estimate of damping is required for the correlation';
		  'calculation. The user can select any value for the damping ';
		  'estimate. It is recommended that a value slightly below the ';
		  'expected average damping value be chosen. The default is 0.02 (2%).';
		  '';
		  'BANDWIDTH OF FIT';
		  'This allows the user to analyze or fit the experimental data over';
		  'a reduced frequency range of interest. These edit boxes are defaulted';
		  'to the frequency range values selected by the user for calculating';
		  'the pseudo-inverse. These frequencies should not exceed the';
		  'frequency range selected by the user in the previous interface.';
		  'For example: if the user selected a modified frequency ';
		  'range of 100 to 1200 Hz for calculating the pseudo-inverse,';
		  'then the bandwidth of fit CAN NOT be below 100 Hz or greater';
		  'than 1200 Hz. Also, the low frequency should always be selected';
		  'greater than 0 Hz ';
		  '';
		  'FREQUENCY LINES';
 		  'The correlation coefficients for SMAC are only calculated';
		  'near each resonant peak to enhance identification of the';
		  'damping coefficient, because the strongest effects of the';
		  'damping are seen near resonance. The NO. OF LINES is used';
		  'to select the number of frequency lines on each side of the';
		  'resonant frequency in making the correlation coefficient calculation.';
		  'For any frequency of interest, the number of lines should extend';
		  'past the half power point frequency.  For lightly damped systems,';
		  'SMAC needs many more lines than this to converge on the damping.';
		  'The default NO. OF LINES is set to 20. Using more lines helps ';
		  'convergence, however, more lines will slightly increase computation times.';
		  '';
		  'PLOT FRFs';
		  'This pushbutton enables the user to view all the existing';
		  'FRF data that has been loaded into SMAC. Any number of data'; 
		  'channels can be overlayed to view response characteristics.';
		  '';
		  'CALCULATE CORRELATION';
		  'This button executes the correlation coefficient calculation';
		  'given the initial parameters and frequency bandwidth. A plot';
		  'window will be initiated showing the calculated correlation';
		  'coefficients versus analysis or bandwidth of fit frequencies.';
		  '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.'; 
		  '';
		  'BACKUP';
		  'The BACKUP button allows the user to move back to the pseudo-';
		  'inverse interface to re-calculate this value if desired.'};

    %--------------
    case 'corre_plot',
	he.str = { ...
          'This interface plots the results of the SMAC Correlation';
		  'Coefficient calculation. The lower plot shows the correlation';
		  'coefficient values plotted versus the user selected fit or';
		  'analysis frequency range.  The plot on the top can be toggled';
          'between NMIF MMIF & CMIF to help visualize resonant frequencies.';
          'The minimum coefficient is used to automatically select peaks';
          'to be included in the initial root list.  The user can either';
          'type in a minimum coefficient or select it graphically  from the';
          'plot.  The selected roots and their associated correlation';
          'coefficient are shown in the listbox on the left.  Selected';
          'roots are highlighted on the plot with a red asterisk (*).';
		  '';
		  'The user can continue to try different minimum correlation';
		  'coefficient values until the list of peak frequencies are';
		  'satisfactory for SMAC fitting. The user must then select the';
		  'INITIATE AUTO SMAC button to continue the analysis.';
		  '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.'; 
		  '';
		  'BACKUP';
		  'The BACKUP button allows the user to move back to the correlation';
		  'coefficient setup interface to re-calculate this value if desired.'};

    %--------------
    case 'sf_hlp',
	he.str = { ...
          'This interface is used to set initial parameter values for';
		  'the SMAC AUTO FIT codes.';
		  '';
	   	  'FREQUENCY RANGE';
		  'For the automated SMAC algorithms to be executed, the user';
		  'must specify a set of frequency and damping ranges over which';
		  'the codes must search to converge on the true root. To specify';
		  'the search frequency bandwidth (BANDWIDTH for INDIVIDUAL ';
		  'FREQUENCY FIT), the user specifies a percentage of the peak';
		  'frequencies determined from the correlation calculations. A';
		  'bandwidth of 0.5 to 3 percent has been used successfully in';
		  'determining frequencies. The default is 1 percent (0.5% below';
		  'the peak frequency and 0.5% above the peak frequency).  You may';
          'also specify the number of frequencies over this range to use in';
		  'this fit.';
		  '';
		  'DAMPING RANGE';
		  'In the LOW % VALUE and HIGH % VALUE edit boxes, the user must';
		  'specify the damping range over which the codes must search to';
		  'converge on an appropriate damping value. The default values';
		  'are LOW % = 0.25 % and HIGH % = 5 %.  The user may also specify';
		  'the number of damping values to use over this range in the fit.';
		  '';
		  'PLOT FRFs';
		  'This pushbutton enables the user to view all the existing';
		  'FRF data that has been loaded into SMAC. Any number of data'; 
		  'channels can be overlayed to view response characteristics.';
		  '';
		  'EXECUTE AUTO SMAC';
		  'This pushbutton executes the automated SMAC algorithms in ';
		  'search of the true roots of the experimental data given the';
		  'specified frequency and damping ranges.';
		  '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.'; 
		  '';
		  'BACKUP';
		  'The BACKUP button allows the user to move back to the ';
		  'correlation coefficient calculation interface to re-calculate';
		  'peak and correlation values if desired.'};

    case 'sh_hlp',
	he.str = { ...
          'The results of the SMAC AUTO FIT will be displayed in the';
		  'listbox on the SMAC Synthesis and Mode Shape Generator';
		  'interface. The root (frequency and damping) that was ';
		  'calculated automatically from SMAC along with its correlation';
		  'value are listed with the original frequency and correlation.';
		  'The roots that could not be fit with the SMAC AUTO FIT ';
		  'algorithms return zeros for the frequency, damping and the';
		  'root correlation values. These can easily be removed from the';
		  'list table by using the pushbutton CONDENSE ROOTS. Also, at';
		  'any time the user can select a root from the list and delete';
		  'it using the DELETE button. A warning box comes up asking the';
		  'user to confirm before the roots are finally removed from';
		  'the list.';
		  '';
		  'The user can select all the roots (highlight list of roots)';
		  'by using the pushbutton ALL. A subset of roots may also be';
		  'selected by using the mouse and the <cntrl> key to highlight';
		  'only the roots of interest. To check the quality of the fit';
		  'the user then can do a SYNTHESIS.  The user first selects the';
          'type of syntheses to perform, then selects the SYNTHESIZE';
          'pushbutton to execute the synthesis.  This is important to help';
          'determine the quality of the fit.  This quality check should';
          'provide a determination whether all the modes of interest in';
		  'a particular frequency band have been properly identified.';
		  'NOTE: Only the roots that are selected in the list will be';
		  'included in the synthesis.';
		  'These quality checks (SYNTHESIS TYPES) include:';
		  '';
		  'NMIF SYNTHESIS';
		  'Normal Mode Indicator Function (NMIF) tends to be a very';
		  'sensitive comparison for both well excited and weakly excited';
		  'modes.';
		  '';
		  'CMIF SYNTHESIS';
		  'Complex Mode Indicator Function (CMIF) is a good tool for ';
		  'Finding multiple roots for multi-reference data and for';
		  'examining the strength of the modes excited in the data.';
		  'However, the CMIF can obscure the effects of the';
		  'weakly excited modes in the data.';
		  '';
		  'MMIF SYNTHESIS';
		  'Multivariate Mode Indicator Function (CMIF) is a good tool for ';
		  'finding multiple roots in consistent linear data sets.  It can';
		  'be confused by small frequency shifts from reference to reference';
		  'for the same mode in inconsistent data sets';
		  '';
		  'FRF SYNTHESIS';
		  'The typical method of checking whether all the modes of in a';
		  'frequency band have been properly identified is by comparisons';
		  'of the synthesized and actual FRF data. When REAL MODE has';
		  'been selected, only the imaginary part of the FRF is synthe-';
		  'sized and compared. When COMPLEX MODE has been selected, the';
		  'magnitudes of the FRF is compared.'; 
		  '';
		  'RESIDUALS';
		  'The residual terms toggle specifies whether to include residual';
          'terms in the synthesis.  The RESIDUALS button next to the toggle';
          'opens a form where you can configure each of the residual';
          'options that are available.  Residual terms account for out-of-';
          'band modes.';
          '';
		  '';
		  'After observing the synthesis of the current set of roots in';
		  'the list to the experimental data, the user can determine if';
		  'additional roots need to be identified. The ADD pushbutton brings';
		  'up a new form that allows the user to manually fit a root';
		  '';
		  'CREATE MODE SHAPES';
		  'This pushbutton executes code to calculate the mode shapes';
		  'and exports those mode shapes (in local coordinates) into an';
		  'ASH FILE FORMAT. The user is prompted through a browser';
		  'for the desired filename for the mode shape file.  The Residual';
		  'Terms toggle setting for the Synthesis determines whether these';
		  'residual terms will also be used in the mode shape generation.';
		  '';
		  'PLOT SHAPES';
		  'This pushbutton uses IMAT to let the user display and animate';
		  'the fit modes shapes. The first time this is selected, the user';
		  'will be prompted for a Univeral file containing the FEM geometry';
		  'on which to display the mode shapes.';
		  '';
		  'PLOT MAC';
		  'This pushbutton displays a Mode Assurance Criteria (MAC) bar';
          'graph.  The MAC provides an assessment of the linear indepen-';
          'dence of the mode shapes.  Diagonals are by definition 100%, so';
          'the off-diagonal terms are of the greatest importance.  Low off-';
          'diagonal terms indicate linear independent mode shapes.';
		  '';
		  'RESYNTHESIZE';
		  'This pushbutton opens a new form where you can resynthesize FRF';
          'and the derived quantities such as NMIF, CMIF, and MMIF from the';
          'mode shapes and residual information.  This is useful for valid-';
          'ating the mode shapes once they have been fit.';
		  '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.'; 
		  '';
		  'BACKUP';
		  'This pushbutton will send the user back to the SMAC AUTOFIT';
		  'interface and allow the user to re-calculate the data fits with';
          'different parameters.'};


    case 'mfit_hlp',

	he.str = { ...
          'The SMAC Manual Add Root GUI allows you to manually fit a';
		  'mode based on correlation coefficient.  This is an iterative';
		  'process, where you will need to alternate between a Frequency';
		  'Fit and a Damping Fit.  The idea is to zoom in on narrower';
		  'and narrower frequency and damping ranges until a parabolic';
		  'fit (shown as a dashed red line) matches very closely to the';
		  'synthesized data (shown as a solid blue line).';
		  '';
		  'When the form is initially displayed, the correlation coef-';
		  'ficients for the full frequency range are displayed.  You can';
		  'use the Select buttons or the edit boxes to select a new frequency';
		  'or damping range.  If your range is too narrow, use the Zooom';
		  'Out buttons to increase the range.  The Best Fit parameters are';
		  'shown in the top right of the form.  These are calculated from';
		  'the peak of the parabolic fit for frequency and damping, res-';
		  'pectively.';
		  '';
		  'The Add Root algorithm works by taking two separate fits of';
		  'the data.  For the Frequency Fit, the algorithm uses the Best Fit';
		  'damping and iterates over a range of frequencies specified by';
		  'the frequency fit.  It synthesizes FRF and computes correlation';
		  'coefficients.  The Damping Fit uses the Best Fit frequency and';
		  'iterates over a range of damping values specified on the form.';
		  'Each time a range changes, both fits are performed again and';
		  'the results displayed.';
		  '';
		  'Once you are satisfied, use the Add Root button to add this';
		  'root to the list.'};

    case 'res_hlp',
	he.str = { ...
          'The RESYNTHESIS capabilities in this function allow you to';
          'resynthesize Normal Mode Indicator Functions (NMIF), Complex';
          'Mode Indicator Functions (CMIF), Multivariate Mode Indicator';
          'Functions (MMIF), and Frequency Response Functions (FRF) from';
          'previously generated mode shapes.  If any residual effects were';
          'used in the original FRF synthesis when generating the mode';
          'shapes, these will be available during the resynthesis provided';
          'the residual information was saved into the ASH file when they';
          'were originally synthesized.';
          '';
          'You can selectively turn on and off any of the residual effects';
          'used in the original synthesis.  You can also control which';
          'modes are used by selecting the ones you wish to use from the';
          'mode shape listbox.';
          '';
          'The Synthesis Options area lets you choose which synthesis';
          'comparisons to make.  The resynthesized FRF will be overlaid on';
          'the experimental FRF you originally read in.  Once you have set';
          'up your resynthesis options, you can see the results by clicking';
          'on the PLOT button.  If you want to close the plot figures, you';
          'can either close them individually or click on the CLOSE PLOTS';
          'button.';
          '';
		  'SAVE';
		  'This pushbutton enables the user to save the existing SMAC';
		  'configuration to a .mat file for later processing.  This uses';
          'the same SMAC structure that is used during the curve fit, but';
          'much of the information is left blank.  As a result, you cannot';
          'use this MAT file for a SMAC curve fit.';
		  '';
		  'BACKUP';
		  'This pushbutton will send the user back to the SMAC SYNTHESIS';
		  'interface and allow the user to re-calculate the data fits with';
          'different parameters.  This button is only available if you';
          'arrived at this form from the SMAC SYNTHESIS form.'};

end

set(help_hand,'String',he.str);

% Resize the help window and the figure so the text fits

ext=get(help_hand,'Extent');
pos=get(help_hand,'Position');
pos(3)=ext(3)*1.2;
set(help_hand,'Position',pos);
pof=get(h0,'Position');
pof(3)=pos(3)+10;
set(h0,'Position',pof);

% Set the figure and control untils to normalized

set(help_hand,'Units','normalized');
set(h1,'Units','normalized');

set(h0,'Units','normalized');

% Now put the figure in the upper left corner of the screen and make it 65%
% of the screen height

pof=get(h0,'Position');
pof([1 2 4])=[0.01 0.96-0.65 0.65];
set(h0,'Position',pof);
