function s=smac_create_data_structure
% SMAC_CREATE_DATA_STRUCTURE  Create default SMAC data structure
%
% SS=SMAC_CREATE_DATA_STRUCTURE
%
% SMAC_CREATE_DATA_STRUCTURE creates a default SMAC structure.  This file
% also documents the structure field contents.

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
%  22-May-2004 / ATA Engineering / Dan Hensley
%    o Initial creation
%
%  03-Jun-2004 / ATA Engineering / Dan Hensley
%    o Added some entries
%
%  15-Jun-2004 / ATA Engineering / Dan Hensley
%    o Added frequency and damping convergence tolerance fields to the .fit
%      structure
%
%  06-Aug-2004 / ATA Engineering / Dan Hensley
%    o Change .pinv to cell array, 1 for each reference
%
%  03-Sep-2004 / ATA Engineering / Dan Hensley
%    o Add field to .corr for counting # repeated roots
%
%  01-Nov-2004 / ATA Engineering / Dan Hensley
%    o Split .mif field into .exp (experimental) and .ana (analytical)
%
%  10-Nov-2004 / ATA Engineering / Dan Hensley
%    o Add field .fit.corr_ref
%
%  03-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add fields to .residuals
%    o New version:  4
%
%  13-Jun-2005 / ATA Engineering / Dan Hensley
%    o Add rootsel field to .shape to remember selected roots
%    o New version:  5
%
%  28-Jun-2005 / ATA Engineering / Dan Hensley
%    o Remove mode inertance and compliance
%    o New version:  6
%
%  29-Jun-2005 / ATA Engineering / Dan Hensley
%    o Update structure for frequency range for low and high mode residuals
%    o New version:  7
%
%  31-Aug-2005 / ATA Engineering / Dan Hensley
%    o Add column to ss.corr.corr, and new field 'refpick' to .shape
%    o New version:  8
%
%  06-Mar-2006 / ATA Engineering / Dan Hensley
%    o Make default frequency convergence 0.05% (was 0.1%)
%    o Add new .shape.psi field and bump version number to 9
%
%==========================================================================

s=struct( ...
	'filename',     '', ...                % Filename from which data was imported
	'fe',           {imat_fn([])}, ...     % Experimental FRF (one column per reference)
	'fa',           {imat_fn([])}, ...     % Analytical FRF (one column per reference)
...
	'ref_coords',	{imat_ctrace([])}, ... % Reference coordinates
...
	'freqrange',	[], ...                % Pseudo-inverse frequency range
...                                        %   2x2 matrix, column 1 is frequency and column 2 is index
	'freqrangecc',	[], ...                % Correlation coefficient frequency range
...                                        %   2x2 matrix, column 1 is frequency and column 2 is index
	'svd',          [], ...                % NOT USED.  Structure containing SVD information
	'pinv',         [], ...                % Pseudo-inverse of experimental FRF matrix (third dimension, 1 per ref)
...
	'mif',          [], ...                % Structure containing MIF functions from all data
	'corr',         [], ...                % Structure containing correlation information
	'fit',          [], ...                % Structure containing curve-fit information
...
	'done',         [], ...                % Structure containing information for how far we've gone
...
	'realcomplex',	1, ...                 % 1=real normal modes, 2=complex modes
	'residuals',	[], ...                % Residuals to use in synthesis
...
	'shape',        [], ...                % Mode shape structure
...
	'version',	9);                		   % Version number of this structure

%-----------------------------------

s.corr=struct( ...
	'corrind',	[], ...             % Values selected from corr
	'corr',		[], ...             % Correlation coefficients (column 1=max value, column 2=ref giving max value)
	'corr_ref',	[], ...             % Correlation coefficients for each reference (column 1-n=corr coef for each reference)
	'nroots',	[], ...             % Number of roots for each selected index
	'zeta',		0.02, ...           % Damping ratio for initial correlation coefficient
	'nl',		20 );               % Number of spectral lines to use

s.fit=struct( ...
	'zeta',		[0.0025 0.05], ...	% 1x2 structure of damping range
	'zetapts',	21, ...             % Number of damping values in the range to use in the auto fit
	'zetaconv',	0.02, ...           % Damping convergence fraction (dcon in smac_autofit)
	'freqp',	1, ...              % Frequency range for autofit
	'freqpts',	11, ...             % Number of frequencies to use in the auto fit
	'freqconv',	0.0005, ...         % Frequency convergence fraction (fcon in smac_autofit)
	'rootlistorig',	[], ...         % nx3 matrix of original roots (freq, damp, corr)
	'rootlist',	[], ...             % nx3 matrix of roots (freq, damp, corr)
	'corr_ref',	[] );               % Matrix of correlation values for each root by reference

s.mif=struct( ...
	'exp',		[], ...             % MIFs from experimental data
	'ana',		[] );               % MIFs from analytical fits

s.mif.exp=struct( ...
	'nmif',		{imat_fn([])}, ...  % NMIF from all data
	'mmif',		{imat_fn([])}, ...  % MMIF from all data
	'cmif',		{imat_fn([])}, ...  % CMIF from all data
	'reftracking',  false);         % Use reference tracking option when evaluating CMIF
s.mif.ana=s.mif.exp;

s.done=struct( ...
	'pinv',		false, ...          % Calculated pseudo-inverse?
	'corr',		false, ...          % Calculated correlation coefficients?
	'selroots',	false, ...          % Selected initial roots?
	'autofit',	false );            % Performed auto SMAC curve-fit?

%-------

s.residuals=struct( ...
	'use',              false, ...  % Flag for whether residuals are used in the synthesis
	'lowmode',          [], ...     % Add an artificial low frequency mode
	'highmode',         [], ...     % Add an artificial high frequency mode
...	'modeinertance',	[], ...     % NOT USED: Add the inertance term around the mode
...	'modecompliance',	[], ...     % Add the compliance term around the mode
	'inertance',        [], ...     % Add inertance over a specified frequency band
	'compliance',       [] );       % Add compliance over a specified frequency band

s.residuals.lowmode=struct( ...
	'active',	false, ...          % Flag for whether this residual should be used
	'freq',	    [], ...             % Frequency at which the low mode is placed
	'damp',		[], ...             % Damping ratio for the low mode
	'frange',	[], ...             % Frequency range over which to fit the low mode
	'residual',	[] );               % Calculated residual value

s.residuals.highmode=struct( ...
	'active',	false, ...          % Flag for whether this residual should be used
	'freq',	    [], ...             % Frequency at which the high mode is placed
	'damp',		[], ...             % Damping ratio for the high mode
	'frange',	[], ...             % Frequency range over which to fit the high mode
	'residual',	[] );               % Calculated residual value

% s.residuals.modeinertance=struct( ...
% 	'active',	false, ...          % Flag for whether this residual should be used
% 	'residual',	[] );               % Calculated residual value
% 
% s.residuals.modecompliance=struct( ...
% 	'active',	false, ...          % Flag for whether this residual should be used
% 	'residual',	[] );               % Calculated residual value

s.residuals.inertance=struct( ...
	'active',	false, ...          % Flag for whether this residual should be used
	'freq',		[], ...             % Frequency range over which these terms are calculated
	'residual',	[] );               % Calculated residual value

s.residuals.compliance=struct( ...
	'active',	false, ...          % Flag for whether this residual should be used
	'freq',		[], ...             % Frequency range over which these terms are calculated
	'residual',	[] );               % Calculated residual value

%-------

s.shape=struct( ...
	'residues',	[], ...             % Store the residues (shape fitting only)
	'psi',	    [], ...             % Reciprocal modal vectors (psi) from the fit
	'shape',	{imat_shp([])}, ... % Curve-fit mode shapes
	'refpick',	[], ...				% Reference from which the final shape was selected
	'rootsel',	[], ...             % List of roots selected in the synth form
	'fem',		[] );               % FEM connectivity
