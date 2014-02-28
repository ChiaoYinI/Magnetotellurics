%   define PKDSAO root data directory 
%      (use environment vbl $PKDSAOdata)
PKDSAOdata = setDataDir('PKDSAOdata');
ADDyear = 0;
%  for new file names set ADDyear = 0; ADDyear = 1 only for
%   old file naming conventions (where *.z** files had names
%   like PKD_159_2003.zmm)
% ADDyear = 1;

%   define time window to plot
YearOne = 2002;
YearEnd = 2002;
dayOne = 1;
dayEnd = 365;

%  set SAVEplots = 1 to save figures as eps files
SAVEplots = 0;
plottype='eps'; %choose from eps, jpeg, eps will select epsc automatically
EQLINE = 0;

%  set output directory for plot files
PlotDir =  [PKDSAOdata num2str(YearOne),'/RHOplots/'];
Fdir = [PKDSAOdata num2str(YearOne),'/f/'];
SITES =  {'PKD','SAO'};

%   then define what is to be plotted:
iSite = 1;       %  plot for which site?  1 = PKD, 2 = SAO
plotBoth = 1;    %  set to 1 to plot both 100 and 200 m dipoles
		 %   (only matters for PKD)
RR =  2;         %  RR = 0,1,2 for single site (*.zss files),
		 %  remote reference (*.zrr files),
		 %  multiple station (*.zmm files) , respectively
BAND = 1;        %  BAND = 1, 40 for 1 hz and 40 hz cases
SUFFIX = '';    %  File naming suffix: 
		 %   For example, PKD_010nr.zss could refer to data processed
		 %   (single site) with robust options turned off
dayStep = 1;     %  1 for 40 hz  ... could be greater for 1 hz
		 % E.G.: dayStep = 10 would allow plotting of 1hz data processed
	 	 %  in 10 day segments
hourStep = 2;    %  only need this defined for 40 hz

%   set subMed = 1 to subtract median over time for each frequency
%    band (i.e., plot deviations from time average rho, phi, etc.)
%       subMed = 0 plots rho phi without median subtracted;
%       could add additional option to read in a long term avg to subtract
subMed = 1;

%  don't plot anything for which fractional error in complex 
%  imepedance exceeds this value
FracErrCutOff = .2;

%   Following are computed, don't need to edit
%   nHour is number of 40 hz segments in 1 day (e.g., 12 for 2 hr)
%   TimeIncDay is length of segment, in days
if BAND == 40
  nHour = 24/hourStep;  %  need this in both cases; could define directly for 1h
  TimeIncDay = 1/nHour;     %  need this in both cases; could define directly for 1hz
else
  TimeIncDay = dayStep;
end
