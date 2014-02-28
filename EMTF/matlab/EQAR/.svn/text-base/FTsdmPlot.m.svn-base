%  This is a driver script to make time-frequency section
%   plots of things that can be derived from the SDM
%   First step: construct a mat file that contains all SDM
%   results from processing a sequence of time segments
%   with multmtrn;   this is done with SDMmerge_40hz.m
%   Here we are processesing 40 hz data in two hour segments
LOAD = 0;
SAVEplots = 1;
MKPLOTS = ones(6,1);

if LOAD
   %  First load the mat file ... don't need to do this if
   %   called immediately after sdmMerge40Hz
   load S0_261_280_40Hz.mat
end

% construct time array for plotting
%   DAYS and HOUR are arrays stored in matfile giving
%    start time for each time segment
nFiles = length(SDMS)
for k = 1:nFiles
   Day(k) = DAYS(k)+HOUR(k)/24;
end

% set "good" indicator to 0 for segments with too few data
%    frac defines fraction of data that must be present to plot
frac = .5;
NFS = [];
for k = 1:nFiles
   if good(k)
      NFS = [NFS SDMS{k}.nf(1)];
   else
      NFS = [NFS 0];
   end
end
NFmin = max(NFS)*frac;
good(NFS<NFmin) = 0;

%   Period for each band in the SDM files .... we are assuming
%   throughout that all time segments are processed with the same
%   time windows and period bands
k1 = min(find(good));
Period = SDMS{k1}.T;

PLOTS = {'SDM Eigenvalues',...
	'Magnetic Field SNR',...
	'Electric Field SNR',...
	'Canonical Coherences: PKD vs. SAO',...
	'Canonical Coherences: H vs. E', ...
	'Average Modes'};

ROOTS = {'EVAL','MSNR','ESNR','CCPS','CCHE','MODE'};

nPLOTS = length(PLOTS);
for kPlot = 1:6
   if MKPLOTS(kPlot)
      plotType = PLOTS{kPlot};
      FTsdmPlotCases
      multiPlot(D,Day,Period,OPTIONS);
      if SAVEplots
         cfile = [ROOTS{kPlot} '_' num2str(f_day) '-' num2str(l_day) '.eps'];
         eval(['print -depsc ' cfile]);
      end
   end
end
