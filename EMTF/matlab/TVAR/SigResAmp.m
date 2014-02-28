function [FTsigAmp,FTresAmp] = SigResAmp(FCfiles,OPTIONS,TF);
%Usage: [FTsigAmp,FTresAmp] = SigResAmp(FCfiles,OPTIONS,TF);
%                          SigResAmp(FCfiles,OPTIONS,TF);
%       [FTsigAmp,FTresAmp] = SigResAmp(FCfiles,OPTIONS);
%                           SigResAmp(FCfiles,OPTIONS);
%  plots residual and/or signal amplitude for selected channels
%  using FC files given in list FCfiles
%
%  An optional TF data structure can be provided as input;
%   if this argument is missing the FCfiles list is used to
%   define an sdm file name (*.S0) and default rules are followed
%   to define the prediction TF (variants might be specified through
%   additional fields in the OPTIONS data structure, but nothing
%   like this is coded at present).
%  OPTIONS has the following fields, used to define which channels 
%   are used for prediction/residuals, band and time averaging,
%   and some plotting details
%   OPTIONS.dir         : root processing output directory 
%				(contains FC, TS, etc.)
%   OPTIONS.Kin, Kout   : channels used for prediction, residuals
%                        these are numbers, following the channel
%                         list constructed from each FC file in
%                         the order given.  For standard case these
%                         are PKD: Hx, Hy, Hz, Ex1, Ey1, Ex2, Ey2
%                             SAO: Hx, Hy, Hz, Ex, Ey
%                         NOTE: Kin, Kout set in optional input TF
%                              file overides OPTIONS for TF/residual plots
%                          Kout is still used to determine plotted signal
%                            channels, and corresponds to channel order if FT
%                            data structure
%  OPTIONS.navg		: number of time segments to average
%  OPTIONS.id		: decimation level
%  OPTIONS.iband	: frequency bands for frequency averaging
%  OPTIONS.plotAmp	: = 1 to plot signal amplitudes
%  OPTIONS.plotRes	: = 1 to plot residual amplitudes
%  OPTIONS.title, NormalizeSig, NormalizeRes, clSig, clRes, linT ... 
%		  control plotting; Idea is that one could omit plotting in
%			  this routine, just return the time/freq
%			  averaged powers and plot these with a
%                         different routine
%   NOTE:  FTsigPwr and FTresPwr are optional outputs, data structures of
%      FT type containing averaged signal and residual power  ... omit
%      these for standard plotting when there is no need to save
%      residual or signal power for other purposes.


nsta = length(FCfiles);
for k = 1:nsta
   ir = max(find(char(FCfiles{k})=='.'));
   root = FCfiles{k}(1:ir-1);
   ir = min(find(root=='_'));
   csta{k} = root(1:ir-1);
   FileNamesFC{k} = [OPTIONS.dir '/FC/' FCfiles{k}];
   FileNamesSysTF{k} = [OPTIONS.dir '/sysTF/' root '.stf'];
end
if nargin == 2
   %   generate TF from default SDM derived from these FC files
   ir = min(find(root=='_'));
   root = root(ir+1:end);
   FileNameSDM = [OPTIONS.dir '/SDM/' root '.S0'];
   [TF] = TFsdmSetX(FileNameSDM,OPTIONS.Kout,OPTIONS.Kin)
   if ~isstruct(TF)
      fprintf(1,'%s\n','Error reading TF file')
      FTsigAmp = [];
      FTresAmp = [];
      return
   end
end
   
%  limit frequencies to those needed for requested frequency bands
if1 = min(min(OPTIONS.iband));
if2 = max(max(OPTIONS.iband));
ifreq = [if1:if2];
%  band definitions are shifted when if1 < 1
iband = OPTIONS.iband - if1+1;
navg = OPTIONS.navg;

%   get Freq/Time array of FCs 
[FT] = FTsetup(FileNamesFC,OPTIONS.id,ifreq,csta);
if ~isstruct(FT)
   %  FTsetup failed (some FC files don't exist/empty)
   % (empty data structure returned on error)
   fprintf(1,'%s\n','Error reading FC files')
   FTsigAmp = [];
   FTresAmp = [];
   return
end 

%   get system transfer function file; add to FT data structure
[sysTF] = sysTFsetup(FileNamesSysTF);
if ~isstruct(sysTF)
   %  sysTFsetup failed (system parameter files don't exist/empty)
   % (empty data structure returned on error)
   fprintf(1,'%s\n','Error reading sysTF files')
   FTsigAmp = [];
   FTresAmp = [];
   return
end 

FT.SysTF = sysTF;

%   check channel compatability, create modified TF structure
[TFft] = TFcompatFT(FT,TF);

%   filter FT array, compute residuals
[FTres] = sdmFilt(FT,TFft);

%   Amplitudes (not power!) of signal, residuals 
[FTsigAmp] = FTavg(FT,OPTIONS.Kout,navg,iband);
chOut = [1:length(TF.Kout)];
[FTresAmp] = FTavg(FTres,chOut,navg,iband);

if (OPTIONS.plotAmp | OPTIONS.plotRes)
  %  plot amplitudes
  SigResPlot(FTsigAmp,FTresAmp,OPTIONS);
end
