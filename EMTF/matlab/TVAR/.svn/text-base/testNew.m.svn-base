path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/TVAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/EQAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/SDM',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/IN',path);
dir = '/home/server/pi/homes/egbert/PKDSAOg/data/2004/';
FileNamesFC = {['FC/PKD_270_00.f7'],['FC/SAO_270_00.f5']};
FileNamesSysTF = {[dir 'sysTF/PKD_270_00.stf'],[dir 'sysTF/SAO_270_00.stf']};
FileNameSDM = [dir 'SDM/270_00.S0'];
id = 1;
%  limit frequencies to those appropriate to decimation level?
ifreq = [5:112];
csta = {'PKD','SAO'};
%   get Freq/Time array of FCs
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
%   get system transfer function file; add to FT data structure
[sysTF] = sysTFsetup(FileNamesSysTF);
FT.SysTF = sysTF;
%   get interstation TF from SDM file
%  Magnetic residuals at Parkfield:
Kout = [1 2 3];
%  Use Hollister only for prediction (omit Hz)
Kin = [1:12];
[TF] = TFsdmSet(FileNameSDM,Kout,Kin);
%   check channel compatability   
[TFft] = TFcompatFT(FT,TF);
[FTres] = sdmFilt(FT,TFft);

%   Now some simple tests
%   Amplitudes  of FCs ... 
Normalize = 0;
%  only compute amplitudes for channels listed in Kout
[FTamp] = FTavg(FT,Kout);
[FTampRes] = FTavg(FTres);
%  Plot PKD mags
[dum,nT,nF] = size(FTamp.data);
D = zeros(nT,nF,3);
for k = 1:3
   if Normalize
      cl = [-.5,.5];
      logSpec = ones(nT,1)*log10(FTamp.avgAmp(k,:));
      D(:,:,k) = squeeze(log10(FTamp.data(k,:,:)))-logSpec;
   else
      cl = [-4,-1];
      D(:,:,k) = log10(FTamp.data(k,:,:));
   end
   OPTIONS{k} = struct('Caxis',cl,...
          'SubTitle',FT.ch_id(k,:),...
          'Title','',...
          'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','Frequency (hz)', ...
          'ColorAxisLabel','log_{10} Amp',...
          'DayLine',0);
end
OPTIONS{1}.Title = 'Magnetic Fields at PKD'
linT = 1;

multiPlot(D,FTamp.t,FTamp.f,OPTIONS,linT);

%  Plot PKD mag residuals
Normalize = 1;
[dum,nT,nF] = size(FTampRes.data);
D = zeros(nT,nF,3);
for k = 1:3
   if Normalize
      %  for residuals normalize by average signal amplitude
      cl = [-30,0];
      logSpec = ones(nT,1)*log10(FTamp.avgAmp(k,:));
      D(:,:,k) = 20*(squeeze(log10(FTampRes.data(k,:,:)))-logSpec);
   else
      cl = [-4,-1];
      D(:,:,k) = log10(FTampRes.data(k,:,:));
   end
   OPTIONS{k} = struct('Caxis',cl,...
          'SubTitle',FT.ch_id(k,:),...
          'Title','',...
          'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','Frequency (hz)', ...
          'ColorAxisLabel','log_{10} Amp',...
          'DayLine',0);
end
if Normalize
   OPTIONS{1}.Title = ....
   'Magnetic Field Residual Power at PKD : Normalized by Signal Power'
   OPTIONS{3}.ColorAxisLabel = 'dB (Power)'
else
   OPTIONS{3}.Title = 'Magnetic Field Residuals PKD'
end
linT = 1;

multiPlot(D,FTampRes.t,FTampRes.f,OPTIONS,linT);
