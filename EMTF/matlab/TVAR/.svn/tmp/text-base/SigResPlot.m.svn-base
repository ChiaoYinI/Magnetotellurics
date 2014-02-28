function SigResPlot(FTsigAmp,FTresAmp,OPTIONS);
%Usage: SigResPlot(FTsigAmp,FTresAmp,OPTIONS);
%  plots residual and/or signal power in FT structures created by SigResAmp.m
%   OPTIONS has these fields:  
%       plotAmp, plotRes = 0/1 to plot or not plot signal residuals
%       NormalizeSig, NormalizeRes 
%       title
%       clSig(2,nch)  :  coloraxis limits for each channel (signal)
%       clRes(2,nch)  :      and residuals  (these are used only for 
%				unnormalized plots)
%       linT	      :  linear or logarithmic frequency axis

[nch,nT,nF] = size(FTsigAmp.data);
Ind = cumsum([1 round(diff(FTsigAmp.t)/FTsigAmp.dt)]);
nTmax = max(Ind);
%nTmax = max(nT,fix((FTsigAmp.t(end)-FTsigAmp.t(1))/FTsigAmp.dt+1));
FILL = nTmax > nT;
D = zeros(nTmax,nF,nch);
if FILL
    % There are missing segments ... need to leave some sections of plot blank
    D = D./D;
    t = FTsigAmp.t(1) + [0:max(Ind)-1]*FTsigAmp.dt;
else
    t = FTsigAmp.t;
end

if OPTIONS.plotAmp
   %  Plot amplitude 
   linT = OPTIONS.linT;
   for k = 1:nch
      if OPTIONS.NormalizeSig
         cAxLab = 'dB (Amplitude)';
         tit2 = ' : Variation in Signal Amplitude'
         cl = OPTIONS.clSigNorm;
         logSpec = ones(nT,1)*log10(FTsigAmp.avgAmp(k,:));
         if FILL
             D(Ind,:,k) = 10*(squeeze(log10(FTsigAmp.data(k,:,:)))-logSpec);
         else
             D(:,:,k) = 10*(squeeze(log10(FTsigAmp.data(k,:,:)))-logSpec);
         end
      else
         cAxLab = 'log_{10} Amp';
         cl = OPTIONS.clSig(:,k);
         if FILL
             D(Ind,:,k) = log10(FTsigAmp.data(k,:,:));
         else
             D(:,:,k) = log10(FTsigAmp.data(k,:,:));
         end
      end
      OPT{k} = struct('Caxis',cl,...
          'SubTitle',[FTsigAmp.sta(k,:) ':' FTsigAmp.ch_id(k,:)],...
          'Title','',...
          'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','Frequency (hz)', ...
          'ColorAxisLabel',cAxLab,...
          'DayLine',0);
      if(~linT) 
         OPT{k}.PeriodAxisLabel = 'log_{10} Frequency (hz)';
      end
   end
   OPT{1}.Title = [ OPTIONS.title tit2];

   multiPlot2(D,t,FTsigAmp.f,OPT,linT);
   if OPTIONS.EQLINE
      if OPTIONS.linT
         plotEQ
      else
         plotEQlog
      end
   end
end   %  if plotAmp

if OPTIONS.plotRes
   %  Plot  residuals
   for k = 1:nch
      if OPTIONS.NormalizeRes
         %  for residuals normalize by average signal amplitude
         cl = OPTIONS.clResNorm;
         logSpec = ones(nT,1)*log10(FTsigAmp.avgAmp(k,:));
         if FILL
             D(Ind,:,k) = 10*(squeeze(log10(FTresAmp.data(k,:,:)))-logSpec);
         else
             D(:,:,k) = 10*(squeeze(log10(FTresAmp.data(k,:,:)))-logSpec);
         end
      else
         cl = OPTIONS.clRes(:,k);
         if FILL
             D(Ind,:,k) = log10(FTresAmp.data(k,:,:));
         else
             D(:,:,k) = log10(FTresAmp.data(k,:,:));
         end
      end
      OPT{k} = struct('Caxis',cl,...
          'SubTitle',[FTresAmp.sta(k,:) ':' FTresAmp.ch_id(k,:)],...
          'Title','',...
          'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','Frequency (hz)', ...
          'ColorAxisLabel','log_{10} Amp',...
          'DayLine',0);
      if(~linT) 
         OPT{k}.PeriodAxisLabel = 'log_{10} Frequency (hz)';
      end
   end
   if OPTIONS.NormalizeRes
      OPT{1}.Title =  [OPTIONS.title ' : Residual Amplitude/Avg. Signal Amplitude']
      OPT{nch}.ColorAxisLabel = 'dB (Amplitude)'
   else
      OPT{nch}.Title = [OPTIONS.title ' : Residuals'];
   end
   linT = OPTIONS.linT;

   multiPlot2(D,t,FTresAmp.f,OPT,linT);
   if OPTIONS.EQLINE
      if OPTIONS.linT
         plotEQ
      else
         plotEQlog
      end
   end
end
