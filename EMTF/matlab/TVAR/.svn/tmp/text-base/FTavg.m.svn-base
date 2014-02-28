function [FTamp] = FTavg(FT,ChUse,nAvg,iband);
%  Usage: [FTamp] = FTavg(FT,ChUse,nAvg,iband);
%   given input complex frequency/time data structure FT
%   compute RMS amplitude averaged over frequeny bands iband
%   and over nAvg consecutive segments for channels numbers
%   listed in ChUse;  Overall RMS time average is added to 
%   FT structure, which contains only amplitudes
%     iband is optional ... default is single frequency bands
%     nAvg is optional  ... default is 1
%     ChUse is optional ... default is all channels
%   This does not presently allow for missing data !!

[nCH,nT,nF] = size(FT.data);
if nargin == 1
  ChUse = [1:nCH];
end
nCHavg = length(ChUse);
if nargin <= 2
   %  just compute abs value of each FC, return in F/T structure
   amp = abs(FT.data(ChUse,:,:));
   avgAmp = squeeze(sqrt(mean(amp.*amp,2)));
   FTamp =  struct('data',amp,'t',FT.t,'f',FT.f,'ch_id',...
	FT.ch_id(ChUse,:),'sta',FT.sta(ChUse,:), ...
	'nch',nCHavg,'amp',1,'avgAmp',avgAmp);
else
   %   time/frequency average
   if nargin == 3
      iband = [1:nF; 1:nF];
   end 
   nTavg = floor(nT/nAvg);
   nTuse = nTavg*nAvg;
   [dum,nFavg] = size(iband);
   amp = zeros(nCHavg,nTavg,nFavg);
   for ib = 1:nFavg
      nb = iband(2,ib)-iband(1,ib)+1;
      fAvg(ib) = mean(FT.f(iband(1,ib):iband(2,ib))); 
      for ich = 1:nCHavg
         temp = FT.data(ChUse(ich),1:nTuse,iband(1,ib):iband(2,ib));
   %      size(temp)
         temp = reshape(abs(temp).^2,[nAvg,nTavg,nb]);
         if(nb == 1)
            amp(ich,:,ib) = sum(temp,1)/nAvg;
         else
            if(nAvg == 1)
               amp(ich,:,ib) = sum(temp,3)/nb;
            else
               amp(ich,:,ib) = sum(squeeze(sum(temp,1)),2)/(nAvg*nb);
            end
         end
      end
   end
   avgAmp = sqrt(mean(amp,2));
   amp = sqrt(amp);
   temp = reshape(FT.t(1:nTuse),nAvg,nTavg);
   dt = FT.dt*nAvg;
   tAvg = mean(temp,1);
   FTamp  =  struct('data',amp,'t',tAvg,'dt',dt,'f',fAvg,'ch_id',...
	FT.ch_id(ChUse,:),'sta',FT.sta(ChUse,:), ...
	'nch',nCHavg,'amp',1,'avgAmp',avgAmp);
end
