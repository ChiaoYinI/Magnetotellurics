function [FTres] = sdmFilt(FT,TFft);
% Usage [FTres] = sdmFilt(FT,TFft);
%  applies filter in structure TFft to FT, yielding 
%  residuals in FTres (if TFft.Resid == 1 ... otherwise
%   predicted)

Kin = TFft.Kall(TFft.KallIn);
Kout = TFft.Kall(TFft.KallOut);
nchOut = length(Kout);
tf = FT.SysTF.TF(Kout,:);
SysTF = struct('TF', FT.SysTF.TF(Kout,:),...
		'freq',FT.SysTF.freq, ...
		'PhysU',FT.SysTF.PhysU);
FTres = struct('t',FT.t,'dt',FT.dt,'f',FT.f,'ch_id',FT.ch_id(Kout,:),...
	'sta',FT.sta(Kout,:),'nch',nchOut,'SysTF',SysTF, ...
        'amp',0);

[nch,nt,nf] = size(FT.data);
FT.data = reshape(FT.data,[nch,nt*nf]);
FTres.data = zeros(nchOut,nt*nf)+i*zeros(nchOut,nt*nf);
nTF = length(TFft.freq);
for ib = 1:nTF
   I = find(TFft.indfData==ib);
   if(~isempty(I))
      i1 = (min(I)-1)*nt+1;
      i2 = max(I)*nt;
      temp = squeeze(TFft.A(:,:,ib))*FT.data(Kin,i1:i2);
      if (TFft.Resid == 1)
         FTres.data(:,i1:i2) = FT.data(Kout,i1:i2) - temp;
      else
         FTres.data(:,i1:i2) = temp;
      end
   end
end
FTres.data = reshape(FTres.data,[nchOut,nt,nf]);
