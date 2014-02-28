function [TFft] = TFcompatFT(FT,TF)
% Usage: [TFft] = TFcompatFT(FT,TF);
%  Constructs  mapping between data channels in FT and TF,
%   and outputs compatible TFft (including correct PhysU, order, etc.
%    + indices for data channels needed for prediction filter
%    and then constructs an actual TF (usable when all channels in KallIn are
%     available
%  NOTE: output structure is NOT the same as input!
%  Output structure: TF.Kall = list of data channels needed and available
%          in FT (in/out)
%    TF.KallIn/Out : elements within Kall to use as input/output
%    TF.V, TF.var TF : info for channels in Kall (same order)
%    TF.SystTF : system parameter structure for output data channels
%    TF.A  : TF matrix mapping input to output channels
%    TF.indfData :  indicies in TF array A to use for each frequency in FT

nKin = length(TF.Kin);
nKout = length(TF.Kout);
Kin = zeros(nKin,1);
Kout = zeros(nKout,1);

if length(TF.dim)==1
   dim = TF.dim;
else
   dim = max(TF.dim);
end

[dum,nCharD] = size(FT.ch_id);
[dum,nCharTF] = size(TF.ch_id);
nCharCH = min(nCharD,nCharTF);
[dum,nCharD] = size(FT.sta);
[dum,nCharTF] = size(TF.sta);
nCharSTA = min(nCharD,nCharTF);
for k = 1:nKin
   for j = 1:FT.nch
      if((FT.ch_id(j,1:nCharCH) == TF.ch_id(TF.Kin(k),1:nCharCH)) & ...
         (FT.sta(j,1:nCharSTA) == TF.sta(TF.Kin(k),1:nCharSTA)))
         Kin(k) = j;
      end
   end
end
for k = 1:nKout
   for j = 1:FT.nch
      if((FT.ch_id(j,1:nCharCH) == TF.ch_id(TF.Kout(k),1:nCharCH)) & ...
         (FT.sta(j,1:nCharSTA) == TF.sta(TF.Kout(k),1:nCharSTA)))
         Kout(k) = j;
      end
   end
end

if(min(Kin>0))
   kTFin = TF.Kin;
%  kTFin should = TF.Kin;
%  this is not necessarily true for kTFout
   kTFout = TF.Kout(Kout>0);
   Kout = Kout(Kout>0);

%  Kout points to output data channels, kTFout to corresponding TF channels;
%   same for inputs

%  make list of all needed data channels (all, input and output)
%   Result: Kall, KallIn, KallOut
   temp = zeros(FT.nch,1);
   temp(Kin) = 1; temp(Kout) = 1;
   Kall = find(temp);
   nchAll = length(Kall);
   KallIn = [];
   KallOut = [];
   KTFall = zeros(size(Kall));
   for k = 1:length(Kin)
     kk = find(Kall == Kin(k));
     kTFall(kk) = kTFin(k);
     KallIn = [KallIn kk];
   end
   for k = 1:length(Kout)
     kk = find(Kall == Kout(k));
     kTFall(kk) = kTFout(k);
     KallOut = [KallOut kk ];
   end

   V = TF.V(kTFall,1:dim,:);
   var = TF.var(kTFall,:);
   f0 = 1./TF.T';
%  if necessary first transform prediction TF
%   for consistency with data TF (physical units or not)
   if((FT.SysTF.PhysU == 0 )  & ( TF.PhysU == 1))
      %    convert TF to measurement units
      for ic = 1:nchAll
         tf = interp1(FT.SysTF.freq,FT.SysTF.TF(Kall(ic),:),f0);
         tf(isnan(tf)) = 0;
         for l = 1:dim
            V(ic,l,:) = (squeeze(V(ic,l,:))).'./tf;
         end
         var(ic,:) = var(ic,:)./(abs(tf).^2);
      end
   else
      if((FT.SysTF.PhysU == 1 )  & ( TF.PhysU == 0))
      %    convert TF to physical units
         for ic = 1:nchAllTF
            tf = interp1(FT.SysTF.freq,FT.SysTF.TF(Kall(ic),:),f0);
            tf(isnan(tf)) = 0;
            for l = 1:dim
               V(ic,l,:) = squeeze(V(ic,l,:)).*tf;
            end
            var(ic,:) = var(ic,:).*(abs(tf).^2);
         end
      end
   end
   
   if( TF.Resid == -1)
      tf = ones(size(FT.SysTF.TF(1:dim,:)));
      SysTF = struct('TF',tf,'freq',FT.SysTF.freq,'PhysU',FT.SysTF.PhysU);
   else
      tf = FT.SysTF.TF(Kout,:);
      SysTF = struct('TF',tf,'freq',FT.SysTF.freq,'PhysU',FT.SysTF.PhysU);
   end
   
   TFft = struct('V',V,'var',var,'Kall',Kall,'dim',TF.dim,...
	   'KallIn',KallIn,'KallOut',KallOut,'SysTF',SysTF,...
	   'freq',f0,'Resid',TF.Resid,'PhysU',FT.SysTF.PhysU,...
           'ch_id',FT.ch_id(Kall,:),'sta',FT.sta(Kall,:),...
	   'UseVout',{TF.UseVout});
   [A,indfData] = InterpTF_FT(TFft,FT);
   TFft.A = A; 
   TFft.indfData = indfData;
else
   TFft = [];
end

