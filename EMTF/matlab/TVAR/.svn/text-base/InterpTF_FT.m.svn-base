function [A,indfData] = InterpTF_FT(TF,FT);
% Sets up actual TF, indicies for nearest-neighbor interpolation of TF
%   to data FC frequencies
% Usage: [A,indfData] = InterpTF_FT(TF,FT)
%   slightly modified from InterTF ... this takes standard FT
%   structure as input, instead of TStemp;
%   also,   "HPf" cut off is not  implemented ... to add this,
%   a cut off frequency would have to be added to TF structure
%   Nearest neighbor interpolation is used
%   NOTE:  this will work with generalized input V produced by
%     CCsigNoiseVec... those vectors associated with incoherent
%     noise in the output channels are ignored;
%     ALSO: variable number of vectors in each freq band is,
%      supported 

[nch,dim,nf] = size(TF.V);
nIn = length(TF.KallIn);
nOut = length(TF.KallOut);
if(TF.Resid==-1)
   nOut = dim;
end
A = zeros(nOut,nIn,nf);
for k = 1:nf
   varInv = diag(1./TF.var(TF.KallIn,k));
   vIn = squeeze(TF.V(TF.KallIn,:,k));
   ndim = min(find(vIn(1,:)==0))-1;
   if(isempty(ndim)) ndim = dim; end
   vIn = vIn(:,1:ndim);
   useVout = find(TF.UseVout{k});
   vOut = squeeze(TF.V(TF.KallOut,useVout,k));
   nVout = length(useVout);
   if(nVout) < ndim
      vOut(:,nVout+1:ndim) = 0;
   end
   temp = varInv*vIn;
   if(TF.Resid >= 0)
      A(:,:,k) = vOut*inv(vIn'*temp)*temp';
   else
      A(:,:,k) = inv(vIn'*temp)*temp';
   end
end

%   Compute frequency indices of TF A for all frequencies
%  f is array of frequencies
fData = FT.f;
[fTF,I] = sort(TF.freq);
%  next two lines are to make sure that the sorted list of
%  frequencies are monotonic (might not be if two TF estiamtion bands
%   (e.g., from two different decimation levels) have the same center
%    frequency)
ii = find(diff(fTF) == 0);
fTF(ii+1) = fTF(ii+1) + fTF(1)*.001;
fTF = [fTF fTF(end)*5]; 

indfTF = [1:nf];
indfTF = indfTF(I);
indfTF = [indfTF 0 ];
%indfData = interp1(fTF,indfTF,fData(fData>TS.HPf),'nearest');
indfData = interp1(fTF,indfTF,fData,'nearest');
