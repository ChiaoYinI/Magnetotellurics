function [CC] = compCC(Sdms,ind1,ind2,locErr)
%   computes canonical covariances and related quantities,
%   for two groups of channels, returns in data structure CC
%   ind1, ind2 give indices used to define two channel groupings,
%   Sdms is SDM data structure
%  u1 and u2 are in SNR coordinates

n1 = length(ind1);
n2 = length(ind2);
nc = min(n1,n2);
T = Sdms.T';
nbt = length(T);

ev1 = zeros(n1,nbt);
ev2 = zeros(n2,nbt);
ccov = zeros(nc,nbt);
ccor = zeros(nc,nbt);

if nargin < 4
  locErr = 0;
end

if locErr
   %   modifications for local error correlation estiamtion: not
   %   tested or used, maybe a bad idea
   %    try replacing this with something that accounts roughly for
   %    statistical significance
   nSigMin = 2;
   nSigMax = 4;
   ccovMin = ones(1,nbt);
   E11 = zeros(n1,n1,nbt);
   E22 = zeros(n2,n2,nbt);
   C11 = zeros(n1,n1,nbt);
   C22 = zeros(n2,n2,nbt);
   U1 = zeros(n1,n1,nbt)+i*zeros(n1,n1,nbt);
   U2 = zeros(n2,n2,nbt)+i*zeros(n2,n2,nbt);
end

for ib = 1:nbt
   S = squeeze(Sdms.S(:,:,ib));
   var = Sdms.var(:,ib);
   sig = 1./sqrt(var);
   S11 = S(ind1,ind1);
   S22 = S(ind2,ind2); 
   S12 = S(ind1,ind2);
   temp = real(eig((S12'/S11)*S12,S22));
   temp = sort(temp);
   ccor(:,ib) = temp(n2:-1:n2-nc+1);
   S11 = diag(sig(ind1))*S11*diag(sig(ind1));
   S22 = diag(sig(ind2))*S22*diag(sig(ind2));
   S12 = diag(sig(ind1))*S12*diag(sig(ind2));
   temp = real(eig(S11));
   temp = sort(temp);
   ev1(:,ib) = temp(n1:-1:1);
   temp = real(eig(S22));
   temp = sort(temp);
   ev2(:,ib) = temp(n2:-1:1);
   [u1,ev,u2] = svd(S12,0);
   [temp,ind] = sort(diag(ev));
   ind = ind(end:-1:end-nc+1);
   ccov(:,ib) = temp(end:-1:end-nc+1);
   u1 = u1(:,ind);
   u2 = u2(:,ind);

   if locErr
      nSig = sum(ccov(:,ib)>ccovMin(ib))
      nSig = max(nSig,nSigMin);
      nSig = min(nSig,nSigMax);
      u1 = u1(:,end:-1:end-nSig+1);
      u2 = u2(:,end:-1:end-nSig+1);
      E11(:,:,ib) = S11 - u1*diag(ccov(1:nSig,ib))*u1';
      E22(:,:,ib) = S22 - u2*diag(ccov(1:nSig,ib))*u2';
      ErrInv = 1./sqrt(diag(squeeze(E11(:,:,ib))));
      C11(:,:,ib) = diag(ErrInv)*squeeze(E11(:,:,ib))*diag(ErrInv);
      ErrInv = 1./sqrt(diag(squeeze(E22(:,:,ib))));
      C22(:,:,ib) = diag(ErrInv)*squeeze(E22(:,:,ib))*diag(ErrInv);
   end
   U1(:,:,ib) = u1;
   U2(:,:,ib) = u2;
end    % ib

if locErr
   CC = struct('T',T,'ccor',ccor,'ccov',ccov,'ev1',ev1,'ev2',ev2, ...
	'u1',U1,'u2',U2,'C11',C11,'C22',C22,'E11',E11,'E22',E22);
else
   CC = struct('T',T,'ccor',ccor,'ccov',ccov,'ev1',ev1,'ev2',ev2, ...
	'u1',U1,'u2',U2);
end
