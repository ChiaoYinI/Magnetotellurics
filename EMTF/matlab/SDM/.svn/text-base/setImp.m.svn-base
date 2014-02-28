function [Zall] = setImp(Sdms)
% Usage:  [Zall] = setImp(Sdms);
%   extracts impedances from first 2 evecs of SDM

nbt = Sdms.Hd.nbt;
nt = Sdms.Hd.nt;
nsta = Sdms.Hd.nsta;
ih = Sdms.Hd.ih;
ie = [ih(2:end)-2 nt-1]; 
Zall = zeros(2,2,nsta,nbt);
for ib = 1:nbt
   sig = sqrt(Sdms.var(:,ib));
   u = diag(sig)*squeeze(Sdms.U(:,1:2,ib));
   for ista = 1:nsta
      Zall(:,:,ista,ib) = u(ie(ista):ie(ista)+1,:)/...
      		u(ih(ista):ih(ista)+1,:);
   end
end
