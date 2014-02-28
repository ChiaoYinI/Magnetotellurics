function [Unorm] = impNorm(Sdms,Z,ivec,ib,rho_ref);

%  J is to rotate E field after transformation to H by admittance
J = [0 1 ; -1 0];
%  following 3 lines are to take out effect of reference rho 
%     scaling in u_pair.m
period = Sdms.T(ib);
escl = sqrt(period/(5*rho_ref));
J = J/escl;

nt = Sdms.Hd.nt;
nsta = Sdms.Hd.nsta;
ih = Sdms.Hd.ih;
ie = [ih(2:end)-2 nt-1];
var = Sdms.var(:,ib);
Unorm = diag(sqrt(var))*squeeze(Sdms.U(:,ivec,ib));
for ista = 1:nsta
   Unorm(ie(ista):ie(ista)+1,:) = ...
      J*(squeeze(Z(:,:,ista,ib))\Unorm(ie(ista):ie(ista)+1,:));
end
