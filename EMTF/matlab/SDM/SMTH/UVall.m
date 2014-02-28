i = sqrt(-1);
Usm = zeros(nt,2,nb)+i*zeros(nt,2,nb);
Vsm = zeros(nt,maxCN,nb)+i*zeros(nt,maxCN,nb);
Uperp = Vsm;

U_snr = Usm;
V_snr = Vsm;
Up_snr = Vsm;
neigmax = maxCN;
Sig_MT = zeros(2,2,nb);
Sig_CN = zeros(neigmax,neigmax,nb);
Neig = zeros(nb,1);
for ib = 1:nb
   UVset;
   [V1,U1,U1perp,V1_snr,U1_snr,U1p_snr,c,SIG_CN,SIG_MT,S_incoh] = ...
                      UVest(S,var,uhat,minCNsig,neigmax);
   neig = length(SIG_CN);
   Neig(ib) = neig;
   Usm(:,:,ib) = U1;
   Vsm(:,1:neig,ib) = V1;
   Uperp(:,1:neig,ib) = U1perp;
   U_snr(:,:,ib) = U1_snr;
   V_snr(:,1:neig,ib) = V1_snr;
   Up_snr(:,1:neig,ib) = U1p_snr;
   Sig_MT(1,1,ib) = SIG_MT(1);
   Sig_MT(2,2,ib) = SIG_MT(2);
   for k = 1:neig
     Sig_CN(k,k,ib) = SIG_CN(k);
   end
end
%   transfer functions relative to first station
TFsm = Usm;
for ib = 1:nb
   U1 = TFsm(1:2,:,ib);
   TFsm(:,:,ib) = TFsm(:,:,ib)/U1;
end
%  normalize U ...
Unorm = sqrt(sum(conj(Usm(1:2,:,:)).*Usm(1:2,:,:),1));
for k = 1:nt
   Usm(k,:,:) = squeeze(Usm(k,:,:)./Unorm);
end
for ib = 1:nb
   for l = 1:2
      Usm(:,l,ib) = chngph(squeeze(Usm(:,l,ib)),[1:2]);
   end
end
