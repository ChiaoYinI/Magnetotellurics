function [Sdms,err] = loadsdms(fid_sdm,irecl,nbt,nt);
% loads nt x nt Spectral density matrices from *.S0
%  file connected to unit fid_sdm 
%  using record length info info in irecl, for all nbt
% frequency bands.  Call sdm_init first
%  USAGE:  [Sdms,err] = loadSdms(fid_sdm,irecl,nbt,nt)

err = 0;
T = zeros(nbt,1);
nf = T;
var = zeros(nt,nbt);
Sig = zeros(nt,nbt);
S = zeros(nt,nt,nbt)+i*zeros(nt,nt,nbt);
lambda = zeros(nt,nbt);
U = zeros(nt,nt,nbt)+i*zeros(nt,nt,nbt);

for ib = 1:nbt
   [T(ib),nf(ib),var(:,ib),S(:,:,ib),err] = ...
      sdm_in(fid_sdm,nt,ib,irecl);
   if( err < 0)
      return
   end
   %  solve generalized eigenvalue problem
   [u1,eval1] = eig(S(:,:,ib),diag(var(:,ib))); 
   eval1 = real(eval1);
   %  scale into eigenvectors of scaled sdm
   u1 = diag(sqrt(var(:,ib)))*u1;
   %  now normalize so that eigenvectors are orthonormal
   normU = sum(conj(u1).*u1,1); normU = 1./sqrt(normU);
   u1 = u1*diag(normU);
   % eigenvectors of scaled sdm ...
   %   u = sqrt(diag(var))*u;
   %  make sure eigenvalues are in correct order ...
   [temp,ind] = sort(diag(eval1));
   ind = ind([nt:-1:1]);
   U(:,:,ib) = u1(:,ind);
   lambda(:,ib) = temp([nt:-1:1]) ;
   Sig(:,ib) = real(diag(squeeze(S(:,:,ib)))); 
end

Sdms = struct('T',T,'nf',nf,'var',var,'S',S,'U',U,...
   'lambda',lambda,'Sig',Sig);
