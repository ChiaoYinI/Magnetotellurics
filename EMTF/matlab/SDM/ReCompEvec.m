function [Sdms] = ReCompEvec(Sdms)
%  Usage : [Sdms] = ReCompEvec(Sdms);
% Given Sdms structure, recompute generalized
%   eigenvectors and eigenvalues

[nt,nbt]  = size(Sdms.var);
lambda = zeros(nt,nbt);
U = zeros(nt,nt,nbt)+i*zeros(nt,nt,nbt);

for ib = 1:nbt
   %  solve generalized eigenvalue problem
   [u1,eval1] = eig(Sdms.S(:,:,ib),diag(Sdms.var(:,ib))); 
   eval1 = real(eval1);
   %  scale into eigenvectors of scaled sdm
   u1 = diag(sqrt(Sdms.var(:,ib)))*u1;
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
end

Sdms.U = U;
Sdms.lambda = lambda;
