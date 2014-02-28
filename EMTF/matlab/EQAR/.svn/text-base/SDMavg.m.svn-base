function [SdmAvg] = SDMavg(SDMS,SDMHD,good,grouping);
%  Usage:  [SdmAvg] = SDMavg(SDMS,SDMHD,good,grouping);
%  given input cell arrays of SDM (and header) data
%   structures, and an optional list of "good" SDMs to use,
%   compute an average SDM structure
%
%  grouping is also optional, defaults to 'standard'
%
%  There is no check for consistency/correctness of SDMS in
%   cell array ... need to check before calling

if nargin < 4
  grouping = 'standard';
end

if nargin < 3
   good = ones(size(SDMS));
end

%  use header info from first "good" SDM 
k1 = min(find(good));
[NT,NBT] = size(SDMS{k1}.var);
S = zeros(size(SDMS{k1}.S));
T = SDMS{k1}.T;

%  average over good SDMS  (downweight unusually large amplitude days)
%  firs tcompute weights to limit effect days with exceptioinally strong signal
[wts] = SDMwts(SDMS,SDMHD,good);

%  weighted average
nf = 0;
nFiles = length(SDMS);
wSum = 0;
for k = 1:nFiles
   if(good(k))
      nf = nf+SDMS{k}.nf;
      w = min(wts(k,:));
      wSum = wSum+w;
%      for ib = 1:NBT
%         S(:,:,ib) = S(:,:,ib)+SDMS{k}.S(:,:,ib)*wts(k,ib);
%      end
      S = S+SDMS{k}.S*w;
   end
end
S = S/wSum;
%for ib = 1:NBT
%   S(:,:,ib) = S(:,:,ib)/sum(wts(:,ib));
%end

%   put results into output SDM data structure (with header
%   info in .Hd field)
SdmAvg = struct('T',T,'nf',nf,'S',S,'Hd',SDMHD{k1});

%  compute incoherent noise variances
[var,sig,IER] = SDMvar(SdmAvg,grouping);

%  compute eigenvalues and eigenvectors
lambda = zeros(NT,NBT);
U = zeros(NT,NT,NBT)+i*zeros(NT,NT,NBT);
for ib = 1:NBT
   %  solve generalized eigenvalue problem
   [u1,eval1] = eig(S(:,:,ib),diag(var(:,ib)));
   %  scale into eigenvectors of scaled sdm
   u1 = diag(sqrt(var(:,ib)))*u1;
   %  now normalize so that eigenvectors are orthonormal
   normU = sum(conj(u1).*u1,1); normU = 1./sqrt(normU);
   u1 = u1*diag(normU);
   %  make sure eigenvalues are in correct order ...
   [temp,ind] = sort(diag(eval1));
   ind = ind([NT:-1:1]);
   u1 = u1(:,ind);
   U(:,:,ib) = u1;
   lambda(:,ib) = real(temp([NT:-1:1]));
end

%   put results into output SDM data structure (with header
%   info in .Hd field)
SdmAvg = struct('T',T,'nf',nf,'var',var,'S',S,'U',U,...
       'lambda',lambda,'Sig',sig,'Hd',SDMHD{k1});
