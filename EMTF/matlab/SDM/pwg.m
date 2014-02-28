function [Vpwg]= pwg(S,SdmHd,var,Options)
%  Usage: [Vpwg]= pwg(S,SdmHd,var,Options);

nEvec = Options.nEvec;
theta = Options.theta;

%  make PW and grad vectors
lats = SdmHd.stcor(1,:);
lons = SdmHd.stcor(2,:);
ih = SdmHd.ih;
indH = [ih ih+1];
ihEnd = [ih ih(end)+SdmHd.nch(end)];
nsta = length(lats);
iHonly = [1:2:2*nsta+1];
[nt,dum] = size(S);
[W,Wmag] = setPWG(lats,lons,ihEnd);

%  make eigenvectors (possibly using weights in var (increase
%   z error variances before calling to focus on Horizontal fields) 
   [u1,eval1] = eig(S,diag(var)); 
   eval1 = real(diag(eval1));
   %  scale into eigenvectors of scaled sdm
   %   NOTE:  THESE ARE NOT IN PHYSICAL UNITS!!!!...  but they are
   %   orthonormal after rescaling:
   N = sqrt(var);
   u1 = diag(N)*u1;
   %  make sure eigenvalues are in correct order ...
   [temp,ind] = sort(eval1);
   ind = ind([nt:-1:1]);
   Up = u1(:,ind);
   lambda = temp([nt:-1:1]) ;
   lambda = 1./lambda(1:nEvec);
   Wp = diag(1./N)*W;

if nEvec > 5
   % LS fit of eigenvectors as linear combinations of dealized PWG vectors
   bHat = Wp(indH,:)\Up(indH,1:nEvec);
   res = Up(indH,1:nEvec)-Wp(indH,:)*bHat;
   resCov = res'*res;
   nu = real(Options.evalWt*trace(resCov)/sum(lambda));
   resCov = resCov + nu*diag(lambda);
   [V,D] = eig(resCov);
   d = diag(real(D));
   [d,I] = sort(d);
   V = V(:,I);
%  need to check this!!!
   Upwg = diag(N)*Up(:,1:nEvec)*V(1:nEvec,1:5);
else
   Upwg = diag(N)*Up(:,1:nEvec);
end
if( theta ~= 0 )
    [W,Wmag] = setPWG(lats,lons,ihEnd,theta);
    [Upwg] = Rotate(Upwg,ih,theta);
end
Vpwg = Upwg*inv(W'*Upwg);
for k=1:5
    Vpwg(:,k)=Vpwg(:,k)*Wmag(k);
end
if (theta ~= 0)
    size(Vpwg)
   [Vpwg] = Rotate(Vpwg,ih,-theta);
end

if Options.TeTm
    Vx = (Vpwg(:,3)+Vpwg(:,4))/2;
    Vy = (Vpwg(:,3)-Vpwg(:,4))/2;
    Vpwg(:,3) = Vx;
    Vpwg(:,4) = Vy;
end
