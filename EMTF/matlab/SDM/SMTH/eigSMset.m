flambda_wts = flambda*diag(fwt);
flambda_wts = flambda_wts.*flambda_scales;
flambda_wts(neig_plt+1:nt,:) = ...
   ones(nt-neig_plt,1)*flambda_wts(neig_plt,:);
flambda_wts = 1./flambda_wts;

%  set up orthogonal polynomials
%   limits in log period for interpolation
Tmin = log10(min(Sdms.T));
Tmax = log10(max(Sdms.T));
Trange = Tmax-Tmin;
Tmin = Tmin - .10*(Trange);
Tmax = Tmax + .10*(Trange);
maxdeg = nb/2;
maxdeg = min(maxdeg,15);
N = 200;
K = ceil(maxdeg);
[Q,SR,V] = orthp2nd(Tmin,Tmax,N,K);

%  Set up X (design matrix) and "data"
d = zeros(nt,nb,2)+i*zeros(nt,nb,2);
t = log10(Sdms.T);
X = zeros(nt,nb,nt-2,K+1,2)+i*zeros(nt,nb,nt-2,K+1,2);
%  orthogonal polynomials evaluated at log10(periods)
P = [];
for k=0:K
   P = [P t.^k];
end
P = P*V;
% loop over frequency bands
for ib = 1:nb
   Ut = Sdms.U(:,:,ib);
   Ut = diag(1./sigma(:,ib))*Ut*diag(flambda_wts(:,ib));
   d(:,ib,:) = -Ut(1:2,:)';
   d(:,ib,1) = d(:,ib,1)*sigma(1,ib);
   d(:,ib,2) = d(:,ib,2)*sigma(2,ib);
   for k = 0:K
      for l = 1:2
         X(:,ib,:,k+1,l) = sigma(l,ib)*Ut(3:nt,:)'*P(ib,k+1);
         X(:,ib,:,k+1,l) = squeeze(X(:,ib,:,k+1,l))*diag(scales(3:nt,l,ib));
      end
   end
end
nr = nt*nb;nc = (nt-2)*(K+1);
X = reshape(X,[nr,nc,2]);
d = reshape(d,[nr,2]);

XX = zeros(nc,nc,2)+i*zeros(nc,nc,2);;
Xd = zeros(nc,2)+i*zeros(nc,2);;
XX(:,:,1) = X(:,:,1)'*X(:,:,1);
XX(:,:,2) = X(:,:,2)'*X(:,:,2);
Xd(:,1) = X(:,:,1)'*d(:,1);
Xd(:,2) = X(:,:,2)'*d(:,2);

%  set up roughness penalty matrix
temp = ones((nt-2),1)*SR';
temp = reshape(temp,[(nt-2)*(K+1),1]);
temp = reshape(temp,[(nt-2)*(K+1) 1]);
Ruff = diag(temp);
