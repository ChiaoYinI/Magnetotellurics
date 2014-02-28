function [alpha,res,misfit,datsz] = eigSMfit(nu,nu_comp)
%   Usage [alpha,res,misfit,datsz] = eigSMfit(nu,nu_comp);

global alpha XX Xd Ruff X d

res = d;
alpha = Xd;
[nn,dum] = size(alpha);
nt2 = length(nu_comp);
K = nn/nt2;
R = ones(K,1)*nu_comp(1,:);
R = reshape(R',nn,1);
R = Ruff*R;
R = diag(R);
XXt = squeeze(XX(:,:,1));
alpha(:,1) = (XXt + nu(1)*R)\Xd(:,1);
res(:,1) = d(:,1) - squeeze(X(:,:,1))*alpha(:,1);

R = ones(K,1)*nu_comp(2,:);
R = reshape(R',nn,1);
R = Ruff*R;
R = diag(R);
XXt = squeeze(XX(:,:,2));
alpha(:,2) = (XXt + nu(2)*R)\Xd(:,2);
res(:,2) = d(:,2) - squeeze(X(:,:,2))*alpha(:,2);
 
datsz = sum(sum(conj(d).*d));
misfit = sum(sum(conj(res).*res));

