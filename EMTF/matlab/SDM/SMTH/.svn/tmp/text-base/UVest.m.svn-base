function [V,U,Uperp,V_snr,U_snr,Up_snr,c,SIG_CN,SIG_MT,S_incoh] = ...
                      UVest(S,var,uhat,eigmin,neigmax)
% Usage [V,U,Uperp,V_snr,U_snr,Up_snr,c,SIG_CN,SIG_MT,S_incoh] = ...
%                      UVest(S,var,uhat,eigmin,neigmax);

[nt,dum] = size(uhat);
if(dum ~= 2)
  fprintf(1,'%s \n','Error : uhat must have two columns')
  return
end
Sigma_N_inv = diag(1./sqrt(var));
Sigma_N = diag(sqrt(var));
Sp = (Sigma_N_inv*S)*Sigma_N_inv;
U = orth(Sigma_N_inv*uhat);
temp = Sp - U*U'*Sp;
temp = temp - temp*U*U';
[V,D] = eig(temp);
D = diag(D);
[D,ind] = sort(D);
neig = sum(D >= eigmin);
neig = min(neig,neigmax)+2;
V = V(:,ind);
Uperp = V(:,nt:-1:nt-neig+3);
W = [Uperp  U ];
SIG = W'*(Sp-eye(nt))*W;
%SIG = W'*Sp*W;
K = neig-2;
SIG11 =SIG(1:K,1:K); 
SIG12 =SIG(1:K,K+1:neig); 
SIG22 =SIG(K+1:neig,K+1:neig); 
c = SIG11\SIG12;
V = Uperp + U*c';
V_snr = orth(V);
A = V_snr'*V;
U_snr = U;
SIG_CN = A*SIG11*A';
[w,d] = eig(SIG_CN);
[temp,ind] = sort(diag(d));
nind = length(ind);
if(nind > 0 ) 
  ind = ind([nind:-1:1]);
  w = w(:,ind);
  V_snr = V_snr*w;
  SIG_CN = temp([nind:-1:1]);
end
SIG_MT = SIG22 - (SIG12'/SIG11)*SIG12;
[w,d] = eig(SIG_MT);
[temp,ind] = sort(diag(d));
ind = ind([2,1]);
w = w(:,ind);
U_snr = U_snr*w;
SIG_MT = temp([2:-1:1]);
S_incoh = Sp - (U_snr*diag(SIG_MT)*U_snr' + V_snr*diag(SIG_CN)*V_snr');
V = Sigma_N*V_snr;
U = Sigma_N*U_snr;
Up_snr = Uperp;
Uperp = Sigma_N*Uperp;
