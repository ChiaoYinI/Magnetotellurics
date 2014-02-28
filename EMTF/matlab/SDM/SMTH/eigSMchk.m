indN = [1:2];
[TF] = eigTF(Sdms.U,Sdms.var,indN);
TF_smth = TF;
nalpha = (K+1)*(nt-2);
for k =3:nt
  TF_smth(k,1,:) = P * alpha(k-2:nt-2:nalpha,1); 
  TF_smth(k,2,:) = P * alpha(k-2:nt-2:nalpha,2); 
end
TF_smth = TF_smth.*scales;

TF_dot = zeros(nt,2,nb);
Smth_dot = zeros(nt,2,nb);
TF_tot = zeros(2,nb);
Smth_tot = zeros(2,nb);
for ib = 1:nb
  SC = diag(1./sqrt(Sdms.var(:,ib)));
  TFsc = SC*squeeze(TF_smth(:,1,ib));
  Smth_tot(1,ib) = TFsc'*TFsc;
  Smth_dot(:,1,ib) = (abs((squeeze(Sdms.U(:,:,ib))'*TFsc)).^2)/Smth_tot(1,ib); 
  TFsc = SC*squeeze(TF_smth(:,2,ib));
  Smth_tot(2,ib) = TFsc'*TFsc;
  Smth_dot(:,2,ib) = (abs((squeeze(Sdms.U(:,:,ib))'*TFsc)).^2)/Smth_tot(2,ib); 
  TFsc = SC*squeeze(TF(:,1,ib));
  TF_tot(1,ib) = TFsc'*TFsc;
  TF_dot(:,1,ib) = (abs((squeeze(Sdms.U(:,:,ib))'*TFsc)).^2)/TF_tot(1,ib); 
  TFsc = SC*squeeze(TF(:,2,ib));
  TF_tot(2,ib) = TFsc'*TFsc;
  TF_dot(:,2,ib) = (abs((squeeze(Sdms.U(:,:,ib))'*TFsc)).^2)/TF_tot(2,ib); 
end
temp = squeeze(Smth_dot(1:2,1,:));
temp = cumsum(temp);
cum_min1 = min(temp(2,:));
temp = squeeze(Smth_dot(1:2,2,:));
temp = cumsum(temp);
cum_min2 = min(temp(2,:));
cum_min1 = min(cum_min1,cum_min2);
cum_min1 = floor(50*cum_min1)/50;
fit_lim = [cum_min1,1];
