function [TF] = eigTF(U,var,indN)
%  given evecs of noramlized SDM + incoherent noise variances
%  returns TFs of all components relative to two componentes given by indN
%     indices of assumed normal field channels

[nt,nb] = size(var);
TF =  zeros(nt,2,nb);
for ib = 1:nb
  T = diag(sqrt(var(:,ib)))*U(:,indN,ib);
  TF(:,:,ib) = T/T(indN,:);
end
