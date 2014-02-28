function [c,s,t] = mvcorr(S,i1,i2,i3)
%  Usage:  [c,s,t] = mvcorr(S,i1,i2,i3);
%   computes multiple coherence between selected channels 
%   Inputs: S = SDM (complex cross-product matrix)
%           i1 = predicted channels
%           i2 = predicting channels
%           i3 is optional; if not empty coherences are 
%              partial multiple coherences, given channels in i3
%   Outputs: c = squared coherences
%	     s = variance of predicted channels
%            t = transfer function between predicting/predicted channels
%                (conditional tf, if there are any "given" channels)

[n,m,nb] = size(S);
n1 = length(i1);
n2 = length(i2);
n3 = length(i3);

c = zeros(n1,nb);
s = zeros(n1,nb);
t = zeros(n1,n2,nb);
if(n3 == 0 )
  for ib = 1:nb
     t(:,:,ib) = S(i1,i2,ib)/S(i2,i2,ib);
     C =  t(:,:,ib)*S(i2,i1,ib);
     c(:,ib) = diag(C)./diag(S(i1,i1,ib));
     s(:,ib) = diag(S(i1,i1,ib));
  end
else
   i0 = [reshape(i1,[n1,1]) ; reshape(i2,[n2,1])];
   for ib = 1:nb
      C1 = S(i0,i0,ib) - S(i0,i3,ib)/S(i3,i3,ib)*S(i3,i0,ib);
      t(:,:,ib) = C1(1:n1,n1+1:n1+n2)/C1(n1+1:n1+n2,n1+1:n1+n2);
      C = t(:,:,ib)*C1(n1+1:n1+n2,1:n1);
      c(:,ib) = diag(C)./diag(C1(1:n1,1:n1));
      s(:,ib) = diag(C1(1:n1,1:n1));
   end
end 
c = real(c);
s = real(s);
