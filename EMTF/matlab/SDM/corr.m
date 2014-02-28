function [c] = mvcorr(S,i1,i2,i3)

%   returns multiple coherence of channels listed in i1 with those in i2
%    (i.e., i2 gives the predicting channels, i1 the predicted)
%   if i3 is not empty coherences are partial multiple coherence
%   given channels in i3

[n,m,nb] = size(S);
n1 = length(i1);
n2 = length(i2);
n3 = length(i3);

c = zeros(n1,nb);
if(n3 == 0 )
  for ib = 1:nb
     C =  S(i1,i2,ib)/S(i2,i2,ib)*S(i2,i1,ib);
     c(:,ib) = diag(C)./diag(S(i1,i1,ib));
  end
else
   i0 = [reshape(i1,[n1,1]) ; reshape(i2,[n2,1])];
   for ib = 1:nb
      C1 = S(i0,i0,ib) - S(i0,i3,ib)/S(i3,i3,ib)*S(i3,i0,ib);
      C = C1(1:n1,n1+1:n1+n2)/C1(n1+1:n1+n2,n1+1:n1+n2)*C1(n1+1:n1+n2,1:n1);
      c(:,ib) = diag(C)./diag(C1(1:n1,1:n1));
   end
end 
     
