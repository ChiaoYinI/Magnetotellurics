function [var,sig,IER] = SDMvar(Sdms,grouping)
% Usage: [var,sig,IER] = SDMvar(Sdms,grouping);

IER = 0;
ih = Sdms.Hd.ih;
chid = Sdms.Hd.chid;
nsta = length(ih);
[dum,nt] = size(chid);
[dum,dum,nb] = size(Sdms.S);
%   find indices where components in list change from E to H or H to E
CUST = 0;

switch grouping
  case 'standard'
    Hnum = 72; Enum = 69;
    temp = fix(chid(1,:));
    ind = (temp == Hnum) - (temp == Enum);
    ind = ind(1:nt-1).*ind(2:nt);
    ind = 1+find(ind == -1);
    for ista = 1:nsta
      if(sum(ind == ih(ista) ) == 0)
         ind = [ind ih(ista)];
      end
    end
    ind = sort(ind);
  case 'all'
    ind = [1:nt];
  case 'sites'
    ind = ih;
  case 'custom'
    CUST = 1;
  otherwise
       fprintf(1,'%s \n','Case Not Coded in evecCbk')
end
  
if CUST
   ngrp = length(ih);
else
   ngrp = length(ind);
   ind = [ind nt+1];
end

var = zeros(nt,nb);
sig = zeros(nt,nb);
dmu = .9;
o = ones(nt,1);
pmax = 1.0001; pmin = .3;
for ib = 1:nb
   TF = zeros(nt,nt);
   v = zeros(nt,1);
   for igrp = 1:ngrp
      if CUST
         II = ih{igrp};
         ind = ones(nt,1);
         ind(II) = 0;
         JJ = find(ind);
      else
         i1 = ind(igrp);
         i2 = ind(igrp+1)-1;
         II = [i1:i2];
         JJ = [1:i1-1 i2+1:nt];
      end

      S11 = Sdms.S(II,II,ib);
      S22 = Sdms.S(JJ,JJ,ib);
      S12 = Sdms.S(II,JJ,ib);

      scl = diag(1./sqrt(real(diag(S22))));
      [V,d] = eig(scl*S22*scl);
      d = diag(d);
      dInv = zeros(size(d));
      dInv(d>.01) = 1./(d(d>.01));
      S22inv = scl*V*diag(dInv)*V'*scl;
      TF(II,JJ) = S12*S22inv;
      S1g2 = S11 - TF(II,JJ)*S12';
      v(II) = real(diag(S1g2));
      sig(II,ib) = real(diag(S11));
   end
   if any(v == 0)
      IER = -1;
      return
   end 
   mu = 1;
   done = 0;
   while 1-done
      X = mu*abs(TF).*abs(TF);
      X = diag(1./v)*X*diag(v) + eye(nt);
      vhat = X\o; 
      if(max(vhat) > pmax | min(vhat) < pmin )
         mu = mu*dmu ;
      else
         done = 1;
      end
      var(:,ib) = vhat.*v;
    end
end
