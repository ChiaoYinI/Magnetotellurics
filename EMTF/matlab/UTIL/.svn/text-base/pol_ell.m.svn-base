function pol_ell(x,y,dxr,dyr,dxi,dyi,clr,ArrSc)

% This function plots polarization ellipse centered
% at (x(n),y(n)) corresponding to the complex vectors
%  ( dxr(n) + i* dxi(n) , dyr(n) + i*dyi(n)); 
%   n = 1:N ;  clr is the color of the line used for 
%  the polarization ellipse  
%  The routine calls ellipse.m to compute the ellipse
%  Scaling of the complex vector into the (x,y) space
%   must be done before calling this routine

ang1 = pi*20/180; 
ang2 = pi*20/180; 
ca = cos(ang1);
sa = sin(ang1);
sa2 = sin(ang2);
hold on
N = length(x);
clrsym = [clr '-'];
for n = 1:N
  u = dxr(n)+i*dxi(n);
  v = dyr(n)+i*dyi(n);
  cuv = ellipse(u,v,x(n),y(n));
  ii = size(cuv);
  i2 = nargin == 8;
  i1 = max(fix(ii(2)/15),1);
  plot(cuv(1,1:end-i2),cuv(2,1:end-i2),clrsym)
  if nargin ==8
     a0x = cuv(1,end-1);
     a0y = cuv(2,end-1);
     aT = cuv(:,end-1)-cuv(:,end-1-i1) ;
     aT = aT/norm(aT);
     ax1 = a0x-ArrSc*(ca*aT(1)-sa*aT(2));
     ay1 = a0y-ArrSc*(sa*aT(1)+ca*aT(2));
     ax2 = a0x-ArrSc*(ca*aT(1)+sa*aT(2));
     ay2 = a0y-ArrSc*(-sa*aT(1)+ca*aT(2));
     ax3 = a0x - ArrSc*aT(1)*(1-sa2);
     ay3 = a0y - ArrSc*aT(2)*(1-sa2);
     pX = [ax1 a0x ax2 ax3]; pY = [ay1 a0y ay2 ay3]; 
     h = patch(pX,pY,clr);
     set(h,'EdgeColor','none')
  end
end
