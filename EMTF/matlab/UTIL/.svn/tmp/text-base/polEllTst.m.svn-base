function polEllTst(x,y,dxr,dyr,dxi,dyi,clr,ArrSc)

% This function plots polarization ellipse centered
% at (x(n),y(n)) corresponding to the complex vectors
%  ( dxr(n) + i* dxi(n) , dyr(n) + i*dyi(n)); 
%   n = 1:N ;  clr is the color of the line used for 
%  the polarization ellipse  
%  The routine calls ellipse.m to compute the ellipse
%  Scaling of the complex vector into the (x,y) space
%   must be done before calling this routine

ang1 = pi*15/180; 
ang2 = pi*15/180; 
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
  plot(cuv(1,1:end-i2),cuv(2,1:end-i2),clrsym,'LineWidth',1.5)
  if nargin ==8
     x = cuv(1,end-2)-cuv(1,end-1); y = cuv(2,end-2)-cuv(2,end-1);
     s = sqrt(x^2+y^2);
     X = [x,y];
     ii = 1;
     while s < ArrSc
        ii = ii + 1;
        x = cuv(1,end-ii-1)-cuv(1,end-ii); 
        y = cuv(2,end-ii-1)-cuv(2,end-ii);
        s = s+sqrt(x^2+y^2);
        if(s <= .8*ArrSc) NarrC = ii; end 
        X = [X;x,y];
     end
     Narr = ii;
     X1 = X*[ca sa; -sa ca];
     X2 = X*[ca -sa; sa ca];
     X1 = [cuv(1,end-1) cuv(2,end-1); X1];
     X2 = [cuv(1,end-1) cuv(2,end-1); X2];
     X1 = cumsum(X1);
     X2 = cumsum(X2);
     pX = [X1(:,1);cuv(1,end-NarrC);flipud(X2(:,1))];
     pY = [X1(:,2);cuv(2,end-NarrC);flipud(X2(:,2))];
%     plot(pX,pY,'k-','LineWidth',2);
     h = patch(pX,pY,clr);
%     set(h,'EdgeColor','none')
  end
end
set(gca,'Box','on');
