function [V,Vmag] = setPWG(lats,lons,ih,theta)
%Usage: [V] = setPWG(lats,lons,ih,theta);

%  make centered x-y station coordinates (in km ... x is north)
a = 6.3708e+6;
lons(lons>180) = lons(lons>180)-360;
kmPdg = a*2*pi/360;
clat = mean(lats);
clon = mean(lons);
y = (lons-clon)*cos(clat*pi/180)*kmPdg;
x = (lats-clat)*kmPdg;
x = x/1000;  y = y/1000;

if nargin >= 4
   %  rotate x/y station coordinates so that x is 
   % oriented theta degrees E of N
    c = cos(theta*pi/180);
    s = sin(theta*pi/180);
    xPrime = c*x+s*y;
    yPrime = -s*s+c*y;
    x = xPrime; y = yPrime;
end

%  Last element of ih must be 1 greater than ntot!
nTot = ih(end)-1;
V = zeros(nTot,5);
iX = ih(1:end-1);
iY = iX+1;
V(iX,1) = 1;
V(iY,2) = 1;
V(iX,3) = x';
V(iY,3) = y';
V(iX,4) = x';
V(iY,4) = -y';
V(iX,5) = y';
V(iY,5) = x';
for k = 1:5
    Vmag(k) = sqrt(real(V(:,k)'*V(:,k)));
    V(:,k) = V(:,k)/Vmag(k);
end