%  plot H and E components of one eigenvector on a single plot

function hfig = evplt_HE(rect,ll_lim,sfac,uH,uE,uZ,H_sta,...
   axlab,ctit,l_ellipse,l_label,l_Hz)

stmonitr
if(axlab(2))
   ym = 'auto';
else
   ym = 'manual';
end
if(axlab(1))
   xm = 'auto';
else
   xm = 'manual';
end
temp = size(uH);
nH = temp(2);
temp = size(uE);
nE = temp(2);
ssH = real(sum(sum(uH(3:4,:).*conj(uH(3:4,:)))));
ss = real(   sum(sum (uH(3:4,:).*conj(uH(3:4,:)) ))+ ...
    sum(sum(uE(3:4,:).*conj(uE(3:4,:)))));
lat_range = ll_lim(2)-ll_lim(1);
lon_range = ll_lim(4) - ll_lim(3);
if( ~ l_ellipse) sfac = .8*sfac ; end;
scale = sfac*sqrt(lat_range*lon_range/ss);
hfig = axes('Position',rect);

% plot complex H vectors
y = real(uH(1,:));x = real(uH(2,:));
dxr = real(uH(4,:));dyr = real(uH(3,:));
dxr = dxr*scale;dyr = dyr*scale;
dxi = imag(uH(4,:));dyi = imag(uH(3,:));
dxi = dxi*scale;dyi = dyi*scale;
if l_ellipse 
  pol_ell(x,y,dxr,dyr,dxi,dyi,'g')
else
  quiver(x,y,dxr,dyr,0,'g');
  hold on
  quiver(x,y,dxi,dyi,0,'g--');
end

%  
%if l_Hz
%   cmp = hsv(64);
%   nH = length(x);
%   scaleH = sqrt(ssH/nH)
%   temp = size(uZ);
%   nZ = temp(2);
%   for k = 1:nZ
%      x1 = x(k);
%      y1 = y(k);
%      Hz_amp = abs(uZ(3,k))
%      Hz_ph = atan2(real(uZ(3,k)),imag(uZ(3,k)));
%      Hz_ph = mod(Hz_ph+2*pi,2*pi);
%      %  size of marker, larger for larger Hz
%      ms = 2 + 60*(Hz_amp/scaleH)
%      %  color gives phase
%      ind = floor(32*Hz_ph/pi)+1
%      mc = cmp(ind,:); 
%      plot(x1,y1,'Marker','o','MarkerSize',ms,...
%	'MarkerFaceColor',mc,'MarkerEdgeColor',mc);
%   end 
%else
%   plot(x,y,'wo');
%end
plot(x,y,'wo');

if(l_Hz)
   nH = length(x);
   scaleH = sqrt(ssH/nH)
   temp = size(uZ);
   nZ = temp(2);
   for k = 1:nZ
      Hz = uZ(3,k)/scaleH;
      HzR = fix(real(Hz*100))/100;
      HzI = fix(imag(Hz*100))/100;
      HzTxtR = num2str(HzR,2);
      HzTxtI = num2str(HzI,2);
      if HzR < 1
         %HzTxtR = HzTxtR(2:end);
      end
      if HzI < 1
         %HzTxtI = HzTxtI(2:end);
      end
      HzTxt = [HzTxtR ',' HzTxtI];
      xx = x(k)+lat_range/15;
      text(xx,y(k),char(HzTxt),'Color',[1,1,1],'FontSize',9);
   end
end

%  label sites
if(l_label & ~l_Hz)
  xx = x+lat_range/10;
  text(xx,y,char(H_sta'),'Color',[1,1,1],'FontSize',11);
end

% plot complex E vectors
y = real(uE(1,:));x = real(uE(2,:));
dxr = real(uE(4,:));dyr = real(uE(3,:));
dxr = dxr*scale;dyr = dyr*scale;
dxi = imag(uE(4,:));dyi = imag(uE(3,:));
dxi = dxi*scale;dyi = dyi*scale;
if(l_ellipse)
  pol_ell(x,y,dxr,dyr,dxi,dyi,'r')
else
  quiver(x,y,dxr,dyr,0,'r');
  hold on
  quiver(x,y,dxi,dyi,0,'r--');
end
fatlines(gca,line_thick)
plot(x,y,'wo')
set(gca,'Ylim',ll_lim(1:2),'Xlim',ll_lim(3:4),'XTickLabelMode',xm, ...
   'YTickLabelMode',ym,'FontWeight','bold','FontSize',11,'Color',[0,0,0]);
title(ctit);
set(get(gca,'title'),'FontWeight','Bold','FontSize',12);
hold off
return
