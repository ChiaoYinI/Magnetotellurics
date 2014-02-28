ctitle1 = [ 'H_x Polarization' ;...
            'H_y Polarization'];
ctitle2 = [ 'Real Part' ;...
            'Imag Part'];
nplt = sum(plot_ch); 
if nplt > 3
  nwide = 2;
else
  nwide = 1;
end
nplt2 = ceil(nplt/nwide);
kk = 0;
x0 = .06;
y0 = .04;
dx = (rect_ch_plot(3)-(nwide)*x0 )/nwide;
dy = (rect_ch_plot(4)-(nplt2+1)*y0)/nplt2;
l = polarization;
for k = 3:nt
  if plot_ch(k)
     kk = kk + 1;
     ii = ceil(kk/nwide);
     jj = kk - nwide*(ii-1);
     x = x0*jj + dx*(jj-1);
     y = y0*(nplt2-ii+1) + dy*(nplt2-ii);
     rect = [ x y dx dy ];

     if ReIm == 'real'
        hax(kk) = axes('Position',rect); 
        semilogx(10.^t,squeeze(real(TF(k,l,:))),'r*');
        hold on ; 
        semilogx(10.^t,squeeze(real(TF_smth(k,l,:))),'b-');
        set(gca,'fontweight','bold','Xlim',Tlim)
     else
        hax(kk) = axes('Position',rect); 
        semilogx(10.^t,squeeze(imag(TF(k,l,:))),'r*');
        hold on ; 
        semilogx(10.^t,squeeze(imag(TF_smth(k,l,:))),'b-');
        set(gca,'fontweight','bold','Xlim',Tlim)
     end
     xtxt = .8; ytxt = .8;
     h = text(xtxt,ytxt,ctitle(k-2,:),...
	'Units','normal',...
	'FontWeight','bold');
     if( kk == 1)
           htxt = text(.1,.8,...
		[ ctitle1(polarization,:) '  :: ' ctitle2(reim,:)],...
                 'Units','normal');
           set(htxt, 'FontSize',12,'FontWeight','Bold','FontAngle','Italic')
     end
  end
end
