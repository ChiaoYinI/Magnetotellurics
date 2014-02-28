%  chuseset sets up the polarization and predicted channels to plot
%   by settin poluse and chuse

lk = get(gco,'UserData');
l = lk(1); k = lk(2);
if( k == 1 )
   poluse = l;
   set(chhand(1,l,k),'Value',1);
   set(chhand(1,3-l,k),'Value',0);
else
   for m = 2:nlines
      for n = 1:2
         set (chhand(1,n,m),'Value',0);
         chuse (n,m) = 0;
      end	
   end
   set (chhand(1,l,k),'Value',1);
   chuse(l,k) = get(chhand(1,l,k),'Value');
end