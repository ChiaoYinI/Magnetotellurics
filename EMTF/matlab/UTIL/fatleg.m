function fatleg(hleg,thick)
%  makes all lines thicker in axis with handle h_ax
hh = get(hleg,'Children');
for k=1:length(hh)
   if(get(hh(k),'Type') == 'line')
      set(hh(k),'LineWidth',thick)
   end
   if(get(hh(k),'Type') == 'text')
      set(hh(k),'FontWeight','bold')
   end
end
