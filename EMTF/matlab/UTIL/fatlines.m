function fatlines(h_ax,thick)
%  makes all lines thicker in axis with handle h_ax
hh = get(h_ax,'Children');
for k=1:length(hh)
   obtype = get(hh(k),'Type');
   if(obtype(1:4) == 'line')
      set(hh(k),'LineWidth',thick)
   end
end
