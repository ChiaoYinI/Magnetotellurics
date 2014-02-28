X = [272.7188,272.7188];
h = get(gcf,'Children');
nh = length(h);
for k = 1:nh
   if(isempty(get(h(k),'Tag')) )
      axes(h(k))
      Y = get(h(k),'Ylim');
      hold on
      plot(X,Y,'--','LineWidth',3,'color',[.0,.0,.0]);
   end
end
