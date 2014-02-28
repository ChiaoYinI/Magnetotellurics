function SetCorCh

%  this just ensures that only one box in a row is checked
kk = get(gcbo,'Value');
k = get(gcbo,'UserData');
hfig = get(gcbo,'Parent');

if kk==1
   set(findobj('Parent',hfig,'style','checkbox','UserData',k),'Value',0);
   set(gcbo,'Value',1);
end
