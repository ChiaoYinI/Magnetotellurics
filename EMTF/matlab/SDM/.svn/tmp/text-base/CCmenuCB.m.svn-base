action = get(gcbo,'Tag');
hCCmenu = get(gcbo,'Parent');
hSDM = get(hCCmenu,'UserData');
UD = get(hSDM,'UserData');

switch action
   case 'PLOT'
      nt = get(gcbo,'UserData');
      Grp1 = zeros(nt,1);
      Grp2 = zeros(nt,1);
      for k = 1:nt
         Grp1(k) = get(findobj('Parent',hCCmenu,'Tag','Group 1',...
			'UserData',k),'Value');
         Grp2(k) = get(findobj('Parent',hCCmenu,'Tag','Group 2',...
			'UserData',k),'Value');
      end
      ind1 = find(Grp1);         
      ind2 = find(Grp2);         
      hCC = SDMsub(UD.Sdms,ind1,ind2);

   case 'Group 1'
      k = get(gcbo,'UserData');
      l = get(gcbo,'Value');
      if(l == 1)
         set(findobj('Parent',hCCmenu,'Tag','Group 2',...
	   'UserData',k),'Value',0);
      end
   case 'Group 2'
      k = get(gcbo,'UserData');
      l = get(gcbo,'Value');
      if(l == 1)
         set(findobj('Parent',hCCmenu,'Tag','Group 1',...
	   'UserData',k),'Value',0);
      end
   case 'QUIT'
     delete(hCCmenu);
end
