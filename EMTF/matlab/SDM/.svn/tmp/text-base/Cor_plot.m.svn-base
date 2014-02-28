stmonitr;

hfig = get(gcbo,'Parent')
h1 = findobj('Tag','IND','Parent',hfig,'Value',1);
i1 = [];
for k = 1:length(h1)
    i1 = [i1;get(h1(k),'UserData')];
end
h2 = findobj('Tag','DEP','Parent',hfig,'Value',1);
i2 = [];
for k = 1:length(h2)
    i2 = [i2;get(h2(k),'UserData')];
end
h3 = findobj('Tag','GIVEN','Parent',hfig,'Value',1);
i3= [];
for k = 1:length(h3)
    i3 = [i3;get(h3(k),'UserData')];
end
n1 = length(i1); n2 = length(i2); n3 = length(i3);
[c,s] = mvcorr(Sdms.S,i1,i2,i3);
figure('Name','Correlations')
ch_names = [ setstr(csta') setstr(chid')];
ctitle = 'Multiple/Partial Coherence (squared)'
CurveLabels = [];
for k1 = 1:n1
   cl = [ ch_names(i1(k1),:) '/'];
   for k = 1:n2
      cl = [ cl ch_names(i2(k),:) '; ' ]; 
   end
   if n3 ~= 0
      cl = [ cl ' | '];
      for k = 1:n3
         cl = [ cl ch_names(i3(k),:) '; '];
      end
   end
   CurveLabels = [ CurveLabels ; cl ];
end
semilogx(Sdms.T,c');
legendstr = ['legend('];
for k = 1:n1
   legendstr = [ legendstr 'CurveLabels(' num2str(k) ',:),' ];
end
legendstr = [ legendstr '0)' ];
eval(legendstr);
hleg = findobj('Tag','legend');
fatleg(hleg(1),line_thick);
leg_pos = get(hleg(1),'Position');
leg_pos(3:4) = leg_pos(3:4)*leg_scale,l_scale(n1);
set(hleg(1),'Position',leg_pos);

title(ctitle);
fatlines(gca,line_thick);
set(gca,'Ylim',[0,1],'FontSize',12,'FontWeight','bold');
xlabel('Period (sec)');
ylabel('Squared Correlation');
set(get(gca,'Title'),'FontSize',12,'FontWeight','bold')
set(findobj('Tag','ADD'),'enable','on');


