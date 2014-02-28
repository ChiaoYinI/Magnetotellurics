function multiPlot(D,T,F,OPTIONS,linT)
%   high level plotting script for making a set time/frequency
%   of pseudocolor plots on a single page.
%   USAGE:  multiPlot(D,T,F,OPTIONS);
%	D(nT,nF,nPlot) : array to plot; each plot is nT x nF
%	T(nT) : time
%	F(nF) : period
%	OPTIONS : cell array of structures (one for each plot)
% 	 giving axis and title labels, limits, etc
%   NOTE: this version assumes F is logarithmically spaced

if nargin < 5
   linT = 0
end
%  colormap : yellow-centric
cmap = colmap;
x = [0:16]; xi = [0:.125:16];
cmapI = zeros(3,length(xi));
for k = 1:3
   cmapI(k,:) = interp1(x,cmap(:,k),xi);
end
cmapI = cmapI';
cmapI = flipud(cmapI);

% array sizes
[nT,nF,nPlot] = size(D);

%  axis postions depend on number of plots ... this makes sense
%   for no more than 4

%  f is fraction of vertical space used for actual plot
f = .8;
yTot =  .85;

yHt = (yTot/nPlot)*f;
ySpace = (yHt/nPlot)*(1-f);
for k = 1:nPlot
   y0 = (1-yTot)+(k-1)*(yHt+ySpace);
   RectAx(nPlot-k+1,:) = [.1,y0,.9,yHt];
end

%  log of period (OR FREQUENCY)
if linT
   lF = F;
else
   lF = log10(F);
end

figure('Position',[100,100,900,600],...
        'PaperPosition',[1,1,10,6.5],...
        'PaperOrientation','Landscape');

yal = fix(nPlot/2);
if (nPlot == 2 | nPlot == 3)
   yal = 2;
end
for k = 1:nPlot
   data = squeeze(D(:,:,k))';
   axes('Position',RectAx(k,:));
   pcolor(T,lF,data); shading flat;
   if(~isempty(OPTIONS{k}.Caxis))
      caxis(OPTIONS{k}.Caxis); 
   end
   colorbar
   colormap(cmapI);
   set(gca,'FontWeight','demi','FontSize',12)
   if(k<nPlot)
      set(gca,'Xtick',[])
      if(k==1)
         title(OPTIONS{k}.Title);
      end
   else
      xlabel(OPTIONS{k}.TimeAxisLabel,'FontSize',14)
   end
   text(.05,.2,OPTIONS{k}.SubTitle,'Units','Normalized',...
		'FontSize',14,'FontWeight','demi')
   if k == yal
      ylabel(OPTIONS{k}.PeriodAxisLabel,'FontSize',14)
   end
   if OPTIONS{k}.DayLine
      Y = [min(lF),max(lF)]; 
      hold on
      for day = min(fix(T)):max(fix(T))
         X = [day day];
         plot(X,Y,'k-')
      end
      hold off
   end
end

hCB = findobj('Tag','Colorbar');
for k = 1:nPlot
   rectCB = get(hCB(k),'Position')
   rectCB(1) = 0;%rectCB(1) = rectCB(1)-rectCB(3)*.5;
   rectCB(3) = rectCB(3)*.3;
   set(hCB(k),'Position',rectCB,'FontWeight','demi')
   if k == nPlot
      text(-1.2,1.2,OPTIONS{k}.ColorAxisLabel,...
		'Parent',hCB(k),'Units','Normalized',...
		'FontSize',12,'FontWeight','demi')
   end
end
