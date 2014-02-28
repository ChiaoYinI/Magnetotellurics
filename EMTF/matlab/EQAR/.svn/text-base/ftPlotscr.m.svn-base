   function [h,tMed,fMed] = ftPlot(S,FT,plotWhat,plotOpt)
   %this script is ftPlot except that 
   
   
%   Usage [h,tMed,fMed] = ftPlot(S,plotWhat,plotOpt)
%   generic frequency/time plotting code
%   takes two cell arrays, and two structures as input: 
%  Inputs:
%   1)  S : cell array of data structures, one cell for each time segment
%   2)  FT : structure containing time info for labelin x axis,
%           frequency info for labeling y axis
%   3)  plotWhat a cell array containing a keyword defining what
%         should be plotted (e.g., rho, phi, etc.) and any parameters
%         needed to compute this from data in S (e.g., indicies of impedance)
%   4)  plotOpt structure giving  plotting options,k
%                     such as specified caxis limits, titles
%                     control parameters for scaling, 
%                     subtracting time mean, etc.
%  Outputs:
%   1)  figure handle; h = -1 on (some) errors
%   2)  tMed = median (over time) of plotted variable as a function of frequency
%   3)  fMed = median (over frequency) of plotted variable as a function of time

cmapI = yellowZero;

nSeg = length(S);
T = FT.T;
[i1,i2] = size(T);
if(i1 == 1) T = T'; end
nFreq = length(T);
for k = 1:length(S)
   if(1-S{k}.miss)
      firstGood = k;
      break;
   end
end
if nFreq ~= length(S{firstGood}.T)
   h = -1
   fprintf(2,'%s','Error: data structure, FT not consistent')
   return
end
f = zeros(nFreq,nSeg)/0;
%  define array to plot (f) using data in "S", instructions in
%        "plotWhat"
switch plotWhat{1}
   case 'rho' 
   %   Apparent resistivity
   %   In this case structure S should be of type Z2x2
      if (S{firstGood}.type ~= 'Z2x2')
         fprintf(2,'%s',...
           'Error: plot type, data structure inconsistent')
         h = -1
         return
       end
       %   for rho, second cell of plotWhat contains indices 
       %    of impedance array : i, j and optional impedance number
       indZ = plotWhat{2};
       if length(indZ) == 2
         indZ(3) = 1;
       end 
       ErrCutOff = plotWhat{3};
       for k = 1:nSeg
          if(S{k}.miss == 0)
             Z = squeeze(S{k}.Z(indZ(1),indZ(2),:,indZ(3)));
             rho =  (abs(Z).^2)/5;
             if ErrCutOff > 0
                Noise = squeeze(S{k}.E(indZ(1),indZ(1),:,indZ(3)));
                Sig = squeeze(S{k}.S(indZ(2),indZ(2),:,indZ(3)));
                fracErr =  2*sqrt(Noise.*Sig)./abs(Z);
                rho(fracErr > ErrCutOff) = NaN;
             end
             [i1,i2] = size(rho);
             if(i1 == 1) rho = rho'; end
             f(:,k) = log10(T.*rho);
           end
       end

   case 'phi' 
   %   Impedance Phase
   %   In this case structure S should be of type Z2x2
      if (S{firstGood}.type ~= 'Z2x2')
         fprintf(2,'%s',...
           'Error: plot type, data structure inconsistent')
         h = -1
         return
       end
       %   for phi, second cell of plotWhat contains indices 
       %    of impedance array : i, j and optional impedance number
       indZ = plotWhat{2};
       if length(indZ) == 2
         indZ(3) = 1;
       end 
       ErrCutOff = plotWhat{3};
       for k = 1:nSeg
          if(S{k}.miss == 0)
             Z = squeeze(S{k}.Z(indZ(1),indZ(2),:,indZ(3)));
             phi =  atan(imag(Z)./real(Z));
             if ErrCutOff > 0
                Noise = squeeze(S{k}.E(indZ(1),indZ(1),:,indZ(3)));
                Sig = squeeze(S{k}.S(indZ(2),indZ(2),:,indZ(3)));
                fracErr =  2*sqrt(Noise.*Sig)./abs(Z);
                phi(fracErr > ErrCutOff) = NaN;
             end
             [i1,i2] = size(phi);
             if(i1 == 1) phi = phi'; end
             f(:,k) = phi*180/pi;
           end
       end
   case 'tfp' 
   %   Phase of system transfer function
   %   In this case structure S should be of type sysTF
      if (S{firstGood}.type(1:4) ~= 'sysT')
         fprintf(2,'%s',...
           'Error: plot type, data structure inconsistent')
         h = -1
         return
       end
       %   for sysTF, second cell of plotWhat contains indices 
       %    of system TF array 
       indTF = plotWhat{2};
       for k = 1:nSeg
          if(S{k}.miss == 0)
             TF = S{k}.TF(indTF,:);
             phi =  atan(imag(squeeze(TF))./real(squeeze(TF)));
             [i1,i2] = size(phi);
             if(i1 == 1) phi = phi'; end
             f(:,k) = phi*180/pi;
           end
       end
   case 'tfa' 
   %   Amplitude of system transfer function
   %   In this case structure S should be of type sysTF
      if (S{firstGood}.type(1:4) ~= 'sysT')
         fprintf(2,'%s',...
           'Error: plot type, data structure inconsistent')
         h = -1
         return
       end
       %   for sysTF, second cell of plotWhat contains indices 
       %    of system TF array 
       indTF = plotWhat{2};
       for k = 1:nSeg
          if(S{k}.miss == 0)
             TF = S{k}.TF(indTF,:);
             amp =  log10(abs(TF));
             [i1,i2] = size(amp);
             if(i1 == 1) amp = amp'; end
             f(:,k) = amp;
           end
       end
        
    otherwise
       fprintf(2,'%s','Error: This option not coded')
       h = -1
       return
    end  %  switch

    %  if an array of "scaling or phase factors" 
    %  (shifts on a log scale for each time segment)
    %  is provided use to adjust plotted vbl before plotting
    %  Idea is that this might be used to correct calibration
    if exist(plotOpt.fMean)
       if length(plotOpt.fMean) == nSeg
          for k = 1:nSeg
             f(:,k) = f(:,k)-plotOpt.fMean(k);
          end
       end
    end
    %  if an array of "scaling or phase factors" 
    %  (shifts on a log scale for each frequency band)
    %  is provided use to adjust plotted vbl before plotting
    %  Idea is that this might be used to subtract an input long
    %  term mean instead of subtracting the time median of the 
    %  plotted vbl
    if exist(plotOpt.tMean)
       if length(plotOpt.tMean) == nFreq
          for k = 1:nFreq
             f(k,:) = f(k,:)-plotOpt.tMean(k);
          end
       end
    end

   %  compute median over time for each frequency band,
   %   and over frequency for each time segment of the array to
   %   be plotted
   for k = 1:nFreq
       temp = f(k,:);
       temp = temp(find(~isnan(temp)));
       if length(temp > 0)
          tMed(k) = median(temp);
       else
          tMed(k) = NaN;
       end
   end 
   fMed = zeros(nSeg,1);
   for k = 1:nSeg
       temp = f(:,k);
       temp = temp(find(~isnan(temp)));
       if length(temp > 0)
          fMed(k) = median(temp);
       else
          fMed(k) = NaN;
       end
   end 
     
   if plotOpt.subMed == 1
      %  subtract time median and plot residuals
      for k = 1:nFreq
         f(k,:) = f(k,:)-tMed(k);
      end
   end
%!pwd
f2001 = f;
save f2 f2001
 %  f
%  Plot array f
   h = figure('Position',plotOpt.Pos, ...
       'PaperPosition',plotOpt.PaperPos,'PaperOrientation','Landscape')
   Y = T;
   ymin = ceil(min(log10(T)));
   ymax = floor(max(log10(T)));
   yt = 10.^[ymin:ymax];
   ystep = Y(2)-Y(1);
   X = [FT.t1:FT.step:(nSeg-1)*FT.step+FT.t1];
   xstep = X(2)-X(1);
   xtStep = fix(50/xstep);
   xt = round(X(2:xtStep:end)/5)*5
   for k = 1:length(xt)
     xl{k} = num2str(mod(xt(k),365));
   end
   f = [f f(:,end)]; X = [X X(end)+xstep];
   f = [f; f(end,:)];Y = [Y;Y(end)+ystep];
   pcolor(X,Y(2:end),f(2:end,:));shading flat;colormap(cmapI);
   set(gca,'Position',[.1,.15,.85,.75]) ;
   set(gca,'Ytick',yt,'YScale','log')
   %set(gca,'Xtick',xt,'XtickLabel',xl)
   axis('xy')
   set(gca,'FontWeight','demi','FontSize',12)
   ylabel('Period (s)')
   xlabel('Day')
   lims = plotOpt.lims;
   if(length(lims) == 0)
      lims = defaultLims(f)
   end
   title(plotOpt.title)
   caxis(lims)
   %colorbar
   %hCB = findobj('Tag','Colorbar','Parent',gcf);
   %rectCB = get(hCB,'Position');
   %rectCB(1) = rectCB(1)-rectCB(3)*.5;
   %rectCB(3) = rectCB(3)*.3;
   %set(hCB,'Position',rectCB,'FontWeight','demi','FontSize',12);
   %text('string',plotOpt.ColorAxisLabel,'Position',plotOpt.ColorAxisLabelPos,...
   %             'Parent',hCB,'Units','Normalized',...
   %             'FontSize',12,'FontWeight','demi',...
   %		'Rotation',90,'HorizontalAlignment','center')

