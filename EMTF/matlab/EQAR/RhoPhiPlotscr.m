%   Edit plotWhat to define what time window to make plots for, etc.
plotWhat;
%   RhoPhiFileList constructs the list of files using parameters
%   set in plotWhat
RhoPhiFileList
%   Edit CaxisLims to change color axis on Time/frequency plots
CaxisLims

%  set up to save all impedances (and other possible vbles)
%   in a single matlab save file
%  set SAVE = 1 to save as mat file
SAVE = 0;
SaveVbls = 'periods Zall Z2x2all LocSite dayOne dayEnd TimeIncDay';
SaveMatFile = ['Z_' root '_' num2str(dayOne) '-' num2str(dayEnd) '.mat'];

%   plotting options ...  
PlotOpt = struct('Pos',[100,100,900,350],...
	'PaperPos',[1,1,9,3.5],...
        'FontSize',14,...
        'lims',[],...
        'title','',...
        'subMed',subMed, ...
        'tMean',[], ...
        'fMean',[]);

%   default FT options ... FT is structure which defines time/frequency
%   axes for the sections to be plotted
FT = struct('t1',dayOne,'t2',dayEnd','step',TimeIncDay);

%   following define "plotWhat" cell arrays used to
%   specify what ftPlot will plot, and which components of 
%   the impedance files correspond to xy/yx, etc.
%   This is very obtuse!
indsXY = [1,2];
indsYX = [2,1];
indsXY2 = [1,2,2];
indsYX2 = [2,1,2];
rhoXY = cell(3,1);
rhoXY{1} = 'rho';
rhoXY{3} = FracErrCutOff;
rhoYX = rhoXY;
rhoYX2 = rhoYX;
rhoXY2 = rhoYX;
rhoXY{2} = indsXY;
rhoYX{2} = indsYX;
rhoYX2{2} = indsYX2;
rhoXY2{2} = indsXY2;
phiXY = cell(3,1);
phiXY{1} = 'phi';
phiXY{3} = FracErrCutOff;
phiYX = phiXY;
phiYX2 = phiYX;
phiXY2 = phiYX;
phiXY{2} = indsXY;
phiYX{2} = indsYX;
phiYX2{2} = indsYX2;
phiXY2{2} = indsXY2;

%  read impedances for all fileNames in list
%   (with missing data flag set for any missing/incomplete files) 
%   TO decide if impedance file is incomplete:  "normal" number
%   of bands is computes (median), and anything below this is incomplete.
%   THIS ONLY WORKS if all data were processed with the same bset file
%   and if only a small fraction are incomplete ... as is normally the case
[Zall,Z2x2all] = ZallRead(fileNamesMT);display 'read Z success'

%  set periods (in sec) in structure FT
%   take periods from the first good impedance structure
for k = 1:length(Zall)
   if Zall{k}.good == 1
      FT.T = Zall{k}.T;
      break
   end
end

if SAVE
   eval([ 'save ' SaveMatFile ' ' SaveVbls])
end
    
%  Now make plots ... this part could be edited to make different
%    plots; add new parameters to PlotOpt, and new code to ftPlot
%    to implement new features in plotting.
   %   plot rho
   PlotOpt.lims = rhoLims{iSite}(1,:);
   PlotOpt.title = ['rhoXY   :  ',LocSite,' : ' ...
           num2str(YearOne) ' : ',...
           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
   if subMed
      sM = 'dev_'
      PlotOpt.ColorAxisLabel = 'Fractional Deviation'
      %PlotOpt.ColorAxisLabelPos = [5,.5];
      PlotOpt.ColorAxisLabelPos = [5,0];
      PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
   else
      sM = ''
      PlotOpt.ColorAxisLabel = 'log_{1} ohm-m'
      PlotOpt.ColorAxisLabelPos = [3.5,.5];
   end
 
    plotFile = [PlotDir LocSite 'rhoXY_' sM ...
                 num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
    [h,tMed_rhoXY,fMed_rhoXY] = ftPlotscr(Z2x2all,FT,rhoXY,PlotOpt);
   if EQLINE
      plotEQlog
   end
   if SAVEplots
      plot2File(plottype, plotFile);    
   end

   %PlotOpt.lims = rhoLims{iSite}(2,:);
   %PlotOpt.title = ['rhoYX   :  ',LocSite,' : ' ...
   %        num2str(YearOne) ' : ',...
   %        num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
   %if subMed
   %   PlotOpt.ColorAxisLabel = 'Fractional Deviation'
   %   PlotOpt.ColorAxisLabelPos = [5,.5];
   %   PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
   %else
   %   PlotOpt.ColorAxisLabel = 'log_{1} ohm-m'
   %   PlotOpt.ColorAxisLabelPos = [3.5,.5];
   %end
  % 
  %  plotFile = [PlotDir LocSite 'rhoYX_' sM ...
  %               num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
  %  [h,tMed_rhoYX,fMed_rhoYX] = ftPlot(Z2x2all,FT,rhoYX,PlotOpt);
  % if EQLINE
  %    plotEQlog
  % end
  % if SAVEplots
   %    plot2File(plottype, plotFile);  
   %end

   %   plot phi
   %PlotOpt.lims = phiLims{iSite}(1,:)
   %PlotOpt.title = ['phiXY   :  ',LocSite,' : ' ...
   %        num2str(YearOne) ' : ',...
   %        num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
   %plotFile = [PlotDir LocSite 'phiXY_' sM ...
   %             num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
   %             
   %             PlotOpt.ColorAxisLabel = 'Degrees'
   %PlotOpt.ColorAxisLabelPos = [3,.5];
   %if(subMed)
   %   PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
   %end
   %[h,tMed_phiXY,fMed_phiXY] = ftPlot(Z2x2all,FT,phiXY,PlotOpt);
   %if EQLINE
   %   plotEQlog
   %end
   %if SAVEplots
   %    plot2File(plottype, plotFile);  
   %end

%   PlotOpt.lims = phiLims{iSite}(2,:)
%   PlotOpt.title = ['phiYX   :  ',LocSite,' : ' ...
%           num2str(YearOne) ' : ',...
%           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
%   PlotOpt.ColorAxisLabel = 'Degrees'
 %  PlotOpt.ColorAxisLabelPos = [3,.5];
%   if(subMed)
%      PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
%  end
%   plotFile = [PlotDir LocSite 'phiYX_' sM  ...
%                num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
%   [h,tMed_phiYX,fMed_phiYX] = ftPlot(Z2x2all,FT,phiYX,PlotOpt);
%   if EQLINE
%      plotEQlog
%   end
%   if SAVEplots
%       plot2File(plottype, plotFile);
%   end

%  if(iSite ==1) & plotBoth
%   %  plot for second set of dipoles at PKD
%   %   plot rho
%   PlotOpt.lims = rhoLims{3}(1,:);
%   PlotOpt.title = ['rhoXY2   :  ',LocSite,' : ' ...
%           num2str(YearOne) ' : ',...
%           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
%   plotFile = [PlotDir LocSite 'rhoXY2_' sM ...
%                num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
%   if subMed
%      PlotOpt.ColorAxisLabel = 'Fractional Deviation'
%      PlotOpt.ColorAxisLabelPos = [5,.5];
%      PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
%   else
%      PlotOpt.ColorAxisLabel = 'log_{1} ohm-m'
%      PlotOpt.ColorAxisLabelPos = [3.5,.5];
%   end
%   [h,tMed_rhoXY2,fMed_rhoXY2] = ftPlot(Z2x2all,FT,rhoXY2,PlotOpt);
%   if EQLINE
%      plotEQlog
%   end
%   if SAVEplots
%       plot2File(plottype, plotFile);
%   end

%   PlotOpt.lims = rhoLims{3}(2,:);
%   PlotOpt.title = ['rhoYX2   :  ',LocSite,' : ' ...
%           num2str(YearOne) ' : ',...
%           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
%   plotFile = [PlotDir LocSite 'rhoYX2_' sM ...
%                num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
%   if subMed
%      PlotOpt.ColorAxisLabel = 'Fractional Deviation'
%     PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
 %     PlotOpt.ColorAxisLabelPos = [5,.5];
 %  else
%      PlotOpt.ColorAxisLabel = 'log_{1} ohm-m'
%      PlotOpt.ColorAxisLabelPos = [3.5,.5];
%   end
%   [h,tMed_rhoYX2,fMed_rhoYX2] = ftPlot(Z2x2all,FT,rhoYX2,PlotOpt);
%   if EQLINE
%      plotEQlog
%   end
 %  if SAVEplots
%       plot2File(plottype, plotFile);
%   end

   %%   plot phi
%   PlotOpt.lims = phiLims{iSite}(1,:)
%   PlotOpt.title = ['phiXY2   :  ',LocSite,' : ' ...
%           num2str(YearOne) ' : ',...
%           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
%   plotFile = [PlotDir LocSite 'phiXY2_' sM ...
%                num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
%   PlotOpt.ColorAxisLabel = 'Degrees'
%   PlotOpt.ColorAxisLabelPos = [3,.5];
%   if(subMed)
%     PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
%   end
%   [h,tMed_phiXY2,fMed_phiXY2] = ftPlot(Z2x2all,FT,phiXY2,PlotOpt);
%   if EQLINE
%      plotEQlog
%   end
%   if SAVEplots
%       plot2File(plottype, plotFile);
%   end

%   PlotOpt.lims = phiLims{iSite}(2,:)
%   PlotOpt.title = ['phiYX2   :  ',LocSite,' : ' ...
%           num2str(YearOne) ' : ',...
%           num2str(dayOne) '-' num2str(dayEnd) ' :' root ]; 
%   plotFile = [PlotDir LocSite 'phiYX2_' sM ...
%                num2str(dayOne) '-' num2str(dayEnd) root '.' plottype];
%   PlotOpt.ColorAxisLabel = 'Degrees'
%   PlotOpt.ColorAxisLabelPos = [3,.5];
%   if(subMed)
%      PlotOpt.title = [ PlotOpt.title ' : Deviation from Median'];
%   end
%   [h,tMed_phiYX2,fMed_phiYX2] = ftPlot(Z2x2all,FT,phiYX2,PlotOpt);
%   if EQLINE
%      plotEQlog
%   end
%   if SAVEplots
%       plot2File(plottype, plotFile);
%   end
% end
