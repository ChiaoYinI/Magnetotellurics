%   root PKDSAO data directory
PKDSAOdata = setDataDir('PKDSAOdata');

%add for use if your ebvironment variables are configured for your own dirs
% PKDSAOdata= '/home/gauss/scratch/PKDSAOg/data/';

%  set ONLY8 = 1 to restrict analysis to 2 horizontal H, 
% and 2 E at each site (8 channels total)
OMIT_CH = 1;
%USE_CH = [1 2 4 5 6 7 8 9 11 12 ];
USE_CH = [1 2 3 4 5 6 7 8 9 10 11 12 ];
USE_BAND = [1:20,22:25];
OMIT_CH = 0;
grouping = 'standard';
% ONLY8 = 1;
%  USE_CH = [ 1 2 4 5 8 9 11 12];

%   define total time window to plot
%  NOTE: YEAR is now given as a character string, denoting the
%   directory to work in for ploting, not a numerical year 
YEAR = '2002';
DAYstart = 1;
DAYend = 365;
%  Number of days to put in one plot (only data for one page is loaded or
%   saved at one time)
% This is default ... full range on one page
DaysPerPlot = DAYend-DAYstart+1;

%  and sampling band
BAND = 1;    %   1 or 40 (hz)

%  SAVE = 1 to output SDM data structure for all segments
SAVE = 1;
%  set LOAD = 1 to load an already saved SDM structure
%   need to run with SAVE = 1, LOAD = 0 for same time window, band before
%   calling with LOAD = 1
%  directory for to save plot results to
SDMplotDir = 'SDMplots/';
plotDIR = [PKDSAOdata YEAR '/' SDMplotDir];

LOAD = 0;

%    frac defines fraction of data that must be present to plot SDM
%   results (this eliminates SDMs for which most data were missing)
frac = .4;

%   set SAVEplots = 1 to print figures to files in eps format;
%   file names are generated automatically
SAVEplots = 1;
EQLINE = 0;

%  ROOTS is a list of all plots that can be made at present ...
%  mkPLOT is an indicator array which defines which plots to actually make
%   ... just make them all
mkPLOT = [1 1 1 0 0 1 0 0 0 0];

%   set up TFTYPE for coherent nosie/residual estiamteoin
setTFTYPE
IER = 0;

%  new plots can be added to this list, and instructions for computations
%  what to plot, labels, etc., added to FTsdmPlotCases
%   See FTsdmPlotCases for definitions of each of the available cases
PLOTS = {'SDM Eigenvalues',...
	'Magnetic Field SNR',...
	'Electric Field SNR',...
	'Canonical Coherences: PKD vs. SAO',...
	'Canonical Coherences: H vs. E', ...
	'Average Modes', ...
	'Coherent Noise',...
        'Coherent Between Sites',...
        'Residuals: PKD',...
        'Residuals: SAO'};

binByHours = 0;

%   ROOTS defines plot file name for each plot type; if new plot types
%   are added to FTsdmPlotCases (and PLOTS list), add to this list also
ROOTS = {'EVAL','MSNR','ESNR','CCPS','CCHE','MODE','CohN','CohS','ResP','ResS'};
SdmAvgComputed = 0;

%  loop over pages; first set beginging and ending day for each page
for dayOne = DAYstart:DaysPerPlot:DAYend-1
   dayEnd = dayOne+DaysPerPlot-1';
   if LOAD
       %  construct default file names, using same rules that SDM merge uses
       %   for output
       DataDir = [PKDSAOdata YEAR '/SDM/'];
       if BAND == 40
          cfile = [DataDir 'S0_' num2str(dayOne) '_' num2str(dayEnd) '_40Hz.mat'];
       else
         cfile = [DataDir 'S0_' num2str(dayOne) '_' num2str(dayEnd) '_1Hz.mat'];
       end
       eval(['load ' cfile])
    else
       %  SDMmerge returns cell arrays SDMS SDMHD containing all available
       %   sdm processing results within selected time window.
       [SDMS,SDMHD,TimeInDays,good] = SDMmerge(dayOne,dayEnd,YEAR,BAND,SAVE);
    end

    % set "good" indicator to 0 for segments with too few data
    %   (this uses frac, defined above)
    nFiles = length(SDMS);
    NFS = [];
    for k = 1:nFiles
       if good(k)
          NFS = [NFS SDMS{k}.nf(1)];
       else
          NFS = [NFS 0];
       end
    end
    NFmin = max(NFS)*frac;
    good(NFS<NFmin) = 0;

    if OMIT_CH
       % subsample the 8 main channels from the full SDM
       %   for all "good" SDMs   ... can change which 8 channels
       %   to use here, but we assume later that these are Hx, Hy, Ex, Ey
       %   first for PKD, then for SAO.  (could change which set of
       %    dipoles are used for Ex/Ey at PKD)
       for k = 1:nFiles
          if good(k)
             [SDMS{k},SDMHD{k},IER] = SDM_ch_subset(SDMS{k},SDMHD{k},...
			USE_CH,USE_BAND,grouping);
             if(IER < 0)
                ['Error at file ', num2str(k)]
                break
             end
          end
       end
    end

    for k = 1:nFiles
       if(good(k))
          if(any(any(SDMS{k}.lambda(1:4)==0)) | ...
       	     any(any(SDMS{k}.var==0)) | any(any(SDMS{k}.Sig==0)))
             good(k) = 0;
          end
       end
    end
  
    %   compute averaged SDM
    %   Use "really good" data for initial estimate of averaged sdm
    %     for 2004 days 155-300 (less bad days).
    %    day numbers to use should be in VGind ... those in this
    %     list flagged already as "bad" (i.e., with good = 0) are not
    %     used
    load SDMavgDays
    veryGood = zeros(size(good));
    veryGood(VGind) = 1;
    veryGood = veryGood.*good;

    [SdmAvg] = SDMavg(SDMS,SDMHD,veryGood,grouping);

    fac = .10;
    nGood = sum(veryGood);
    [Nch,Nbt] = size(SdmAvg.var);
    SNRavg = SdmAvg.Sig./SdmAvg.var;
    lowSNR = ones(nGood,Nch,Nbt);
    kk = 0;
    for k = 1:nFiles
       if(veryGood(k))
          kk = kk + 1;
          SNR = (SDMS{k}.Sig./SDMS{k}.var);
          lowSNR(kk,:,:) = (SNR < fac*SNRavg);
          ChSkip = (sum(squeeze(lowSNR(kk,:,:)),2) > Nbt/3);
          if(any(ChSkip)) 
              %fprintf(1,'%d,%s',k,' bad')
              veryGood(k) = 0;
          end
       end
    end
    [SdmAvg] = SDMavg(SDMS,SDMHD,veryGood,grouping);
    SdmAvgComputed = 1;
    Sdms = SdmAvg;
    SdmHd = SdmAvg.Hd;
    cfile = ['SDMAVG_' num2str(dayOne) '-' num2str(dayEnd) '.mat'];
    eval(['save ' cfile ' Sdms SdmHd']);
    
    %   Period for each band in the SDM files .... this assumes
    %   that all time segments are processed with the same
    %   time windows and period bands (this is checked in SDMmerge)
    k1 = min(find(good));
    Period = SDMS{k1}.T;

    nPLOTS = length(PLOTS);
    %  loop over plots
    for kPlot = 1:nPLOTS
       if mkPLOT(kPlot)
          if kPlot >= 7  
             %  first call TFsdmPlotCases for calculation of coherent
             %    noise and residuals
             plotType = 'Coherent Noise Setup';
             FTsdmPlotCases
          end
          %   FTsdmPlotCases does all of the computing and plot set up
          %   for coherent noise cases there is a second call to FTsdmPlotCases
          plotType = PLOTS{kPlot}
          FTsdmPlotCases
          %   multiplot does the actual plotting
          if BAND == 1
             for k = 1:length(OPTIONS)
                OPTIONS{k}.DayLine = 0;
             end
          end
          multiPlot(D,TimeInDays,Period,OPTIONS);
          if EQLINE
              plotEQ
          end
          if SAVEplots
             %  print figure to eps file
             cfile = [ROOTS{kPlot} '_' num2str(dayOne) '-' num2str(dayEnd) '.eps'];
             eval(['print -depsc ' cfile]);
%             cfile = [plotDIR ROOTS{kPlot} '_' num2str(dayOne) '-' num2str(dayEnd) '.jpg'];
%             eval(['print -djpeg90 ' cfile]);
          end
       end
    end
end
