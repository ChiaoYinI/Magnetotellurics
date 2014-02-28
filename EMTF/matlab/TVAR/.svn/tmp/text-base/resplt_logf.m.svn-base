function [SigAll,ResAll] = resplt_logf(DATE,site,EH,bsFile,Pl)
%  Usage:  resplt_logf(DATE,site,EH,bsFile)
%    Inputs: DATE structure (wtih fields .year .day, 
%            site (PKD or SAO),
%            EH  = 'E' or 'H'
%            bsFile = band set up file
%            Pl : optional argument ... = 1 (default) for plotting
%   Plots signal amplitiude and normalized residuals for one 
%   day of 2 hour time segments, averaged into logarithmically
%   spaced frequency bands/decimation levels specified in bsFile
%
%   This script is a simple driver, with many parameters hard-coded
%
%   Before use set : PKDSAOdata to root PKDSAO data directory

if nargin < 5
    Pl = 1;
end
PKDSAOdata = '/home/server/pi/homes/egbert/PKDSAOg/data/';
dir = [PKDSAOdata num2str(DATE.year)];

%  get averaging bands and decimation levels from band setup file
fid = fopen(bsFile,'r');
fline = fgetl(fid);
nBand = sscanf(fline,'%d');
ID = zeros(1,nBand);
iband = zeros(2,nBand);
for ib = 1:nBand
   fline = fgetl(fid);
   temp = sscanf(fline,'%d',3);
   ID(ib) = temp(1);
   IBAND(:,ib) = temp(2:3);
end
fclose(fid);
nd = max(ID);
NAVG = [16 4 1 1];

%  set up channels for site, EH
if site == 'PKD'
    %  use only SAO for prediction?
    Kin = [8 9 11 12];
    %  or use all channels?  
    % Kin = [1:12];
    if EH == 'H'
        Kout = [1 2 3];
        ctitle = 'Magnetic fields at Parkfield';
        clSig = [-3,-3,-3.5;0,0,-.5];
        clRes = [-4,-4,-4.5;-1,-1,-1.5];
    else
        Kout = [4 5 6 7];
        ctitle = 'Electric fields at Parkfield';
        clSig = [-2.5;.5]*ones(1,4);
        clRes = [-2.5;.5]*ones(1,4);
    end
else
    % use only PKD for prediction
    Kin = [1:7];
    if EH == 'H'
        Kout = [8 9 10];
        ctitle = 'Magnetic fields at Hollister';
        clSig = [-3,-3,-3.5;0,0,-.5];
        clRes = [-4,-4,-4.5;-1,-1,-1.5];
    else
        Kout = [11 12];
        ctitle = 'Electric fields at Hollister';
        clSig = [-2;1]*ones(1,4);
        clRes = [-2;1]*ones(1,4);
    end
end

%  loop over decimation levels
freq = [];
avgAmpSig = [];
avgAmpRes = [];
for id = 1:nd
   idL1 = min(find(ID==id)); 
   idL2 = max(find(ID==id)); 
   navg = NAVG(id);
   iband = IBAND(:,idL1:idL2);

   %  set options for residual averaging
   OPTIONS = struct('dir',dir,'Kin',Kin,'Kout',Kout,...
    'navg',navg,'id',id,'iband',iband,'plotAmp',0,'plotRes',0);

   % loop over hours, storing time averaged residual and signal powers (no plotting
   %    of individual 2 hr sections
   k = 0;
   k1 = 1;
   good = ones(12,1);
   notInitialized = 1;
   for hr = 0:2:22
       k = k+1;

       % find file names from DATE: assumes standard PKD/SAO array from later years
       DayStr = num2str(DATE.day);
       if DATE.day < 10
          DayStr = ['00' DayStr];
       elseif DATE.day < 100
          DayStr = ['0' DayStr];
       end
       HRstr = num2str(hr);
       if (hr < 10)
          HRstr = ['0' HRstr];
       end

       FCfiles{1} = ['PKD_' DayStr '_' HRstr '.f7'];
       FCfiles{2} = ['SAO_' DayStr '_' HRstr '.f5'];

       % check to see that both files are in FC directory
       cfile1 = [dir '/FC/' FCfiles{1} ];
       fid1 = fopen(cfile1,'r');
       if(fid1 > 0)
          fclose(fid1);
       end
       cfile2 = [dir '/FC/' FCfiles{2} ];
       fid2 = fopen(cfile2,'r');
       if(fid2 > 0)
          fclose(fid2);
       end
       if (fid1 < 0 | fid2 < 0 )
          good(k) = 0;
          fprintf(1,'%s\n',['FC files '  FCfiles{1} ...
		', ' FCfiles{2} ' not found'])
       else
          [Sig,Res] = SigResAmp(FCfiles,OPTIONS);
          if ~isstruct(Sig)
             fprintf(1,'%s \n','Error reading FC files')
             good(k) = 0;
          else
             [nch,nT,nF] = size(Sig.data);
             if notInitialized
              %   initialize structures for the full day of signal and residual powers
                 if id == 1
                    SigAll = Sig;
                    ResAll = Res;
                 end
                 s = zeros(nch,12*nT,nF);
                 t = [];
                 %  is this needed?
                 %s = s./s;
                 r = s;
                 freq = [freq Sig.f];
                 avgAmpSig = [avgAmpSig squeeze(Sig.avgAmp)];
                 avgAmpRes = [avgAmpRes squeeze(Res.avgAmp)];
                 notInitialized = 0;
              end
              k2 = k1 + nT-1;
              s(:,k1:k2,:) = Sig.data;
              r(:,k1:k2,:) = Res.data;
              t = [t Sig.t];
              k1 = k2+1;
           end
       end
    end
    if id == 1
       % copy accumultated time + signal and residual amplitudes 
       %       into SigAll, ResAll
       time = t;
       SigAll.data = s(:,1:k2,:);
       SigAll.t = t;
       ResAll.data = r(:,1:k2,:);
       ResAll.t = t;
    else
       %  interpolate signal and residual amplitudes to decimation
       %     level 1 times
       nT1 = length(time);
       temp = zeros(nch,nT1,idL2);
       temp(:,:,1:idL1-1) = SigAll.data(:,:,:);
       for ich = 1:nch
          temp(ich,:,idL1:idL2) = ...
  		interp1(t,squeeze(s(ich,1:k2,:)),time,'nearest'); 
       end
       SigAll.data = temp;
       temp = zeros(nch,nT1,idL2);
       temp(:,:,1:idL1-1) = ResAll.data(:,:,:);
       for ich = 1:nch
          temp(ich,:,idL1:idL2) = ...
		interp1(t,squeeze(r(ich,1:k2,:)),time,'nearest'); 
       end
       ResAll.data = temp;
    end
end       
SigAll.f = freq;
ResAll.f = freq;
SigAll.avgAmp = avgAmpSig;
ResAll.avgAmp = avgAmpRes;
   
if Pl
   OPTIONS = struct('dir',dir,'Kin',Kin,'Kout',Kout,...
    'navg',navg,'id',id,'iband',iband,'plotAmp',1,'plotRes',1,...
    'NormalizeSig',0,'NormalizeRes',1,'title',ctitle,...
    'clRes',clRes,'clSig',clSig,'linT',0,...
    'clSigNorm',[-.5,.5],'clResNorm',[-30,0]);

   SigResPlot(SigAll,ResAll,OPTIONS);
end
