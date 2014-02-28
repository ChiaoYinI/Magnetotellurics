function [SigDay,ResDay] = resPlt1hz10Day(DATE,site,EH,id,Pl)
%  Usage:  resplt1hz10Day(DATE,site,EH,id,Pl)
%          [SigDay,ResDay] =  resplt1hz10Day(DATE,site,EH,id,Pl)
%    Inputs: DATE structure (wtih fields .year .day, 
%            site (PKD or SAO),
%            EH  = 'E' or 'H'
%            id = decimation level
%            Pl : optional argument ... = 1 (default) for plotting
%   Plots signal amplitiude and normalized residuals for 10 days
%    of 1hz data processed in 1 day segments
%      for the E or H channels at one site
%   This script is a simple driver, with many parameters hard-coded

%   edit directory: $PKDSAOdata

%  could make number of days any input option ...
%    ... hard coded here for now
Ndays = 10;
if nargin < 5
    Pl = 1;
end
PKDSAOdata = '/home/server/pi/homes/egbert/PKDSAOg/data/';
dir = [PKDSAOdata num2str(DATE.year)];

%   just average over time ... no frequency band averaging
id = min(id,3);
NAVG = [16 4 1];
navg = NAVG(id);
iband = [1:32;1:32];

%  set up channels for site, EH
if site == 'PKD'
    %  use only SAO for prediction?
    Kin = [8 9 11 12];
    %  or use all channels?  
    % Kin = [1:12];
    if EH == 'H'
        Kout = [1 2 3];
        ctitle = 'Magnetic fields at Parkfield';
        clSig = [-3,-3,-3.5;1,1,.5];
        clRes = [-3,-3,-3.5;1,1,.5]/2;
    else
        Kout = [4 5 6 7];
        ctitle = 'Electric fields at Parkfield';
        clSig = [-2.5;.5]*ones(1,4);
        clRes = [-2.5;.5]*ones(1,4)/2;
    end
else
    % use only PKD for prediction
    Kin = [1:7];
    if EH == 'H'
        Kout = [8 9 10];
        ctitle = 'Magnetic fields at Hollister';
        clSig = [-4,-4,-4.5;-1,-1,-1.5];
        clRes = [-4,-4,-4.5;-1,-1,-1.5]/2;
    else
        Kout = [11 12];
        ctitle = 'Electric fields at Hollister';
        clSig = [-2;1]*ones(1,4);
        clRes = [-2;1]*ones(1,4)/2;
    end
end
    
OPTIONS = struct('dir',dir,'Kin',Kin,'Kout',Kout,...
    'navg',navg,'id',id,'iband',iband,'plotAmp',0,'plotRes',0)

%   loop over days, storing time averaged residual and 
%      signal powers (no plotting of individual processed segments
k = 0;
k1 = 1;
good = ones(Ndays,1);
notInitialized = 1;
for day = 1:Ndays
    k = k+1;

   %   find file names from DATE: assumes standard PKD/SAO array from later years
   DAY = DATE.day+day-1;
   DayStr = num2str(DAY);
   if  DAY < 10
      DayStr = ['00' DayStr];
   elseif DAY < 100
      DayStr = ['0' DayStr];
   end

   FCfiles{1} = ['PKD_' DayStr '.f7'];
   FCfiles{2} = ['SAO_' DayStr '.f5'];

   %  check to see that both files are in FC directory
   cfile = [dir '/FC/' FCfiles{1} ];
   fid1 = fopen(cfile,'r');
   cfile = [dir '/FC/' FCfiles{2} ];
   fid2 = fopen(cfile,'r');
   if (fid1 < 0 | fid2 < 0 )
      good(k) = 0;
   else
      [Sig,Res] = SigResAmp(FCfiles,OPTIONS);
      if ~isstruct(Sig)
         % fprintf(1,'%s','Error reading TF file')
         good(k) = 0;
      else
         [nch,nT,nF] = size(Sig.data);
         if notInitialized
          %   initialize structures for the full day of signal and residual powers
             SigDay = Sig;
             ResDay = Res;   
             s = zeros(nch,12*nT,nF);
             t = [];
             s = s./s;
             r = s;
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
SigDay.data = s(:,1:k2,:);
SigDay.t = t;
ResDay.data = r(:,1:k2,:);
ResDay.t = t;

if Pl
   OPTIONS = struct('dir',dir,'Kin',Kin,'Kout',Kout,...
    'navg',navg,'id',id,'iband',iband,'plotAmp',1,'plotRes',1,...
    'NormalizeSig',0,'NormalizeRes',1,'title',ctitle,...
    'clRes',clRes,'clSig',clSig,'linT',1,...
    'clSigNorm',[-10,10],'clResNorm',[-15,0]);
   SigResPlot(SigDay,ResDay,OPTIONS);
end
