function resPlt1(DATE,site,EH,id)
%  Usage:  resPlt1(DATE,site,EH,id)
%    Inputs: DATE structure (wtih fields .year .day ), 
%            site (PKD or SAO),
%            EH  = 'E' or 'H'
%            id = decimation level
%   plots signal amplitiude and normalized residuals for a 1 day segment
%   given by DATE.year DATE.day, for the E or H channels at one site
%     This script is a simple driver, with many parameters hard-coded

%   edit directory: $PKDSAOdata

PKDSAOdata = '/home/server/pi/homes/egbert/PKDSAOg/data/';
dir = [PKDSAOdata num2str(DATE.year)];

%   plot raw unaveraged FC amplitudes; no averaging
navg = 1;
%  might need to change this
iband = [1:32; 1:32];

DayStr = num2str(DATE.day);
if DATE.day < 10
   DayStr = ['00' DayStr];
elseif DATE.day < 100
   DayStr = ['0' DayStr];
end

%   find file names from DATE: assumes standard PKD/SAO array from later years
FCfiles{1} = ['PKD_' DayStr '.f7'];
FCfiles{2} = ['SAO_' DayStr '.f5'];

%  check to see that both files are in FC directory
cfile = [dir '/FC/' FCfiles{1} ];
fid1 = fopen(cfile,'r');
fclose(fid1);
cfile = [dir '/FC/' FCfiles{2} ];
fid2 = fopen(cfile,'r');
fclose(fid2);
if (fid1 < 0 | fid2 < 0 )
   fprintf(1,'Error: missing FC files')
   return
end

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
        clSig = [-1,-1,-1.5;-4,-4,-4.5];
        clRes = [-1,-1,-1.5;-4,-4,-4.5]/2;
    else
        Kout = [11 12];
        ctitle = 'Electric fields at Hollister';
        clSig = [-2;1]*ones(1,2);
        clRes = [-2;1]*ones(1,2)/2;
    end
end
    
OPTIONS = struct('dir',dir,'Kin',Kin,'Kout',Kout,...
    'navg',navg,'id',id,'iband',iband,'plotAmp',1,'plotRes',1,...
    'NormalizeSig',0,'NormalizeRes',1,'title',ctitle,...
    'clSig',clSig,'clRes',clRes,'linT',1,...
    'clSigNorm',[-10,10],'clResNorm',[-15,0]);

SigResAmp(FCfiles,OPTIONS);
