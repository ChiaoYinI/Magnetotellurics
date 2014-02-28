function [FT] = FTsetup(FileNamesFC,id,ifreq,csta)
%  Usage: [FT] = FTsetup(FileNamesFC,id,ifreq,csta);
%   returns Fourier coefficients for array defined by
%   the list of FC files in cell array FileNamesFC for decimation
%   level id
%
%   CALLS:   fc_open.m, mk_start_freqs.m, mk_isets.m, fc_get.m
%    which depend on: unpack.m 

ids = [];start_decs = [];  nch = [];orient = [];fids = [];chid = [];
sta = [];
nsta = length(FileNamesFC);
%  NOTE: fracOffset is fractional offset between sets ...
%    this is required to correctly convert set number to time,
%    and is not in the fc files (should be!)
%   the value hard-coded here is the what has been used for
%   essentially all previous use of EMTF/dnff
fracOffset = 0.75;

%  open FC files, read headers
for ista = 1:nsta
  %  fprintf(1,'%s %s \n','File ',FileNamesFC{ista})
  [fid,nd,nf,nch1,chid1,orient1,drs,stdec,decs,start_dec] ...
                    = fc_open(FileNamesFC{ista});
  if (fid < 0 )
     fprintf(1,'%s \n','File not found');
     FT = [];
     return ;
  end

  fids = [ fids fid ];
  start_decs = [start_decs ; start_dec ] ;
  chid = [chid; chid1];
  [n1,n2] = size(chid1);
  sta1 = reshape(blanks(n1*n2),n1,n2);
  nChar = length(char(csta{ista}));
  for ich = 1:nch1
     sta1(ich,1:nChar) = char(csta{ista});
  end
  sta = [sta ; char(sta1)];
  nch = [nch  nch1 ];
  orient = [orient orient1(1,:) ];
  if(ista == 1)
       nfreq = max(ifreq);
       freq = ifreq/(decs(2,id)*drs(id)); 
       offset = (fracOffset*decs(2,id));
  end
end
ich = 0;
nt = sum(nch);

%   set up array of pointers to start of each frequency in each FC file
[start_freqs] = mk_start_freqs(id,fids,nfreq,nch,start_decs);

%   find set numbers available for all sites
[isets,isets_pt] =  mk_isets(fids,start_decs,nch,nfreq,id);
nsets = length(isets);

%  time expressed in days (1 Jan 0:00 UT = 1.0 )
%   NOTE:  The 2/3 centers the time mark on the FT window
%    ... ASSUMING the offset between windows is 0.75 of window length 
%   NOTE: dnff asigns set # 1 to window beginging at initial time
%     (1 Jan 0:00 UT for PKD/SAO processing)
fac = 1-(1-fracOffset)/fracOffset;
dt = drs(id)*offset/86400;
time = 1+dt*(isets(1:nsets)+fac-1);

%  Fourier coefficients
[fc] = fc_get(fids,start_freqs,nch,isets_pt,id,ifreq);
for k = 1:length(fids)
   fclose(fids(k));
end

% at present we only allow for the case where both data and
%   TF are in physical units ... the actual system TF is not
%   stored in field sysTF, but the SysTF.PhysU flag is set
%   BUT: see SysTFsetup.m
SysTF.PhysU = 1;
FT = struct('data',fc,'t',time,'dt',dt,'f',freq,'ch_id',chid,...
	'sta',sta,'nch',nt,'SysTF',SysTF,'amp',0);
