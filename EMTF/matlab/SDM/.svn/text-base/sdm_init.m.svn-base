%************************************************************************
%sdm_init : %Initializes S# file containing SDMs from multiple station program
%Usage:  [fid,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,...
%	decl,sta,chid,csta,orient,periods,ndf] = sdm_init(cfile,system);
%   Input:   cfile = file name
%   Returns: fid = file id
%            irecl = record length
%            nbt,nt,nsta,nsig = # of bands, components, stations, evecs
%            nch(nsta),ih(nsta+1),stcor(2,nsta),decl(nsta),sta(nsta),
%            orient(nsta)  = the usual
%            periods(nbt)  = periods
%            ndf(nbt)  = number of FC vectors averaged for each band

%   fid is negative if there is any read error

function [fid,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,...
         decl,sta,chid,csta,orient,periods,ndf] = sdm_init(cfile,system)
%  'b' option in fopen assumes binary Pw-file was written on UNIX system
%  change to 'l' for files written on PC
if nargin < 2
% try to define system files were written on
   lCOMP=computer; 
   fid = fopen(cfile,'r'); % native
   if fid < 0
      return
   end
   [irecl,count] = fread(fid,1,'long');
   if count == 0
      fclose(fid);
      fid = -fid;
      return
   end
     
   [a,count]=fread(fid,3,'long');
   if count < 3
      fclose(fid);
      fid = -fid;
      return
   else
      nsta=a(3);
   end
   fclose(fid);
   opt='n';
   if irecl<0 | nsta<0 | nsta>100,
       if lCOMP(1:2)=='PC',
           opt='b';% files are in big endian, read on PC
       else
           opt='l';% files are in little endian, read in Unix
       end                     
   end
   fid = fopen(cfile,'r',opt);
   if fid < 0
      return
   end
else
   if system == 'U'
      fid = fopen(cfile,'r','b');
      if fid < 0
         return
      end
   else
      fid = fopen(cfile,'r','l');
      if fid < 0
         return
      end
   end
end
[irecl,count] = fread(fid,1,'long');
if count == 0
   fclose(fid);
   fid = -fid;
   return
end
[nbt,count]   = fread(fid,1,'long');
if count == 0
   fclose(fid);
   fid = -fid;
   return
end
[nt,count]    = fread(fid,1,'long');
if count == 0
   fclose(fid);
   fid = -fid;
   return
end
[nsta,count]  = fread(fid,1,'long');
if count == 0
   fclose(fid);
   fid = -fid;
   return
end
[nsig,count]  = fread(fid,1,'long');
if count == 0
   fclose(fid);
   fid = -fid;
   return
end
csta = [];
%   First assume this is a "newest" version sdm file: 4 character station ids
%    and 6 character channel ids
nCsta = 4;
sta = zeros(nCsta,nsta);
csta = [];
for k=1:nsta
  [temp,count] = fread(fid,1,'long');
  if count == 0
     fclose(fid);
     fid = -fid;
     return
  else
     nch(k) = temp;
  end
  [temp,count]  = fread(fid,1,'long');
  if count == 0
     fclose(fid);
     fid = -fid;
     return
  else
     ih(k) = temp;
  end
  [temp,count] = fread(fid,2,'float');
  if count < 2
     fclose(fid);
     fid = -fid;
     return
  else
     stcor(1:2,k) = temp;
  end
  [temp,count]  = fread(fid,1,'float');
  if count == 0
     fclose(fid);
     fid = -fid;
     return
  else
     decl(k) = temp;
  end
  [temp,count]   = fread(fid,nCsta,'char');
   if count == 0
      fclose(fid);
      fid = -fid;
      return
   else
      sta(1:nCsta,k) = temp;
   end
  for l=1:nch(k)
    csta = [csta [sta(1:nCsta,k)]];
  end
end
% Read first two channel orientations ...
[temp,count] = fread(fid,2,'float');
if count == 0
   fclose(fid);
   fid = -fid;
   return
else
   orient(1:2,1) = temp;
end
%if (orient(2,1) ~= 0) 
%   this must be a 3 character station id file (NOT newest format)
%   fseek(fid,20,'bof');
%   nCsta = 3;
%   sta = zeros(nCsta,nsta);
%   csta = [];
%   for k=1:nsta
%     [temp,count] = fread(fid,1,'long');
%     if count == 0
%        fclose(fid);
%        fid = -fid;
%        return
%     else
%        nch(k) = temp;
%     end
%     [temp,count]  = fread(fid,1,'long');
%     if count == 0
%        fclose(fid);
%        fid = -fid;
%        return
%     else
%        ih(k) = temp;
%     end
%     [temp,count] = fread(fid,2,'float');
%     if count < 2
%        fclose(fid);
%        fid = -fid;
%        return
%     else
%        stcor(1:2,k) = temp;
%     end
%     [temp,count]  = fread(fid,1,'float');
%     if count == 0
%        fclose(fid);
%        fid = -fid;
%        return
%     else
%        decl(k) = temp;
%     end
%     [temp,count]   = fread(fid,nCsta,'char');
%     if count < nCsta
%        fclose(fid);
%        fid = -fid;
%        return
%     else
%        sta(1:nCsta,k) = temp;
%     end
%     for l=1:nch(k)
%       csta = [csta [sta(1:nCsta,k)]];
%     end
%   end
%   [temp,count] = fread(fid,2,'float');
%   if count < 2
%      fclose(fid);
%      fid = -fid;
%      return
%   else
%      orient(1:2,1) = temp;
%   end
%end
[temp,count] = fread(fid,6,'char');
if count == 0
   fclose(fid);
   fid = -fid;
   return
else
   chid(1:6,1) = temp;
end

%%%  Continue, assuming that this is at least a "NEW" sdm files 
%%%%     ... with 6 character channel IDs
for l = 2:nt
   [temp,count] = fread(fid,2,'float');
   if count < 2
      fclose(fid);
      fid = -fid;
      return
   else
      orient(1:2,l) = temp;
   end
   [temp,count] = fread(fid,6,'char');
   if count < 6
      fclose(fid);
      fid = -fid;
      return
   else
      chid(1:6,l) = temp;
   end
end
%  Check to see if all orientations are reasonable:
%if any( orient(2,:) ~= 0 )  
%    assume that this is an old file ... get rid of this if you
%   have no old sdm files, and channels with non=zero tilts
%fseek(fid,-nt*14,'cof')
%clear chid
%%%  This is for OLD sdm files ... 2 character channel IDs
%%   for l = 1:nt
%      [temp,count] = fread(fid,2,'float');
%      if count < 2
%         fclose(fid);
%         fid = -fid;
%         return
%      else
%         orient(1:2,l) = temp;
%      end
%      [temp,count] = fread(fid,2,'char');
%      if count < 2
%         fclose(fid);
%         fid = -fid;
%         return
%      else
%         chid(1:2,l) = temp;
%      end
%   end
%end

% go through file and get all of the periods ...
for ib = 1:nbt
   status = fseek(fid,irecl*ib,'bof');
   [temp,count] = fread(fid,1,'float') ;
   if count == 0
      fclose(fid);
      fid = -fid;
      return
   else
      periods(ib) = temp;
   end
   [temp,count] = fread(fid,1,'long') ;
   if count == 0
      fclose(fid);
      fid = -fid;
      return
   else
      ndf(ib) = temp;
   end
end
ndf = ndf';

%  File header written by this:
%       write(ivunits(ix),rec=1) irecl,nbt,nt,nstau,nsig,
%     &      (nchu(k),ih(k),stcor(1,k),stcor(2,k),decl(k),sta(k),
%     &      k=1,nstau),(orient(1,l),orient(2,l),l=1,nt)


%  Record for each band written by this
%      subroutine wrt_uev(iouev,irec,s,ns,var,nt,nf,period)
%      complex s(ns)
%      real var(nt),period
%      integer nf
%      write(iouev,rec=irec) period,nf,var,s
%      return
%      end
