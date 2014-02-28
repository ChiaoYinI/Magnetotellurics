%************************************************************************
%Usage:  [fid,recl,nbt,nt,nsta,nsig,nch,ih,stcor,...
%              decl,chid,csta,sta,orient] = Pw_hd(cfile,system);
%   Input:   cfile = file name
%		system = 'U' if PW files were written on Unix, 'P'
%			if on PC; optional
%   Returns: fid_uev = file id
%            recl = record length
%            nbt,nt,nsta,nsig = # of bands, components, stations, evecs
%            nch(nsta),ih(nsta+1),stcor(2,nsta),decl(nsta),sta(nsta),
%            orient(nsta)  = the usual
%
%*********************************************************************

function [fid,recl,nbt,nt,nsta,nsig,nch,ih,stcor,...
              decl,chid,csta,sta,orient]= Pw_hd(cfile,system)
%  'b' option in fopen assumes binary Pw-file was written on UNIX system
%  change to 'l' for files written on PC
lCOMP=computer; 
fid = fopen(cfile,'r'); % native
irecl = fread(fid,1,'long'); % =16 if right written fmt==native
fclose(fid);
opt='n';
if irecl~=16,
     if lCOMP(1:2)=='PC',
      warndlg(['SDM file written on Unix, read on PC']);
      opt='b';% files are in big endian, read on PC
     else
      warndlg(['SDM file written on PC, read on Unix']);
      opt='l';% files are in little endian, read in Unix
     end                     
end
fid = fopen(cfile,'r',opt);
%
fseek(fid,4,'cof');
nt    = fread(fid,1,'long');
nsta  = fread(fid,1,'long');
nsig  = fread(fid,1,'long');
nbt   = fread(fid,1,'long');

for k=1:nsta
  fseek(fid,4,'cof');
  nbytes = fread(fid,1,'long');
  csta_length = nbytes - 20;
  nch(k) = fread(fid,1,'long');
  ih(k)  = fread(fid,1,'long');
  stcor(1:2,k) = fread(fid,2,'float');
  decl(k)  = fread(fid,1,'float');
  sta(1:csta_length,k)   = fread(fid,csta_length,'char');
end
fseek(fid,4,'cof');
nbytes = fread(fid,1,'long');
chid_length = nbytes - 8 - csta_length;
for l = 1:nt
   orient(1:2,l) = fread(fid,2,'float');
   chid(1:chid_length,l) = fread(fid,chid_length,'char');
   csta(1:csta_length,l) = fread(fid,csta_length,'char');
   if(l < nt) fseek(fid,8,'cof'); end
end
%recl(1) = length of header block ...
recl(1) = 24+nsta*(28+csta_length)+nt*(16+chid_length+csta_length);
ns = nt*(nt+1)/2;
recl(2) = 64+8*(ns+2*nt);
