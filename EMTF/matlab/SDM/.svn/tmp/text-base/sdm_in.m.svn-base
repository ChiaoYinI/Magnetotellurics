%************************************************************************
%uev_in : Reads in sdm etc. for a single band (S0**** file)
% Usage: [period,nf,var,S,err] = sdm_in(fid_uev,nt,ib,irecl)
% Inputs:
%     fid_uev = uev file id
%     nt, = number of components
%     ib = period band sought
%     irecl = band record length
% Outputs:
%     period = actual period
%     nf = # of data vectors
%     ev = eigenvalues
%     var = error variances
%     u = eigenvectors
%     err = -1 for read error

function [period,nf,var,S,err] =uev_in(fid,nt,ib,irecl)

err = 0;
status = fseek(fid,ib*irecl,'bof');
[period,count] = fread(fid,1,'float');
if count == 0
   err = -1;
   return
end
 
[nf,count] = fread(fid,1,'long');
if count == 0
   err = -1;
   return
end
ns = (nt*(nt+1))/2;
[var,count] = fread(fid,nt,'float');
if count < nt
   err = -1;
   return
end
[s,count] = fread(fid,[2,ns],'float');
if count < 2*ns
   err = -1;
   return
end
s = s(1,:)+i*s(2,:);
%   convert S to full storage mode ... save programming
%    hassles later
col = [];
row = [];
for k=1:nt
  row = [ row , k*ones(1,k) ]; col = [ col  , [1:k]];
end
S = sparse(row,col,s) ;
S = full(S);  S =  S + S' - diag(diag(S));

%  Record for each band written by this
%      subroutine wrt_uev(iouev,irec,s,ns,var,nt,nf,period)
%      complex s(ns)
%      real var(nt),period
%      integer nf
%      write(iouev,rec=irec) period,nf,var,s
%      return
%      end
