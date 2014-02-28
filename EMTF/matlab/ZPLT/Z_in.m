%************************************************************************
%Z_in(cfile) : Reads in One Z_***** file  (TFs and signal/error cov.
%USAGE:  [z,sig_s,sig_e,periods,ndf,stcde,orient,nch,nche,nbt] = Z_in(cfile);
% nch = total # of channels ; nche = nch-2 = # of predicted channels ; nbt = # of bands
%   (NOTE: First two channels are always the "predictors"
% z(2,nche*nbt) = complex TFs
%   NOTE: Z(1,1:nche) corresponds to response to Hx sources for first band,
%         Z(2,1:nche) is Hy for for first band,
%         Z(1,nche+1:2*nche) corresponds to response to Hx sources for second band,
%         Z(2,nche+1:2*nche) is Hy for for second band,  ETC. ...
% sig_s(2,2*nbt) = complex inverse signal covariance
% sig_e(nche,nche*nbt) = complex residual error covariance
% stdec(3) = station coordinates, declination
% periods(nbt) = periods in seconds
% orient(2,nch) = orientation (degrees E of geomagnetic N) for each channel

function [z,sig_s,sig_e,periods,ndf,stdec,orient,nch,nche,nbt,chid,csta] = Z_in(cfile)
rnan=-999.;
fid = fopen(cfile,'r');
if(fid < 0) 
   'error: cannot open file ',cfile
end
% skip text on first four lines
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
cline = fgets(fid);
% skip text in coordinate line
cline = fgets(fid);
stdec = sscanf(cline,'coordinate %f %f declination %f');
cline = fgets(fid);
nchnbt = sscanf(cline,'number of channels %d number of frequencies %d');
nch = nchnbt(1);
nche = nch-2;
nbt = nchnbt(2);
cline = fgets(fid);
for k=1:nch
   idum = fscanf(fid,'%d',1);
   orient(1:2,k) = fscanf(fid,'%f',2);
   cline = fgets(fid);
   ind = find(cline ~= ' ');
   ind = ind(1:end-1);
   i1 = ind(1);
   ii = ind(diff(ind)-1>0);
   i2 = ii(1);
   csta(k,1:i2-i1+1) = cline(i1:i2);
   i2 = ind(end);
   i1 = ind(min(find(ind>ii(1))));
   chid(k,1:i2-i1+1) = cline(i1:i2);
end
cline = fgets(fid);
z = zeros(2,nche*nbt) + i*zeros(2,nche*nbt);
sig_e = zeros(nche,nche*nbt) + i*zeros(nche,nche*nbt);
sig_s = zeros(2,2*nbt) + i*zeros(2,2*nbt);
ndf = zeros(nbt,1);
for ib = 1:nbt
  cline = fgets(fid);
  tmp=sscanf(cline,'period : %f');
  if isempty(tmp)==0,
   periods(ib) = tmp;
  else
   periods(ib) = NaN;
  end
  cline = fgets(fid);
  tmp=sscanf(cline,'number of data point %d');
  if isempty(tmp)==0,
   ndf(ib) = tmp;
  else
   ndf(ib) = NaN;  
  end
  k1 = nche*(ib-1) + 1;
  k2 = nche*ib;
  cline = fgets(fid);
  ztemp = fscanf(fid,'%e',[4,nche]);
  ztemp(find(ztemp==rnan))=NaN;
  if size(ztemp)>=[4,nche],
    z(1:2,k1:k2) = ztemp(1:2:3,:)+i*ztemp(2:2:4,:);   
  else
   z(1:2,k1:k2) =   NaN; 
  end
  chead = fgets(fid);
  chead = fgets(fid);
  stemp = fscanf(fid,'%e',[2,3]);
  stemp(find(stemp==rnan))=NaN;
  ncht = 2;
  for k = 1:ncht
    for l = 1:ncht
      if(l < k ) 
        kl = (k*(k-1))/2+l;
        if size(stemp)>=[ncht,ncht*(ncht+1)/2],
         sig_s(k,2*(ib-1)+l) = stemp(1,kl)+i*stemp(2,kl);
        else
         sig_s(k,2*(ib-1)+l) = NaN; 
        end
      else
        kl = (l*(l-1))/2+k; % max is ncht*(ncht-1)/2+ncht
        if size(stemp)>=[ncht,ncht*(ncht+1)/2],
            sig_s(k,2*(ib-1)+l) = stemp(1,kl)-i*stemp(2,kl);
        else
            sig_s(k,2*(ib-1)+l) = NaN;
        end
      end
    end
  end
  chead = fgets(fid);
  chead = fgets(fid);
  nse = (nche*(nche+1))/2;
  stemp = fscanf(fid,'%e',[2,nse]);
  stemp(find(stemp==rnan))=NaN;

   for k = 1:nche
    for l = 1:nche
      if(l < k ) 
        kl = (k*(k-1))/2+l;
        if size(stemp)>=[2,nche*(nche+1)/2],
         sig_e(k,nche*(ib-1)+l) = stemp(1,kl)+i*stemp(2,kl);
        else
         sig_e(k,nche*(ib-1)+l) = NaN;   
        end
      else
        kl = (l*(l-1))/2+k;
        if size(stemp)>=[2,nche*(nche+1)/2],
            sig_e(k,nche*(ib-1)+l) = stemp(1,kl)-i*stemp(2,kl);
        else
            sig_e(k,nche*(ib-1)+l) = NaN;
        end
      end
    end
   end
  cline = fgets(fid);
end  
