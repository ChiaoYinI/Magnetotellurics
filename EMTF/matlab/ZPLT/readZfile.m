function [Z] = readZfile(cfile)

%  Usage: Z = readZfile(cfile);
%
%  reads standard Z file, returning everything in structure Z
%  calling writeZfile subsequently will result in essentially the
%  same file, in the same format (differences in real number formatting
%  conventions only)

rnan=-999.;
fid = fopen(cfile,'r');
if(fid < 0) 
   'error: cannot open file ',cfile
end
% skip text on first four lines
cline = fgetl(fid);
cline = fgetl(fid);
chead = fgetl(fid);
cline = fgetl(fid);
temp = strsplit(':',cline);
stname = temp{2};

% skip text in coordinate line
cline = fgetl(fid);
stdec = sscanf(cline,'coordinate %f %f declination %f');
cline = fgetl(fid);
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
cline = fgetl(fid);
z = zeros(2,nche*nbt) + i*zeros(2,nche*nbt);
sig_e = zeros(nche,nche*nbt) + i*zeros(nche,nche*nbt);
sig_s = zeros(2,2*nbt) + i*zeros(2,2*nbt);
ndf = zeros(nbt,1);
cform120 = ['period : %f    decimation level %d'...
                '   freq. band from %d to %d'];
cform125 = 'number of data point %d sampling freq. %f Hz'
periodi = zeros(1,nbt);
level = zeros(1,nbt);
ibandlim = zeros(2,nbt);
sampRate = zeros(1,nbt);
for ib = 1:nbt
  cline = fgetl(fid);
  tmp=sscanf(cline,cform120);
  if isempty(tmp)==0,
   periods(ib) = tmp(1);
   level(ib) = tmp(2);
   ibandlim(:,ib) = tmp(3:4);
  else
   periods(ib) = NaN;
  end
  cline = fgetl(fid);
  tmp=sscanf(cline,cform125);
  if isempty(tmp)==0,
   ndf(ib) = tmp(1);
   sampRate(ib) = tmp(2);
  else
   ndf(ib) = NaN;  
  end
  k1 = nche*(ib-1) + 1;
  k2 = nche*ib;
  cline = fgetl(fid);
  ztemp = fscanf(fid,'%e',[4,nche]);
  ztemp(find(ztemp==rnan))=NaN;
  if size(ztemp)>=[4,nche],
    z(1:2,k1:k2) = ztemp(1:2:3,:)+i*ztemp(2:2:4,:);   
  else
   z(1:2,k1:k2) =   NaN; 
  end
  cline = fgetl(fid);
  cline = fgetl(fid);
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
  cline = fgetl(fid);
  cline = fgetl(fid);
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
  cline = fgetl(fid);
end  

for k = 1:nch
   CH{k} = char(chid(k,:));
   STA{k} = char(csta(k,:));
end
size(sig_e)
sig_e = reshape(sig_e,[nche,nche,nbt]);
Z = reshape(z,[2,nche,nbt]);
sig_s = reshape(sig_s,[2,2,nbt]);
Z = struct('TF',Z,'SIG_E',sig_e,'SIG_S',sig_s,...
	'T',periods,'ndf',ndf,'stcor',stdec(1:2),...
        'theta0',stdec(3),'orient',orient,...
        'chead',chead,'stname',stname,...
        'ibandlim',ibandlim,'level',level,'sampRate',sampRate,...
	    'Nch',nch,'Nche',nche,'nbt',nbt,'chid',{CH},'sta',{STA});
