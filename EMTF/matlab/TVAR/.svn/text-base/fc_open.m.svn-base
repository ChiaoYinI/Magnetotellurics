function [fid,nd,nf,nch,chid,orient,drs,stdec,decs,start_dec] = fc_open(cfile);

%    opens one FC file, reads header, saves some useful info about file
%   
%    THIS VERSION ASSUMES FC FILES ARE PACKED !!!!
%
% Usage: [fid,nd,nf,nch,chid,orient,drs,stdec,decs,start_dec] = fc_open(cfile);

cform1 = 'nch %d nd %d nfmx %d';
nhdrec = 20;

fid = fopen(cfile,'r');
if(fid < 0 )
   fprintf(1,'file not found/n')
   return
end
line = fgets(fid);
if(line == -1)
   fid = -1
   return
end
nn = sscanf(line,cform1);
nch = nn(1);nd = nn(2);nf = nn(3);
i0 = 28 + 1;
i1 = 28+8*nd;
decs = sscanf(line(i0:i1),'%d',[2,nd]);
i0 = i1+5;
i1 = i0 + nd*12 - 1;
drs = sscanf(line(i0:i1),'%e',nd);
i1 = i1+8;
chid = []; orient = [];
for ich = 1:nch
  i0 = i1 + 1 ;
  i1 = i0+3; 
  chid = [chid;line(i0:i1)];
  i0 = i1 + 1;
  i1 = i0 + 15;
  orien = sscanf(line(i0:i1),'%f',2);
  orient = [orient , orien ];
end
i0 = i1 + 5; i1 = i0+30;
stdec = sscanf(line(i0:i1),'%f %f decl %f '); 

irecl = nch+1;
start_dec(1) = nhdrec*4*irecl;
fseek(fid,start_dec(1),'bof');
i0 = 0;
for k=2:nd
   start_dec(k) = start_dec(k-1) + i0;
   i0 = 0;
   for l = 1:nf
      head = fread(fid,irecl,'long');
      if(head(1) > k-1 )
         nskip = 4*irecl*(head(3) + 1 );
         fseek(fid,nskip,'cof');
         i0 = nskip + 4*irecl;
         break
      end
      nskip = 4*irecl*(head(3) + 1 );
      fseek(fid,nskip,'cof');
      start_dec(k) = start_dec(k) + nskip + 4*irecl;
   end
end
