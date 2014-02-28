%  makes list of starting points in FC files for each frequency
%    at a fixed decimation level (for all stations);  
%   call after fc_open, set id = desired decimation level

function [start_freqs] = mk_start_freqs(id,fids,nf,nch,start_decs)

nhdrec = 20;
nn = size(start_decs);nsta = nn(1);nd = nn(2);
start_freqs = zeros(nf,nsta);
for ista = 1:nsta
  irecl = nch(ista)+1; 
  start_freqs(1,ista) = start_decs(ista,id);
  fseek(fids(ista),start_decs(ista,id),'bof') ;
  for l = 2:nf 
    head = fread(fids(ista),irecl,'long'); 
    if(head(1) > id ) 
       break 
    end 
    nskip = 4*irecl*(head(3) + 1 ); 
    start_freqs(l,ista) = start_freqs(l-1,ista)+nskip + 4*irecl;
    fseek(fids(ista),nskip,'cof'); 
  end 
end 
