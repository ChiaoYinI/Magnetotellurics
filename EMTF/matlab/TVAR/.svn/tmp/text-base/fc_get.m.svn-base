function [fc] = fc_get(fids,start_freqs,nch,isets_pt,id,ifreq)
%  gets FCs for one decimation level, all frequencies in list
%   ifreq, one set of files, all stations
%   NEW Version: Dec 27, 2004 .... returns fc array with 3 dimensions

nn = size(isets_pt); nsets = nn(1); nsta = nn(2); ncht = sum(nch);
nFreqT = length(ifreq);

fc = zeros(ncht,nsets,nFreqT);
ch2 = 0;
for ista = 1:nsta
   ch1 = ch2+1; ch2 = ch2 + nch(ista);
   irecl = nch(ista)+1;
   for ib = 1:nFreqT
      ib1 = ifreq(ib);
      fseek(fids(ista),start_freqs(ib1,ista),'bof');
      head = fread(fids(ista),irecl,'long');
      nsets = head(3);
      head = fread(fids(ista),irecl,'float');
      scales = head(1:nch(ista));
      scales = (scales/1000.);
      ifc = fread(fids(ista),[irecl,nsets],'long');
      fc(ch1:ch2,:,ib) =  ...
            unpack(ifc(2:nch(ista)+1,isets_pt(:,ista)),scales);
   end
end
