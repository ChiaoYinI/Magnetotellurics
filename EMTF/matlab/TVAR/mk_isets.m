%   mk_isets : for decimation level id, makes a list of set #s which
%   are available for all stations.  For all available sets  
%   the relative location in the FC files is returned in sets_pt 
%   USAGE: [isets,isets_pt] =  mk_isets(fids,start_decs,nch,nf,id);
%
function [isets,isets_pt] =  mk_isets(fids,start_decs,nch,nf,id);
nsetsmx = 50000;
nm = size(start_decs); nd = nm(2); nsta = nm(1);

isets_pt = zeros(nsetsmx,nsta);
for ista = 1:nsta
   isets = get_isets(fids(ista),start_decs(ista,id),nch(ista));
   nsets(ista) = length(isets);
   isets_pt(1:nsets(ista),ista) = isets';
end
first_set = min(isets_pt(1,:));
for ista = 1:nsta
  isets = isets_pt(1:nsets(ista),ista) - first_set + 1;
  isets_pt(:,ista) = zeros(nsetsmx,1);
  isets_pt(isets,ista) = [1:nsets(ista)]';
end
%   find all indices where all stations have data
inds = find( ~ any( isets_pt' == 0 ) );
isets_pt = isets_pt(inds,:);
isets =  isets(isets_pt(:,nsta)) + first_set - 1;
