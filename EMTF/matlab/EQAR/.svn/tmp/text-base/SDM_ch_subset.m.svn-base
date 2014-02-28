function [SDMs,SDMHDs,IER] = SDM_ch_subset(SDM,SDMHD,ch,band,grouping)
%  Usage: [SDMs,SDMHDs] = SDM_ch_subset(SDM,SDMHD,ch,band,grouping);
%  Usage: [SDMs,SDMHDs] = SDM_ch_subset(SDM,SDMHD,ch,band);
%  Usage: [SDMs,SDMHDs] = SDM_ch_subset(SDM,SDMHD,ch);
%  given input SDM data structure for a set of nt channels,
%   compute eigenvectors, etc. for a subset of channels given
%   by the integer array ch
%
%   optional argument band is used to select a subset of freq bands
%   optional argument grouping is used to select channel grouping
%     for incoherent noise estimation


IER = 0;
if nargin < 5
   grouping = 'standard';
end
if nargin < 4
   band = [1:SDMHD.nbt];
end
if max(ch) > SDMHD.nt
   fprintf(1,'%s\n','Error in SDM_ch_subset')
   return
end

%  extract header appropriate to the
%   selected channels (NOTE: we are assuming (in our definition
%   of the output ih) that the first channel selected for each
%   site is Hx; if this is not true, use of SDMHDs.ih in
%   some routines called subsequently could result in errors
SDMHDs = SDMHD;
nt = length(ch);
nbt = length(band);
StaNum = zeros(SDMHD.nt,1);
ih = [SDMHD.ih SDMHD.nt+1];
for k = 1:SDMHD.nsta
   StaNum(ih(k):ih(k+1)-1) = k;
end
ChUse = zeros(size(StaNum));
ChUse(ch) = 1;
StaNumUse = StaNum(ch);
nch = [];
ihs = [];
StUse = zeros(SDMHD.nsta,1);
for k = 1:SDMHD.nsta
   nch1 = sum(ChUse(ih(k):ih(k+1)-1));
   if nch1 > 0
      nch = [nch nch1];
      StUse(k) = 1;
      ih1 = min(find(StaNumUse == k));
      ihs = [ihs ih1];
   end
end
SDMHDs.nbt = nbt;
SDMHDs.nt = nt;
SDMHDs.nsta = sum(StUse);
SDMHDs.nch = nch;
SDMHDs.ih = ihs;
ind = find(StUse);
SDMHDs.stcor = SDMHD.stcor(:,ind);
SDMHDs.decl = SDMHD.decl(ind);
SDMHDs.chid = SDMHD.chid(:,ch);
SDMHDs.sta = SDMHD.sta(:,ind);
SDMHDs.orient = SDMHD.orient(:,ch);
SDMHDs.ch_name = SDMHD.ch_name(ch,:);
   
%  now extract appropriate parts of SDMs, compute
%   noise variances with standard approach,
%   and find eigenvectors
SDMtemp = SDM;
SDMtemp.S = SDMtemp.S(ch,ch,band);
SDMtemp.var = SDMtemp.var(ch,band);
SDMtemp.Sig = SDMtemp.Sig(ch,band);
SDMtemp.T = SDMtemp.T(band);
SDMtemp.nf = SDMtemp.nf(band);
SDMs = SDMtemp;
SDMtemp.Hd = SDMHDs;
[SDMs.var,SDMs.Sig,IER] = SDMvar(SDMtemp,grouping);
[SDMs] = ReCompEvec(SDMs);
