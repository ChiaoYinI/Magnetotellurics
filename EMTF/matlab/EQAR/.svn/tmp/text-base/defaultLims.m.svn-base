function [lims] = defulatLims(f)

[n,m] = size(f);
f = reshape(f,n*m,1);
f = f(~isnan(f));
fMed = median(f)
fMad = median(abs(f-fMed));
if(fMad > 3)
  fMad = ceil(fMad/5)*5;
end
ll = floor(fMed-2*fMad);
ul = ceil(fMed+2*fMad);
if(ll == ul)
  ul = ll + 1;
end
lims = [ll,ul];
