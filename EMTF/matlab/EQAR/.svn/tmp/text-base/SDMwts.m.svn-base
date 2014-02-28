function [wts] = SDMwts(SDMS,SDMHD,good)

n = length(SDMS);
k1 = min(find(good));
nGood = sum(good);
maxFrac = 5/nGood;
tol = .1;
nbt = length(SDMS{k1}.T);
Hpwr = zeros(n,nbt);
ih = SDMHD{k1}.ih;
inds = [];
for k = 1:length(ih)
   inds = [inds ih(k) ih(k)+1];
end

for k =  1:n
   if(good(k))
      for ib = 1:nbt
         Hpwr(k,ib) = sum(SDMS{k}.Sig(inds,ib));
      end
   end
end

wts = ones(size(Hpwr));
for ib = 1:nbt
   notDone = 1;
   while notDone
      temp = Hpwr(:,ib).*wts(:,ib);
      frac = temp./sum(temp(find(good)));
      inds = find(frac > maxFrac);
      wts(inds,ib) = wts(inds,ib).*(maxFrac./frac(inds));
      notDone = max(frac) > maxFrac*(1+tol);
   end
end
   
