function sigLevel = SDMsig(nf);
nbt = length(nf);
sig = 3;
      %   sigLevel should depend on nf (at least!)
sigLevel = sig*ones(nbt,1);
temp = sig./sqrt(nf(nf<100)/100);
sigLevel(nf<100) = temp;