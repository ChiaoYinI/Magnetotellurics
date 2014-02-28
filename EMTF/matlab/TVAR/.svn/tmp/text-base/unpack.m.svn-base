%  unpacks integers read in from FC files into complex

function [fc] = unpack(ifc,scales)

nn = size(ifc); nsets = nn(2);
M = 65536; Mo2 = M/2;
it2 = floor(ifc/M);
it1 = ifc - M*it2 - Mo2;
k1 = it1 - 8*floor(it1/8);
k2 = it2 - 8*floor(it2/8);
it1 = (it1-k1)/Mo2 + .000122;
it2 = (it2-k2)/Mo2 + .000122;
fc = it1.*(10.^k1) + i* it2.*(10.^k2);
fc = fc.*(scales* ones(1,nsets));
