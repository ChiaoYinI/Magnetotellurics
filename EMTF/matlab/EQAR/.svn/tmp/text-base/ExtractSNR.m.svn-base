function [x] = ExtractSNR(SDMs,good,ib,ich)

x = zeros(size(SDMs));
for k = 1:length(SDMs)
   if(good(k))
      x(k) = SDMs{k}.Sig(ich,ib)/SDMs{k}.var(ich,ib);
   else
      x(k) = NaN;
   end
end
