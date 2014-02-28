function [lims] = pwstlims(T,amp,ph,tr,ti,plt_type)
%  given plot: 
%      REIM = Real and Imaginary parts of TF, linear scale
%      AMPH = Amplitude and phase of TF, linear scale
%      LAMP = Amplitude and phase of TF, log scale

qMaxAmp = .8;
qMaxPh = .9;
qMinPh = 1.-qMaxPh;
ampMin = .5;
lims = zeros(6);
Tmin = floor(2*log10(min(T)))/2;
Tmin = 10^Tmin;
Tmax = ceil(2*log10(max(T)))/2;
Tmax = 10^Tmax;
lims(1:2) = [Tmin Tmax];
if(plt_type == 'REIM')
   % real and imaginary parts
   nm = size(tr); nm = nm(1)*nm(2);
   temp = [ reshape(abs(tr),nm,1); reshape(abs(ti),nm,1)];
   temp = sort(temp);
   trmax = 1.2*temp(round(qMaxAmp*2*nm)); 
   trexp = floor(log10(trmax));
   trexp = 10^trexp;
   trmax = trmax/trexp;
   if(trmax > 5)
      amplim = 10*trexp;
   elseif(trmax > 2)
      amplim = 5*trexp;
   elseif(trmax > 1.5)
      amplim = 2*trexp;
   elseif(trmax > 1)
      amplim = 1.5*trexp;
   else
      amplim = trexp;
   end
   amplim = max(amplim,ampMin);
   lims(3:6) = [ -amplim amplim -amplim amplim];
elseif(plt_type == 'AMPH')
   % Linear Amplitude and Phase
   nm = size(amp); nm = nm(1)*nm(2);
   temp = reshape(amp,nm,1);
   temp = sort(temp);
   trmax = temp(round(qMaxAmp*nm)); 
   trexp = floor(log10(trmax));
   trexp = 10^trexp;
   trmax = trmax/trexp;
   if(trmax > 5)
      amplim = 10*trexp;
   elseif(trmax > 2)
      amplim = 5*trexp;
   elseif(trmax > 1.5)
      amplim = 2*trexp;
   elseif(trmax > 1)
      amplim = 1.5*trexp;
   else
      amplim = trexp;
   end
   amplim = max(amplim,ampMin);
   lims(3:4) = [  0 amplim ];
else
   % Log amplitude and phase
   nm = size(amp); nm = nm(1)*nm(2);
   temp = reshape(amp,nm,1);

   ampmin = floor(2*log10(min(temp)))/2;
   ampmin = 10^ampmin;
   ampmax = ceil(2*log10(max(temp)))/2;
   ampmax = 10^ampmax;
   lims(3:4) = [ampmin ampmax];
end

if(plt_type == 'AMPH' | plt_type == 'LAMP')
   phs_lims = [10. 22.5 45 60. 90 180];
   nm = size(amp); nm = nm(1)*nm(2);
   temp = reshape(ph,nm,1);
   phmax = max(temp); phmin = min(temp);
   temp = sort(temp);
   phmax = temp(round(qMaxPh*nm)); 
   phmin = temp(round(qMinPh*nm)); 
   if(phmax > 0 & phmin < 0 )
      phmax = max([phmax -phmin]);
      for k = 1:length(phs_lims)
         if(phmax < phs_lims(k))
            lims(5:6) = [ -phs_lims(k) phs_lims(k) ];
            break
         end
      end
   elseif(phmin > 0 )
      for k = 1:length(phs_lims)
         if(phmax < phs_lims(k))
            lims(5:6) = [ 0 phs_lims(k) ];
            break
         end
      end
  elseif(phmax < 0 )
      for k = 1:length(phs_lims)
         if(-phmin < phs_lims(k))
            lims(5:6) = [ -phs_lims(k) 0 ];
            break
         end
      end
   end
end
