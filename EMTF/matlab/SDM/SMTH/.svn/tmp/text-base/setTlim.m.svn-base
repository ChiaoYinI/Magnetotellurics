function [Tmin,Tmax] = setTlim(periods)
%   set nicer period limits for logartihmic period scale plots

  x_min = min(periods);
  x_max = max(periods);
  Tmin = 10^(floor(log10(x_min)*2)/2);
  if ((log10(x_min)-log10(Tmin)) < 0.15)
          Tmin = 10^(log10(Tmin)-0.3);
  end
  Tmax = 10^(ceil(log10(x_max)*2)/2);
  if ((log10(Tmax)-log10(x_max)) < 0.15)
          Tmax = 10^(log10(Tmax)+0.3);
  end
