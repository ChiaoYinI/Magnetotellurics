function [ind] = TOD2hr(t,iBin)

%  Usage: ind] = TOD2hr(t,iBin);
%   given time t (in days) return indicies for times in 2 hr bin
%    number iBin (0-2 UT = bin 1, etc.)

tBin = ceil((t-floor(t))*12);
ind = find(tBin == iBin);