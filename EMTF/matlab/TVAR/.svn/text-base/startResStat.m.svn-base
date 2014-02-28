%   script for computing quantiles of residual distribution

%  load residual or signal data structure
[filename, pathname] = uigetfile('*.mat', 'Signal/Residual mat file');
cfile = [pathname filename];
eval(['load ' cfile]);

Pval = [.01:.01:.99];
TimeDiv = struct('nDiv',12,'func','TOD2hr');

[RD] = resStat(RES,Pval,TimeDiv);