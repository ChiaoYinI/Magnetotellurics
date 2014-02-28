%   script for plotting signal and residual amplitudes

%  load residual or signal data structure
[filename, pathname] = uigetfile('*.mat', 'Signal/Residual mat file');
cfile = [pathname filename];
eval(['load ' cfile]);

OPTIONS.plotAmp = 1;
OPTIONS.plotRes = 0;
OPTIONS.kSeg = 1;
nSeg = length(Sig);
for k = 1:nSeg
    segments{k} = num2str(Sig{k}.t(1),5);
end
OPTIONS.segments = segments;
ResSig = 'Sig';

plotAmp(ResSig,Sig{1},Res{1},OPTIONS);
