path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/TVAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/EQAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/SDM',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/IN',path);
dir = '/home/server/pi/homes/egbert/PKDSAOg/data/2004/';
FileNamesFC = {[dir 'FC/PKD_270_00.f7'],[dir 'FC/SAO_270_00.f5']};
FileNamesSysTF = {[dir 'sysTF/PKD_270_00.stf'],[dir 'sysTF/SAO_270_00.stf']};
FileNameSDM = [dir 'SDM/270_00.S0'];
id = 1;
%  limit frequencies to those appropriate to decimation level?
ifreq = [1:108];
csta = {'PKD','SAO'};

id = 1;
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
FTamp = FTavg(FT,1);
%temp = abs(squeeze(FT.data(1,:,:)));
%amp = sqrt(mean(temp.*temp,1))
%freq = 40*ifreq/256;
loglog(FTamp.f,FTamp.avgAmp,'b')
f1 = FTamp.f;
hold on

ifreq = [5:112];
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
FTamp = FTavg(FT,1);
loglog(FTamp.f,FTamp.avgAmp,'b+')
loglog(f1,FTamp.avgAmp,'r+')

id = 2;
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
temp = abs(squeeze(FT.data(1,:,:)));
amp = sqrt(mean(temp.*temp,1))
freq = freq/4;
hold on;
loglog(freq,amp,'r');

FTamp = FTavg(FT,1);
hold on
loglog(FTamp.f,FTamp.avgAmp,'r*')

id = 3;
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
temp = abs(squeeze(FT.data(1,:,:)));
amp = sqrt(mean(temp.*temp,1))
freq = freq/4;
hold on;
loglog(freq,amp,'g');

FTamp = FTavg(FT,1);
hold on
loglog(FTamp.f,FTamp.avgAmp,'gd')

id = 4;
[FT] = FTsetup(FileNamesFC,id,ifreq,csta);
temp = abs(squeeze(FT.data(1,:,:)));
amp = sqrt(mean(temp.*temp,1))
freq = freq/4;
hold on;
loglog(freq,amp,'k');
FTamp = FTavg(FT,1);
hold on
loglog(FTamp.f,FTamp.avgAmp,'ko')




