%path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/TVAR',path);
%path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/EQAR',path);
%path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/SDM',path);
%path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/IN',path);

path('/home/gauss/scratch/PKDSAOg/matlab/EMTF/TVAR',path);
path('/home/gauss/scratch/PKDSAOg/matlab/EMTF/EQAR',path);
path('/home/gauss/scratch/PKDSAOg/matlab/EMTF/SDM',path);
path('/home/gauss/scratch/PKDSAOg/matlab/EMTF/IN',path);
path('/home/gauss/scratch/PKDSAOg/matlab/EMTF/TVAR/OLD',path);
path('/home/gauss/scratch/PKDSAOg/matlab',path);

bsFile = ...
  '/home/server/pi/homes/egbert/PKDSAOg/data/work/CF/bs_test.cfg'
DATE = struct('year',2004,'day',140);
id = 1;
site = 'PKD';
EH = 'H';
%resPlt1hz(DATE,site,EH,id)
Ndays = 20;
[SigAll,ResAll] = resPlt1hz_logf(DATE,site,EH,bsFile,Ndays)
%resPlt1hz10Day(DATE,site,EH,1)
