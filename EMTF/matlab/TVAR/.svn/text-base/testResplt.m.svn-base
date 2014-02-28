path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/TVAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/EQAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/SDM',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/IN',path);
dir = '/home/server/pi/homes/egbert/PKDSAOg/data/2004/';
bsFile = ...
  '/home/server/pi/homes/egbert/PKDSAOg/data/work/CF/bs_40hzC.cfg'
DATE = struct('year',2004,'day',272,'hour',16);
id = 1;
site = 'PKD';
EH = 'H';
resplt2hr(DATE,site,EH,id)
%resplt2hr_Day(DATE,site,EH,id)
%resplt_logf(DATE,site,EH,bsFile)

print -f3 -djpeg90 PKDmag272log.jpg
print -f4 -djpeg90 PKDmagRes272log.jpg
%delete(1); delete(2);

site = 'PKD';
EH = 'E'
%resplt2hr_Day(DATE,site,EH,id)
%resplt2hr(DATE,site,EH,id)
resplt_logf(DATE,site,EH,bsFile)
print -f1 -djpeg90 PKDelec272log.jpg
print -f2 -djpeg90 PKDelecRes272log.jpg
%delete(1); delete(2);

%site = 'SAO'
%EH = 'H'
%resplt2hr_Day(DATE,site,EH,id)
%resplt2hr(DATE,site,EH,id)
%print -f1 -depsc SAOmag270.eps
%print -f2 -depsc SAOmagRes270.eps
%delete(1); delete(2);

%site = 'SAO'
%EH = 'E'
%resplt2hr_Day(DATE,site,EH,id)
%print -f1 -depsc SAOelec270.eps
%print -f2 -depsc SAOelecRes270.eps
