path('/home/server/homes/pi/egbert/PKDSAOg/matlab',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/TVAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/EQAR',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/SDM',path);
path('/home/server/homes/pi/egbert/PKDSAOg/matlab/EMTF/IN',path);
dir = '/home/server/pi/homes/egbert/PKDSAOg/data/2004/';
bsFile = ...
  '/home/server/pi/homes/egbert/PKDSAOg/data/work/CF/bs_40hzC.cfg';
day1 = 245;
day2 = 285;
year = 2004;
site = 'PKD';
fclose('all')
EH = 'H';
saveFile = 'H_PKD_245-285_40hz.mat'
OPTIONS = struct('resplt','resplt_logf','EH',EH,...
	'site',site,'bsFile',bsFile,'saveFile',saveFile);
resPlotBat(day1,day2,year,OPTIONS)
fclose('all')

site = 'SAO';
EH = 'H';
saveFile = 'H_SAO_245-285_40hz.mat'
OPTIONS = struct('resplt','resplt_logf','EH',EH,...
	'site',site,'bsFile',bsFile,'saveFile',saveFile);
resPlotBat(day1,day2,year,OPTIONS)
fclose('all')

site = 'PKD';
EH = 'E';
saveFile = 'H_PKD_245-285_40hz.mat'
OPTIONS = struct('resplt','resplt_logf','EH',EH,...
	'site',site,'bsFile',bsFile,'saveFile',saveFile);
resPlotBat(day1,day2,year,OPTIONS)
fclose('all')

site = 'SAO';
EH = 'E';
saveFile = 'E_SAO_245-285_40hz.mat'
OPTIONS = struct('resplt','resplt_logf','EH',EH,...
	'site',site,'bsFile',bsFile,'saveFile',saveFile);
resPlotBat(day1,day2,year,OPTIONS)
fclose('all')
