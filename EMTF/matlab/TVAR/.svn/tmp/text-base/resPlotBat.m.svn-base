function resPlotBat(day1,day2,year,OPTIONS)
% Usage: resPlotBat(day1,day2,year,OPTIONS)
%  runs resplt script for all days between day1:day2
%  inclusive for given year. 
%  Does not plot results, just saves to a mat file
%  with signal and residual FT structures for each
%  day saved in  2 cell arrays (one each for sig and res)
%
%  OPTIONS data structure:
%     OPTIONS.resplt :  provides resplt script name
%     OPTIONS.saveFile : file name for output mat file
%     OPTIONS.*  any needed input parameters for the
%       given resplt script; field names should match
%       resplt input parameter names in "Usage:"
%       (e.g., site, EH, bsFile, id)  ... DATE is
%       always constructed
%  NOTE:  only coded for some resplt scripts

k = 0;
Pl = 0;
for day = day1:day2
   DATE = struct('day',day,'year',year);
   fprintf(1,'%s ','Day = ');
   fprintf(1,'%d \n',day);
   k = k+1;

   switch OPTIONS.resplt
      case 'resplt_logf'
         [Sig{k},Res{k}] = resplt_logf(DATE,OPTIONS.site,...
		OPTIONS.EH,OPTIONS.bsFile,Pl);
      case 'resplt2hr_Day'
         [Sig{k},Res{k}] = resplt2hr_Day(DATE,OPTIONS.site,...
		OPTIONS.EH,OPTIONS.id,Pl);
   end
end
eval(['save ' OPTIONS.saveFile ' Res Sig OPTIONS day1 day2 year']);
