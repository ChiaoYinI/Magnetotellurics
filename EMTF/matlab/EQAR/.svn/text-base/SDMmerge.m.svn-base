function [SDMS,SDMHD,TimeInDays,good] = ...
	SDMmerge(dayOne,dayEnd,YEAR,BAND,SAVE)
%
%   Usage: [SDMS,SDMHD,Day,good] = ...
%        SDMmerge(dayOne,dayEnd,YEAR,BAND,SAVE);
%
%   makes a list of SDM files (1 hz, 1 day files, or 40 hz, 2 hr filed)
%   in the range dayOne,dayEnd for YEAR
%   NOTE:  YEAR is now a character string giving the directory
%      name, not a numerical year
%   and then loads these into cell arrays SDMS, SDMHD;
%   also returns array Day giving start time (in decimal days) for each file,a
%   and a 0/1 indicator array for "good" data
%
%   Calls : sdm_init, loadsdms
%

PKDSAOdata = setDataDir('PKDSAOdata');
%data dir may neep to be set explicitly:
%PKDSAOdata='/home/gauss/scratch/PKDSAOg/data/';

DataDir = [PKDSAOdata  YEAR '/SDM/'];

if BAND == 40
   outFileAll = ...
   [DataDir 'S0_' num2str(dayOne) '_' num2str(dayEnd) '_40Hz.mat'];
else
   outFileAll = ...
   [DataDir 'S0_' num2str(dayOne) '_' num2str(dayEnd) '_1Hz.mat'];
end

%   make list of sdm files
    Hours = {'00','02', '04', '06', '08', '10', '12', ...
		 '14', '16', '18', '20', '22'};
    nHours =  length(Hours);

    k = 0;
    for day = dayOne:dayEnd
       if(day<10)
          Day = ['00' num2str(day)];
       elseif (day<100)
          Day = ['0' num2str(day)];
       else
          Day = num2str(day);
       end

       if BAND == 1
          k = k+1;
          SDMfiles{k} = [DataDir Day '.S0'];
          TimeInDays(k) = day;
       else
          for hr = 1:nHours
            k = k+1;
            SDMfiles{k} = [DataDir Day '_' Hours{hr} '.S0'];
            HOUR = str2num(Hours{hr});
            TimeInDays(k) = day+HOUR/24;
          end
       end
    end

%   try to read in all files in list, make cell array of SDM structure
    nFiles = length(SDMfiles);
    NBTall = zeros(nFiles,1);
    NTall = zeros(nFiles,1);
    good = ones(nFiles,1);

    for k = 1:nFiles
       k;
       SDMS{k} = [];
       SDMHD{k} = [];
       % see if the file exists
       if(good(k))
          fid = fopen (SDMfiles{k},'r');
          if fid < 0 
          % 	if not, continue with next file
             good(k) = 0;
             SDMfiles{k};
               ' doesn''t exist';
             % 	if it exists, close file, call sdm_init
          else
             fclose(fid);
	     [fid_sdm,irecl,nbt,nt,nsta,nsig,nch, ...
		   ih,stcor,decl,sta,chid,csta, ...
		   orient,periods] = sdm_init(SDMfiles{k});
             if fid_sdm < 0 
                good(k) = 0;
             else
                col = ':' * ones(nt,1);
                ch_name = [char(chid(1:2,:)'),char(col),char(csta'), ...
		       char(col),num2str(orient(1,:)')];
                SdmHd = struct('nbt',nbt,'nt',nt,'nsta',nsta,'nch',nch, ...
		   'ih',ih,'stcor',stcor,'decl',decl,'chid',chid, ...
		   'sta',sta,'orient',orient,'ch_name',ch_name);
                [Sdms,err] = loadsdms(fid_sdm,irecl,nbt,nt); 
                fclose(fid_sdm);
                if ( err < 0 )
                   good(k) = 0;
                else
                   SDMS{k} = Sdms;
                   SDMHD{k} = SdmHd;
                   NBTall(k) = SdmHd.nbt;
                   NTall(k) = SdmHd.nt;
                end
                if any(SDMS{k}.Sig==0)
                   good(k) = 0;
                end
             end
          end    % if (fid)
       end
    end       	%  for k = 1:nFiles

%   now define "normal" array size
    NBTnormal = median(NBTall(NBTall~=0));
    NTnormal = median(NTall(NTall~=0));
%   mark as bad all that are not normal size
    for k = 1:nFiles
      if(NBTall(k) ~= NBTnormal | NTall(k) ~= NTnormal)
         good(k) = 0;
      end
    end
    if SAVE
       eval(['save ' outFileAll ' SDMS SDMHD TimeInDays good dayOne dayEnd'])
    end
