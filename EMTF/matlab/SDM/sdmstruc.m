function [Sdms,SdmHd,err] = sdmstruc(s0_file)
%  USAGE:  [Sdms,SdmHd,err] = sdmstruc(s0_file)
%   err = -1 if there is a read error

global SYST

err = 0;
if (isempty(SYST))
   [fid_sdm,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,decl,sta,chid,csta, ...
                                orient,periods] = sdm_init(s0_file);
else
   [fid_sdm,irecl,nbt,nt,nsta,nsig,nch,ih,stcor,decl,sta,chid,csta, ...
                                orient,periods] = sdm_init(s0_file,SYST);
end
if(fid_sdm < 0)
   err = -1;
   return
end
                               
[Sdms,err] = loadsdms(fid_sdm,irecl,nbt,nt);
fclose(fid_sdm);
if err < 0 
   return
end

ch_name = chid(1:2,:);

SdmHd = struct('nbt',nbt,'nt',nt,'nsta',nsta,'nch',nch,'ih',ih,...
   'stcor',stcor,'decl',decl,'chid',chid,'sta',sta,'orient',orient,...
   'ch_name',ch_name);
