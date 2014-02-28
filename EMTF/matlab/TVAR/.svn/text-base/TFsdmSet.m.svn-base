function [TFres] = TFsdmSet(cfile,Kout,Kin,DIM);
% Usage: [TFres] = TFsdmSet(cfile);
%        [TFres] = TFsdmSet(cfile,Kout);
%        [TFres] = TFsdmSet(cfile,Kout,Kin);
%        [TFres] = TFsdmSet(cfile,Kout,Kin,DIM);
%   reads a standard sdm (*.S0) file
%   and sets up residual TF data structure 
%     NOTE: eigenvalues of SDM added to TF data structure 
%
%   Kout, Kin are output and input channel lists
%    these are optional arguments; if not present
%    simplest default forms (all channels in/out)
%    is assumed.
%   DIM  ... yet to be developed data structure which
%    will control how the dimension of the TF is determined
%    to start with this is just a fixed number ...
%
%   TFsdmSet returns empty TFres (not a structure!) 
%   if there are any errors in reading 

TFtype = 'SDM';
TFfileName = '';
if nargin < 4 
   DIM.nEvec = 2;
end

[Sdms,SdmHd,err] = sdmstruc(cfile);
if err < 0
   TFres = [];
   return
end

if nargin < 3
   Kin = [1:SdmHd.nt];
end
if nargin < 2
   Kout = [1:SdmHd.nt];
end

sta = [];
for ista = 1:SdmHd.nsta
   for k = 1:SdmHd.nch(ista)
      sta = [sta ; char(SdmHd.sta(:,ista)')];
   end
end
V = zeros(size(Sdms.U));
for k = 1:SdmHd.nbt
   V(:,:,k) = diag(sqrt(Sdms.var(:,k)))*squeeze(Sdms.U(:,:,k));
end

TFres = struct('type',2,'PhysU',1,...
	'FileType',TFtype,'FileName',cfile,'Name',TFfileName,...
	'dim',DIM.nEvec,'nf',SdmHd.nbt,'nch',SdmHd.nt,...
	'ch_id',char(SdmHd.chid'),'sta',sta, ...
	'Kin',Kin,'Kout',Kout,...
	'T',Sdms.T,'V',V,'var',Sdms.var,'Resid',1,'nTseg',1,...
        'lambda',Sdms.lambda,...
	'CLEAN',0,'ClnX',0);
