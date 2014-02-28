function [TFres] = TFsdmSet(cfile,Kout,Kin,TFTYPE);
% Usage: [TFres] = TFsdmSet(cfile);
%        [TFres] = TFsdmSet(cfile,Kout);
%        [TFres] = TFsdmSet(cfile,Kout,Kin);
%        [TFres] = TFsdmSet(cfile,Kout,Kin,TFTYPE);
%   reads a standard sdm (*.S0) file
%   and sets up residual TF data structure 
%     NOTE: eigenvalues of SDM added to TF data structure 
%
%   Kout, Kin are output and input channel lists
%    these are optional arguments; if not present
%    simplest default forms (all channels in/out)
%    is assumed.
%   TFTYPE  provides control over exactly how the TF is derived
%    from the SDM file.  The idea is that different "TFtypes"
%    can be devised, and coded here; each will in general require
%    some control parameters which will also be specified as
%    fields in the data structure TFTYPE
%    Diffferent TFtypes may produce different outputs, within limits
%    IN ALL CASES: the output should have the TFres structure, although
%     additional fields could be added to TFres.  The main fields in
%     TFres that are necessary are complete channel identifiers for
%     (at least) all input and output channels, and a set of vectors
%     in TFres.V for each frequency band, one for each signal/noise mode
%     to fit for that band.  These vectors are to be in actual
%     (as opposed to SNR) units.  In the simplest application 
%     (TFtype = SDM) the vectors are just the two dominant 
%     eigenvectors from the SDM ... a more complex approac, using
%     the SDM and canonical coherence analysis is also implemented.
%
%   TFsdmSet returns empty TFres (not a structure!) 
%   if there are any errors in reading files

TFfileName = '';
if nargin < 4 
   TFTYPE.TFtype = 'SDM';
   TFTYPE.nEvec = 2;
end

if cfile(end-3:end) == '.mat'
   eval(['load ' cfile])
   err =  exist(Sdms) & exist(SdmHd) -1;
else
   [Sdms,SdmHd,err] = sdmstruc(cfile);
end
if err < 0
   TFres = [];
   return
end
Sdms.Hd = SdmHd;

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
[nt,nbt] = size(Sdms.var);

if TFTYPE.TFtype == 'SDM'
   %   simplest approach ... evecs of SDM
   %   In this case evecs for all channels are output in the original
   %    order
   nv = TFTYPE.nEvec;
   
   V = zeros(nt,nv,nbt);
   for ib = 1:SdmHd.nbt
      V(:,:,ib) = diag(sqrt(Sdms.var(:,ib)))*squeeze(Sdms.U(:,1:nv,ib));
   end
   for k = 1:nbt
      UseVout{k} = [1:nv];
   end

   TFres = struct('type',2,'PhysU',1,...
	'FileType',TFTYPE.TFtype,'FileName',cfile,'Name',TFfileName,...
	'dim',nv,'nf',nbt,'nch',nt,...
	'ch_id',char(SdmHd.chid'),'sta',sta, ...
	'Kin',Kin,'Kout',Kout,...
	'T',Sdms.T,'V',V,'var',Sdms.var,...
	'Resid',1,'nTseg',1,'UseVout',{UseVout});
elseif TFTYPE.TFtype == 'CC1'
   %   experimental approach: canonical coherence analysis
   %   to estimate number of modes that are coherent between sites
   %      (actually between input and output channels) +
   %        coherent noise at predicting site 
   %   TFTYPE contains also these fields:
   %      TFTYPE.UseAllCh = 1 to use (potentially) all channels 
   %          (Kin+Kout)for prediction of coherent part of data 
   %          (otherwise use only channels declared in Kin)
   %      TFTYPE.CohUseIn(3) = 1 to use intersite, input, output
   %          coherent components in prediction input.  To use output,
   %          UseAllCh must be set to 1.  For intersite the first
   %          two (signal) modes are always used, even for CohUseIn(1)=0
   %          Setting CohUseIn(k) = 0 elimantes any columns of V from
   %          coherence group k 
   %      TFTYPE.CohUseOut(3) = 1  as for CohUseIn, but for output.
   %          Note that CohUseIn(k) = 1 is required for CohUseOut(k)=1
   %          Setting CohUseOut(k) = 0 when CohUseIn(k) = 1 means
   %          that the input components of V will be used in Vin, but 
   %          the corresponding components will be zeroed in Vout.
   %          To keep track of this a flag will be set, and checked in
   %          InterpTF_FT, where the actual prediction TF is assembled.
   %
   %   In this case only channels Kin Kout are output, reordered
   %    as in the Kin and Kout lists (Kin first)
   %    initial input/output groups must be disjoint, but if UseAllCh
   %      is 1, Kin will be reset to Kall.

   %    first estimate coherent noise vectors
   [V,DIM,CH,B11,B22,Sdms] = CCsigNoiseVec(Sdms,Kin,Kout,TFTYPE.SigLevFunc);

   if isempty(V)
      fprintf(1,'%s \n','Error in CCsigNoiseVec')
      TFres = [];
      return
   else
      [nt,ndim,NBT] = size(V);
      %  reset/sort sta, ch_id, var to correspond to order 
      %    of channels in V
      Kall = [Kin Kout];
      nt = length(Kall);
      sta = sta(Kall,:);
      ch_id = char(SdmHd.chid(:,Kall)');
      var = Sdms.var;

      nIn = length(Kin);
      nOut = length(Kout);
      nAll = length(Kall);
      Kout = [nIn+1:nAll];
      if TFTYPE.UseAllCh == 1
         Kin = [1:nAll];
      else
         Kin = [1:nIn];
         TFTYPE.CohUseIn(3) = 0;
      end
      for k = 1:3
         if TFTYPE.CohUseIn(k) == 0
            TFTYPE.CohUseOut(k) = 0;
         end
      end
      for ib = 1:NBT
         useVin = zeros(ndim,1);
         useVout = zeros(ndim,1);
         useVin(1:2) = 1;
         useVout(1:2) = 1;
         i2 = 2;
         if DIM(1,ib)>2
            i1 = i2+1; i2 = DIM(1,ib);
            useVin(i1:i2) = TFTYPE.CohUseIn(1);
            useVout(i1:i2) = TFTYPE.CohUseOut(1);
         end
         if DIM(2,ib)>0
            i1 = i2+1; i2 = i2+DIM(2,ib);
            useVin(i1:i2) = TFTYPE.CohUseIn(2);
            useVout(i1:i2) = TFTYPE.CohUseOut(2);
         end
         if DIM(3,ib)>0
            i1 = i2+1; i2 = i2+DIM(3,ib);
            useVin(i1:i2) = TFTYPE.CohUseIn(3);
         end
         indUse = find(useVin);
         nUse = length(indUse);
         dim(ib) = nUse;
         V(:,1:nUse,ib) = V(:,indUse,ib);
         UseVout{ib} = useVout(indUse);
         if nUse<ndim
            V(:,nUse+1:ndim,ib) = 0;
         end
      end

      TFres = struct('type',3,'PhysU',1,...
	'FileType',TFTYPE.TFtype,'FileName',cfile,'Name',TFfileName,...
	'dim',dim,'DIM',DIM,'nf',nbt,'nch',nt,...
	'ch_id',ch_id,'sta',sta, ...
	'Kin',Kin,'Kout',Kout,...
	'T',Sdms.T,'V',V,'var',var,'Resid',1,...
	'nTseg',1,'UseVout',{UseVout});
      end
   end
end
