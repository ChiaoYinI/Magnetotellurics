function [TFres] = TFsdmSet(cfile,Kout,Kin,TYPE);
% Usage: [TFres] = TFsdmSet(cfile);
%        [TFres] = TFsdmSet(cfile,Kout);
%        [TFres] = TFsdmSet(cfile,Kout,Kin);
%        [TFres] = TFsdmSet(cfile,Kout,Kin,TYPE);
%   reads a standard sdm (*.S0) file
%   and sets up residual TF data structure 
%     NOTE: eigenvalues of SDM added to TF data structure 
%
%   Kout, Kin are output and input channel lists
%    these are optional arguments; if not present
%    simplest default forms (all channels in/out)
%    is assumed.
%   TYPE  provides control over exactly how the TF is derived
%    from the SDM file.  The idea is that different "TFtypes"
%    can be devised, and coded here; each will in general require
%    some control parameters which will also be specified as
%    fields in the data structure TYPE
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
   TYPE.TFtype 'SDM';
   TYPE.nEvec = 2;
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
[nt,nbt] = size(Sdms.var);

if TYPE.TFtype == 'SDM'
   %   simplest approach ... evecs of SDM
   %   In this case evecs for all channels are output in the original
   %    order
   nv = TYPE.nevec;
   DIM = [nv*ones(nbt,1),zeros(nbt,2)];
   
   V = zeros(nt,nv,nbt));
   for ib = 1:SdmHd.nbt
      V(:,:,ib) = diag(sqrt(Sdms.var(:,ib)))*squeeze(Sdms.U(:,1:nv,ib));
   end

   TFres = struct('type',2,'PhysU',1,...
	'FileType',TYPE.TFtype,'FileName',cfile,'Name',TFfileName,...
	'dim',DIM,'nf',nbt,'nch',nt,...
	'ch_id',char(SdmHd.chid'),'sta',sta, ...
	'Kin',Kin,'Kout',Kout,...
	'T',Sdms.T,'V',V,'var',Sdms.var,'Resid',1,'nTseg',1)
elseif TYPE.TFtype == 'CC1'
   %   experimental approach: combine evec, canonical coherence analysis
   %   to : estimate number of modes that are coherent between sites
   %      (actually between input and output channels) +
   %        coherent noise at predicting site 
   %
   %   In this case only channels Kin Kout are output, reordered
   %    as in the Kin and Kout lists (Kin first)

   %   First reorder channels and eliminate any channels not in
   %     Kin/Kout lists
   %
   %   groups are assumed disjoint for now ... check this first:
   K1 = Kin'*ones(size(Kout));
   K2 = ones(size(Kin'))*Kout;
   if any(any(K1==K2))
      fprintf(1,'%s \n',...
	'For TFype = CC1 Kin and Kout must be disjoint lists')
      return
   end

   %  now reorder
   Kall = [Kin Kout];
   Sdms.var = Sdms.var(Kall,:);
   Sdms.S = Sdms.S(Kall,Kall,:);
   Sdms.U = Sdms.S(Kall,Kall,:);
   Sdms.Sig = Sdms.Sig(Kall,:);
   SdmHd.nt = length(Kall);
   SdmHd.chid = SdmHd.chid(:,Kall);  
   Sdms.Hd = SdmHd;
   sta = sta(Kall);   
   nIn = length(Kin);
   Kin = [1:nIn]; 
   nOut = length(Kout);
   Kout = [nIn+1:nIn+nOut];
   nT = nIn+nOut
   DIM = zeros(nbt,3);
   %   NOTE: the reordering of channels is not completely general
   %   or correct.  SdmHd.ih is not correct ... ih is used in
   %   SDMvar if grouping = 'standard' ... but this is the default
   %   grouping so for this case omit the next lines which ...
   %   Recompute noise variance using single channel grouping:
   grouping = 'all';
   [Sdms.var,sig] = SDMvar(Sdms,grouping);
   [Sdms] = ReCompEvec(Sdms);
   nf = Sdms.nf;
   nbt = length(nf)
   %  set cuttoff to detrmine number of total modes from evals of full SDM
   %   functions which determine cutoffs as a function of sample size
   %    are passed as character strings in TYPE
   eval(['sigLevelSDM = ' TYPE.SDMsig '(nf);']);
   eval(['sigLevelCcov = ' TYPE.CcovSig '(nf);']);
   eval(['sigLevelCcoh = ' TYPE.CcohSig '(nf);']);

   %  dimSDM : not used?
   for ib = 1:nbt
      dimSDM(ib) = max(dimSDM(ib),sum(Sdms.lambda(:,ib)>sigLevelSDM(ib)));
   end
   % set number of inter-site coherent modes using canonical
   %     covariances of input/output channels  ... not SDM
   [CC] = compCC(Sdms,Kin,Kout);
   dimMax = max(nIn,nOut);
   V = zeros(nT,dimMax,nbt);
   %  dimCC will be number of modes coherent between sites used for
   %   intersite prediction
   dimCC = 2*ones(nbt,1);
   for ib = 1:nbt
      dimCC(ib) = max(dimCC(ib),sum(CC.ccov(:,ib)>sigLevelCcov(ib)));
      n = dimCC(ib);
      %  construct coherent modes from canonical covs (SNR coordinates here)
      U1 = CC.u1(:,1:n,ib);
      U2 = CC.u2(:,1:n,ib);
      V(:,1:n,ib) = [U1;U2];
      sigInv = 1./sqrt(Sdms.var(Kin,ib));
      S11 = diag(sigInv)*Sdms.S(Kin,Kin,ib)*diag(sigInv);
      Q = eye(nIn)-U1/(U1'*U1)*U1';
      [u,e] = eig(Q*S11*Q);
      e = diag(real(e));
      ind = find(e>sigLevelCcov(ib));
      nI = length(ind);
      if nI>0
         U1p = u(:,ind); 
         C = U1'*S11*U1p/(U1p'*S11*U1p-eye(nI));
         V(:,n+1:n+nI,ib) = [U1p+U1*C;zeros(nOut,nI)];
         DIM(ib,2) = nI;
      end
      DIM(ib,1) = n;
      sig = sqrt(Sdms.var(:,ib));
      V(:,:,ib) = diag(sig)*V(:,:,ib);
   end
   TFres = struct('type',2,'PhysU',1,...
	'FileType',TYPE.TFtype,'FileName',cfile,'Name',TFfileName,...
	'dim',DIM,'nf',nbt,'nch',nt,...
	'ch_id',char(SdmHd.chid'),'sta',sta, ...
	'Kin',Kin,'Kout',Kout,...
	'T',Sdms.T,'V',V,'var',Sdms.var,'Resid',1,'nTseg',1)
end
