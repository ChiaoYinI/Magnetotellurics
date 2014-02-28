k1 = min(find(good));
if BAND == 1
   dl = 0;
else
   dl = 1;
end

switch plotType
   case 'SDM Eigenvalues'
     %   plot 4 largest eigenvalues
     nEvec = 4;
     %   color axis limits: hard wired at present
     cl1 = [0,40;0,40;0,40;0,40];
     cl40 = [0,35;0,35;-5,15;-5,5];

     % construct data array
     nFiles = length(SDMS);
     NBT = SDMHD{k1}.nbt;
     NT = SDMHD{k1}.nt;
     D = zeros(nFiles,NBT,nEvec);
     for k = 1:nFiles
        if good(k)
           D(k,:,:) = 10*log10(abs(SDMS{k}.lambda(1:nEvec,:))');
        else
           D(k,:,:) = D(k,:,:)/0; 
        end
     end

     %   set OPTIONS array: titles, axis labels, limits
     for k = 1:nEvec
       OPTIONS{k} = struct('Caxis',cl1(k,:),...
	  'SubTitle',['# ' num2str(k)],...
          'Title','',...
	  'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	  'ColorAxisLabel','SNR (dB)',...
	  'DayLine',dl);
         if BAND == 40
            OPTIONS{k}.Caxis = cl40(k,:);
         end
      end
      OPTIONS{1}.Title = 'SDM Eigenvalues 1-4';
      OPTIONS{1}.LabelColor = [0 0 0];
   case 'Magnetic Field SNR'
     %   plot Horizontal Magnetic Field SNR
     iComp = find(char(SDMHD{k1}.chid(1,:)=='H'));
     nComp = length(iComp);

     %   color axis limits: hard wired at present
     cl = ones(nComp,1)*[0,40];

     % construct data array
     nFiles = length(SDMS);
     NBT = SDMHD{k1}.nbt;
     D = zeros(nFiles,NBT,nComp);
     for k = 1:nFiles
        if good(k)
           D(k,:,:) = 10*log10( ...
		abs(SDMS{k}.Sig(iComp,:)'./SDMS{k}.var(iComp,:)'));
        else
           D(k,:,:) = D(k,:,:)/0; 
        end
     end

     %   set OPTIONS array: titles, axis labels, limits
     for k = 1:nComp
       OPTIONS{k} = struct('Caxis',cl(k,:),...
	  'SubTitle',SDMHD{k1}.ch_name(iComp(k),:),...
          'Title','',...
	  'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	  'ColorAxisLabel','SNR (dB)',...
	  'DayLine',dl);
      end
      OPTIONS{1}.Title = 'Magnetic Field SNR';
      OPTIONS{1}.LabelColor = [0 0 0];
   case 'Magnetic Field Signal Variation'
     %   plot Horizontal Magnetic Field Signal Power Variations
     %    relative to mean power
     iComp = find(char(SDMHD{k1}.chid(1,:)=='H'));
     nComp = length(iComp);

     %   color axis limits: hard wired at present
     cl = ones(nComp,1)*[-10,10];

     % construct data array
     nFiles = length(SDMS);
     NBT = SDMHD{k1}.nbt;
     D = zeros(nFiles,NBT,nComp);
     for k = 1:nFiles
        if good(k)
           D(k,:,:) = 10*log10( ...
		SDMS{k}.Sig(iComp,:)'./SdmAvg.Sig(iComp,:)');
        else
           D(k,:,:) = D(k,:,:)/0; 
        end
     end

     %   set OPTIONS array: titles, axis labels, limits
     for k = 1:nComp
       OPTIONS{k} = struct('Caxis',cl(k,:),...
	  'SubTitle',SDMHD{k1}.ch_name(iComp(k),:),...
          'Title','',...
	  'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	  'ColorAxisLabel','dB',...
	  'DayLine',dl);
      end
      OPTIONS{1}.Title = 'Magnetic Field Signal Power Variations';
      OPTIONS{1}.LabelColor = [0 0 0];
   case 'Electric Field SNR'
     %   plot Electric Field SNR
     %   Just like magnetic field,  but different components, title

     iComp = find(char(SDMHD{k1}.chid(1,:)~='H'));
     nComp = length(iComp);

     %   color axis limits: hard wired at present
     cl = ones(nComp,1)*[0,35];

     % construct data array
     nFiles = length(SDMS);
     NBT = SDMHD{k1}.nbt;
     D = zeros(nFiles,NBT,nComp);
     for k = 1:nFiles
        if good(k)
           D(k,:,:) = 10*log10( ...
		abs(SDMS{k}.Sig(iComp,:)'./SDMS{k}.var(iComp,:)'));
        else
           D(k,:,:) = D(k,:,:)/0; 
        end
     end

     %   set OPTIONS array: titles, axis labels, limits
     for k = 1:nComp
       OPTIONS{k} = struct('Caxis',cl(k,:),...
	    'SubTitle',SDMHD{k1}.ch_name(iComp(k),:),...
          'Title','',...
	    'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	    'ColorAxisLabel','SNR (dB)',...
	    'DayLine',dl);
     end
     OPTIONS{1}.Title = 'Electric Field SNR';
      OPTIONS{1}.LabelColor = [0 0 0];
      
   case 'Canonical Coherences: PKD vs. SAO'
      ind1 = [1:SDMHD{k1}.ih(2)-1];
      ind2 = [SDMHD{k1}.ih(2):SDMHD{k1}.nt];
      nCC = 4;
      
     %   color axis limits: hard wired at present
     cl = ones(nCC,1)*[0,1];
       
      % construct data array
      nFiles = length(SDMS);
      NBT = SDMHD{k1}.nbt;
      D = zeros(nFiles,NBT,nCC);
      for k = 1:nFiles
         if good(k)
            [CC] = compCC(SDMS{k},ind1,ind2);
            D(k,:,:) = CC.ccor(1:nCC,:)';
         else
            D(k,:,:) = D(k,:,:)/0; 
         end
      end
      %   set OPTIONS array: titles, axis labels, limits
      for k = 1:nCC
        OPTIONS{k} = struct('Caxis',cl(k,:),...
	     'SubTitle',['CC #' num2str(k)],...
          'Title','',...
	     'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	     'ColorAxisLabel','R^2',...
	     'DayLine',dl);
      end
      OPTIONS{1}.Title = 'Canonical Coherences: PKD-SAO';
      OPTIONS{1}.LabelColor = [0 0 0];
       
   case 'Canonical Coherences: H vs. E'
      ind1 = find(char(SDMHD{k1}.chid(1,:)=='H'))
      ind2 = find(char(SDMHD{k1}.chid(1,:)=='E'))
      nCC = 4;
      
     %   color axis limits: hard wired at present
     cl = ones(nCC,1)*[0,1];
       
      % construct data array
      nFiles = length(SDMS);
      NBT = SDMHD{k1}.nbt;
      D = zeros(nFiles,NBT,nCC);
      for k = 1:nFiles
         if good(k)
            [CC] = compCC(SDMS{k},ind1,ind2);
            D(k,:,:) = CC.ccor(1:nCC,:)';
         else
            D(k,:,:) = D(k,:,:)/0; 
         end
      end
      %   set OPTIONS array: titles, axis labels, limits
      for k = 1:nCC
        OPTIONS{k} = struct('Caxis',cl(k,:),...
	     'SubTitle',['CC #' num2str(k)],...
          'Title','',...
	     'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	     'ColorAxisLabel','R^2',...
	     'DayLine',dl);
      end
      OPTIONS{1}.Title = 'Canonical Coherences: H-E';
      OPTIONS{1}.LabelColor = [0 0 0];
   % remaining cases may not work any more !!!
   case 'Average Modes'
     nCH_miss = 3;
     fracLowSNR = 1/3;
     minSNR = 2;
     %  Compute average SDM over "good" files:
     %   use standard grouping for "local noise" definition
     if ~SdmAvgComputed
        grouping = 'standard';
        [SdmAvg] = SDMavg(SDMS,SDMHD,good,grouping);
        SdmAvgComputed = 1;
        Sdms = SdmAvg;
        SdmHd = SdmAvg.Hd;
        cfile = ['SDMAVG_' num2str(dayOne) '-' num2str(dayEnd) '.mat'];
        eval(['save ' cfile ' Sdms SdmHd']);
     end
     SNRavg = SdmAvg.Sig./SdmAvg.var;

     %   plot power in 4 dominant modes (of average SDM)
     %    as a function of data segment
     nEvec = 4;
     %   color axis limits: hard wired at present
     cl40 = [0,35;0,35;-5,15;-5,15];
     cl1 = [0,40;0,40;0,40;0,40];

     % construct data array
     nFiles = length(SDMS);
     NBT = SDMHD{k1}.nbt;
     NT = SDMHD{k1}.nt;
     D = zeros(nFiles,NBT,4);
     ih = SdmAvg.Hd.ih;
     ind1 = [1:ih(2)-1];
     ind2 = [ih(2):NT];
     for k = 1:nFiles
        if good(k)
           SNR = (SDMS{k}.Sig./SDMS{k}.var);
           lowSNR = (SNR < fac*SNRavg | SNR < minSNR);
           ChUse = find(sum(lowSNR,2) < NBT*fracLowSNR);
           PKDind = find(ChUse<ih(2));
           SAOind = find(ChUse>=ih(2));
           nUse = length(ChUse);
           nUsePKD = length(PKDind);
           nUseSAO = length(SAOind);
           %  could also omit days with two many missing channels
           %   at one site 
           if((nUse >= NT-nCH_miss)&(nUsePKD>=2)&(nUseSAO>=2))
              for ib = 1:NBT
                 % u1 is in avg SNR units ... first transform to physical
                 u1 = squeeze(SdmAvg.U(:,1:nEvec,ib));
           %      modified to use incoherent noise scaling from
           %       averaged SDM (so scaling stays invariant with time)
           %      sigAvg = diag(sqrt(SdmAvg.var(:,ib)));
           %      sigAvgInv = diag(sqrt(1./SdmAvg.var(:,ib)));
           %      u1 = sigAvg*u1;
                 sigInv = diag(1./sqrt(SdmAvg.var(ChUse,ib)));
                 v1 = u1(ChUse,:);
                 %   also transform segment SDM to segment SNR units
                 S = sigInv*SDMS{k}.S(ChUse,ChUse,ib)*sigInv;
                 %  now we essentially fit data (in SNR units) as 
                 %  linear combination of columns of u1, and compute variance
                 %  of regression coeficients, averaged over segment ...
                 w = v1/(v1'*v1);
                 s = real(diag(w'*S*w));
                 s(1) = s(1)+s(2);
                 s(2) = s(3)+s(4);
                 %%  residuals
                 H = eye(nUse) - w*v1';
                 %   residuals at PKD, SAO
                 R = real(diag(H*S*H));
                 s(3) =  sum(R(PKDind))/max(nUsePKD-2,1);
                 s(4) =  sum(R(SAOind))/max(nUseSAO-2,1);
                 D(k,ib,:) = 10*log10(s);
              end
           else
              D(k,:,:) = D(k,:,:)/0; 
           end
        else
           D(k,:,:) = D(k,:,:)/0; 
        end
     end

     %   set OPTIONS array: titles, axis labels, limits
     for k = 1:nEvec
       OPTIONS{k} = struct('Caxis',cl1(k,:),...
	  'SubTitle',['# ' num2str(k)],...
          'Title','',...
	  'TimeAxisLabel','Days', ...
          'PeriodAxisLabel','log_{10} Period (s)', ...
	  'ColorAxisLabel','SNR (dB)',...
	  'DayLine',dl);
         if BAND == 40
            OPTIONS{k}.Caxis = cl40(k,:);
         end
      end
      OPTIONS{1}.SubTitle = 'Modes 1+2';
      OPTIONS{1}.Caxis = [0,45];
      OPTIONS{2}.SubTitle = 'Modes 3+4';
      OPTIONS{2}.Caxis = [0,25];
      OPTIONS{3}.SubTitle = 'PKD residual';
      OPTIONS{3}.Caxis = [-10,10];
      OPTIONS{4}.SubTitle = 'SAO residual';
      OPTIONS{4}.Caxis = [-10, 10];
      OPTIONS{1}.Title = 'Power in Modes of Averaged SDM';
      OPTIONS{1}.LabelColor = [0 0 0];
   case 'Coherent Noise Setup'
     %   CALL this before specific Coherent Noise or Residual
     %           computations ... this is effected by TFTYPE
     %  This case is experimental ...  not described in any
     %    publication  ... maybe not a good idea!

     % next we use canonical coherence analysis to find 
     %   (up to) three sets of vectors:
     %   coherent between sites
     %   coherent between channels at first, second sites 
     %      (need to be clear here about channel groupings for
     %         definition of incoherent noise levels!)
     k1 = min(find(good));
     [nt,NBT] = size(SDMS{k1}.var);
     ih = SdmAvg.Hd.ih;
     Kin = TFTYPE.Kin;
     Kout = TFTYPE.Kout;
     SigLevFunc{1} = 'CcovSig';
     SigLevFunc{2} = 'SDMSig';
     SigLevFunc{3} = 'SDMSig';

     if binByHours
        Signal = zeros(nFiles,NBT,3);
        CohNoise = zeros(nFiles,NBT,2);
        Res = zeros(nFiles,NBT,2);
        ResCH = zeros(nFiles,NBT,nt);
   
        hour = round((TimeInDays-floor(TimeInDays))*24);
        DIM = zeros(3,NBT);
        for hr = 0:2:22
           ind = find(hour == hr);
           %  First compute average SDM over "good" files:
           %   use standard grouping for "local noise" definition
           grouping = 'standard';
           for k = 1:length(ind)
              SDMS1{k} = SDMS{ind(k)};
              SDMHD1{k} = SDMHD{ind(k)};
           end
           [SdmAvg] = SDMavg(SDMS1,SDMHD1,good(ind),grouping);
           [V,dim,CH,B11,B22,SdmAvg] = ...
		   CCsigNoiseVec(SdmAvg,Kin,Kout,SigLevFunc);
           %  Fit the vectors to data, and compute average power in
           %   each signal group for each processed segment 
           [Signal1,CohNoise1,Res1,ResCH1] =  ...
		CC_sigNoise(SDMS1,good(ind),dim,V,TFTYPE);
           DIM = max(DIM,dim);
           Signal(ind,:,:) = Signal1;
           CohNoise(ind,:,:) = CohNoise1;
           Res(ind,:,:) = Res1;
           ResCH(ind,:,:) = ResCH1;
        end
     else
     %   average over all time
        if ~SdmAvgComputed
           %  First compute average SDM over "good" files:
           %   use standard grouping for "local noise" definition
           grouping = 'standard';
           [SdmAvg] = SDMavg(SDMS,SDMHD,good,grouping);
           SdmAvgComputed = 1;
        end
     
        [V,DIM,CH,B11,B22,SdmAvg] = ...
		CCsigNoiseVec(SdmAvg,Kin,Kout,SigLevFunc);

        %  Fit the vectors to data, and compute average power in
        %   each signal group for each processed segment 
        [Signal,CohNoise,Res] = ...
		 CC_sigNoise2(SDMS,good,DIM,V,TFTYPE);
%        [Signal,CohNoise,Res,ResCH] = ...
%		 CC_sigNoise(SDMS,good,DIM,V,TFTYPE);
     end

  case 'Coherent Noise'
     %   set up plotting for coherent noise
     CohNoise(CohNoise <=1 ) = 1;
     Res(Res <=1 ) = 1;
     CohNoise = 10*log10(CohNoise);
     Res = 10*log10(Res);
     plotLabel{1} = 'Coherent PKD';
     plotLabel{2} = 'Coherent SAO';
     plotLabel{3} = 'Residual PKD';
     plotLabel{4} = 'Residual SAO';
     cl1 = [0,15;0,15;0,15;0,15];
     cl40 = [0,15;0,15;0,15;0,15];
     for k = 1:4
        OPTIONS{k} = struct('Caxis',cl1(k,:),...
       'SubTitle',plotLabel{k},...
        'Title','',...
       'TimeAxisLabel','Days', ...
       'PeriodAxisLabel','log_{10} Period (s)', ...
       'ColorAxisLabel','SNR (dB)',...
       'DayLine',dl);
       if BAND == 40
          OPTIONS{k}.Caxis = cl40(k,:);
       end
     end
     if binByHours
        OPTIONS{1}.Title = ['Coherent and Incoherent Noise '...
	': Binned by Hour : ' TFTYPE.title ];
     else
        OPTIONS{1}.Title = ['Coherent and Incoherent Noise : '...
		 TFTYPE.title ];
     end
     OPTIONS{1}.LabelColor = [.8 .8 .8];
     D = zeros(nFiles,NBT,4);
     D(:,:,1:2) = CohNoise;
     D(:,:,3:4) = Res;
  case 'Coherent Between Sites'
     %   set up plotting for signals
     Signal(Signal <=1 ) = 1;
     % convert to dB
     Signal = 10*log10(Signal);

     cohNoise = max(DIM(1,:))>2;
     plotLabel{1} = 'Signal #1';
     plotLabel{2} = 'Signal #2';
     if cohNoise
        plotLabel{3} = 'Coh between PKD/SAO';
     end
     cl1 = [0,40;0,40;0,40];
     cl40 = [0,35;0,35;0,25];
     for k = 1:2+cohNoise
        OPTIONS{k} = struct('Caxis',cl1(k,:),...
        'SubTitle',plotLabel{k},...
        'Title','',...
        'TimeAxisLabel','Days', ...
        'PeriodAxisLabel','log_{10} Period (s)', ...
        'ColorAxisLabel','SNR (dB)',...
        'DayLine',dl);
        if BAND == 40
           OPTIONS{k}.Caxis = cl40(k,:);
        end
     end
     D = zeros(nFiles,NBT,2+cohNoise);
     D(:,:,1:2+cohNoise) = Signal;
     if binByHours
        OPTIONS{1}.Title = ['Signal and Coherent Noise '...
	': Binned by Hour : ' TFTYPE.title ];
     else
        OPTIONS{1}.Title = ['Signal and Coherent Noise : '...
	 TFTYPE.title ];
     end
     OPTIONS{1}.LabelColor = [.8 .8 .8];

  case 'Residuals: PKD'
     chPKD = [1:7];
     % convert to dB
     D = 10*log10(ResCH(:,:,chPKD));
     cl = [-5,10];
     for k = 1:length(chPKD)
        OPTIONS{k} = struct('Caxis',cl,...
        'SubTitle',char(SDMHD{k1}.chid(:,chPKD(k)))',...
        'Title','',...
        'TimeAxisLabel','Days', ...
        'PeriodAxisLabel','log_{10} Period (s)', ...
        'ColorAxisLabel','SNR (dB)',...
        'DayLine',dl);
     end
     if binByHours
        OPTIONS{1}.Title = [ 'Residuals : PKD : TFTYPE',...
		TFTYPE.title ' : Binned by Hour'];
     else
        OPTIONS{1}.Title = [ 'Residuals : PKD : TFTYPE',...
		TFTYPE.title ];
     end

  case 'Residuals: SAO'
     chSAO = [8:12];
     % convert to dB
     D = 10*log10(ResCH(:,:,chSAO));
     cl = [-5,10];
     for k = 1:length(chSAO)
        OPTIONS{k} = struct('Caxis',cl,...
        'SubTitle',char(SDMHD{k1}.chid(:,chSAO(k)))',...
        'Title','',...
        'TimeAxisLabel','Days', ...
        'PeriodAxisLabel','log_{10} Period (s)', ...
        'ColorAxisLabel','SNR (dB)',...
        'DayLine',dl);
     end
     if binByHours
        OPTIONS{1}.Title = [ 'Residuals : SAO : TFTYPE',...
		TFTYPE.title ' : Binned by Hour'];
     else
        OPTIONS{1}.Title = [ 'Residuals : SAO : TFTYPE',...
		TFTYPE.title ];
     end
end
