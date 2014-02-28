switch TFvariant
   case 'Variance grouping: channels'
      grouping = 'all';
      [Sdms.var,sig] = SDMvar(Sdms,grouping);
      [Sdms] = ReCompEvec(Sdms);
   case 'CC 1'
      nf = Sdms.nf;
      nbt = length(nf)
      %  set number of total modes (from evals of full SDM)
      sig = 3;
      CcovSig = 2;
      CohSig = 2;
      
      %   sigLevel should depend on nf (at least!)
      sigLevel = sig*ones(nbt,1);
      temp = sig./sqrt(nf(nf<100)/100);
      sigLevel(nf<100) = temp;
      dim = 2*ones(nbt,1);
      %  using significant eigenvectors form singular SDM (in physical, not SNR, units)
      S = zeros(size(Sdms.S));
      for ib = 1:nbt
         dim(ib) = max(dim(ib),sum(Sdms.lambda(:,ib)>sigLevel(ib)));
         u = diag(sqrt(Sdms.var(:,ib)))squeeze(Sdms.U(:,1:dim(ib),ib));
         ev = Sdms.lambda(1:dim(ib),ib);
         S(:,:,ib) = u*diag(ev)*u';
      end
      Sdms.S = S;
      % set number of inter-site coherent modes (using canonical 
      %     covariances of input/output channels, projected into
      %     space of dominant modes
      [CC] = compCC(Sdms,Kin,Kout); 
      dimMax = max(dim);
      nIn = length(Kin);
      nOut = length(Kout);
      nT = nIn+nOut
      V = zeros(nT,dimMax,nbt);
      
      sigLevel = CcovSig*ones(nbt,1);
      temp = CcovSig./sqrt(nf(nf<100)/100);
      sigLevel(nf<100) = temp;
      CCdim = 2*ones(nbt,1);
      
      for ib = 1:nbt
         CCdim(ib) = max(CCdim(ib),sum(CC.ccov(:,ib)>sigLevel(ib)));
         n = CCdim(ib);
         m = dim(ib);
         
         %  construct coherent modes from canonical covs (SNR coordinates here)
         V(:,1:n,ib) = [CC.u1(:,1:n,ib);CC.u2(:,1:n,ib)];
         if n < m
             S11 = Sdms.S(Kin,Kin,ib);
             sigInv = 1./sqrt(Sdms.var(Kin,ib));
             S11 = diag(sigInv)*S11*diag(sigInv);
             Z = CC.u1(:,n+1:m,ib)'*S11(CC.u1(n+1:m,ib);
             [u,e] = eig(Z);
             ind = find(e>CohSig);
             nI = length(ind);
             if nI>0
                v = u(:,ind)'*CC.u1(:,n+1:m,ib);
                V(:,n+1:n+nI,ib) = [v;zeros(nOut,nI)];
             end
             DIM(ib) = n+nI;
         else
             DIM(ib) = n; 
         end
         sig = sqrt(Sdms.var([Kin Kout],ib));
         V(:,:,ib) = diag(sig)*V(:,:,ib);
      end
      
