action = get(gcbo,'Tag');
UD = get(gcf,'UserData')
Sdms = UD.Sdms;

switch action
   case 'Signal Power'
      sig = [Sdms.T';Sdms.Sig];
      hfig_sig = plt_sig(sig,0,chid,stn,sta,title_sig);
   case 'Noise Power'
      noise = [Sdms.T';Sdms.var];
      hfig_noise = plt_sig(noise,0,chid,stn,sta,title_noise);
   case 'SNR'
      sig = [Sdms.T';Sdms.Sig];
      noise = [Sdms.T';Sdms.var];
      hfig_snr = plt_snr(sig,noise,1,chid,stn,sta,arrayid);
   case 'Correlation'
      corrmenu(chid,csta)
   case 'Canonical Correlation'
      [h] = CCmenu(chid,csta,gcf);
   case 'Print'
      sdmPlotFile = [UD.arrayid, '.eps'];
      eval(['print -depsc ' sdmPlotFile ]);
   case 'Quit'
      clearSDMfig
   case 'Eigenvectors'
       evecmenu

   case 'Group: Sites'
      grouping = 'sites';
      [Sdms.var,sig] = SDMvar(Sdms,grouping);
      [Sdms] = ReCompEvec(Sdms);
      set(gcf,'UserData',Sdms);
      Uplt = Sdms.U;
      for ib = 1:nbt
         N = sqrt(Sdms.var(:,ib));
         Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
      end
      E = [Sdms.T';Sdms.lambda];
      clearSDMfig
      sdmPlotSet

   case 'Group: Channels'
      grouping = 'all';
      [Sdms.var,sig] = SDMvar(Sdms,grouping);
      size(Sdms.var)
      [Sdms] = ReCompEvec(Sdms);
      set(gcf,'UserData',Sdms);
      Uplt = Sdms.U;
      for ib = 1:nbt
         N = sqrt(Sdms.var(:,ib));
         Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
      end
      E = [Sdms.T';Sdms.lambda];
      clearSDMfig
      sdmPlotSet

   case 'Group: Sites+E/H'
      grouping = 'standard';
      [Sdms.var,sig] = SDMvar(Sdms,grouping);
      [Sdms] = ReCompEvec(Sdms);
      set(gcf,'UserData',Sdms);
      Uplt = Sdms.U;
      for ib = 1:nbt
         N = sqrt(Sdms.var(:,ib));
         Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
      end
      E = [Sdms.T';Sdms.lambda];
      clearSDMfig
      sdmPlotSet

   case 'Group: Custom'
      grouping = 'custom';
%     need to open window, etc.
%      [Sdms.var,sig] = SDMvar(Sdms,grouping);
%      [Sdms] = ReCompEvec(Sdms);
%      set(gcf,'UserData',Sdms);
%      Uplt = Sdms.U;
%      for ib = 1:nbt
%         N = sqrt(Sdms.var(:,ib));
%         Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
%      end
%      E = [Sdms.T';Sdms.lambda];
   end
      
