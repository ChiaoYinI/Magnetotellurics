cbk = get(gcbo,'Tag');

switch cbk
   case 'GETIB'
     % finds band which is closest to frequncy where slider is released,
     %  then displays in small box above plot button

     per = get(gcbo,'Value');
     per = 10^per;
     ib = find_ib(nbt,periods,per);
     h_ib = findobj('Tag','IB Box','Parent',gcf); 
     if(~isempty(h_ib))
        delete(h_ib);
        clear h_ib;
     end
     h_ib = axes('Position',[.05,.34,.05,.04],...
        'XTickLabelMode','manual', ...
        'YtickLabelMode','manual',...
        'box','on',...
        'Xlim',[0,1],...
        'Ylim',[0,1],'Tag','IB Box'); 
     text('Parent',h_ib,'Position',[.18,.32],'String',num2str(ib),...
        'FontWeight','bold');
  case 'PLOT'
     UD = get(gcf,'UserData');
     per = get(findobj('Tag','GETIB','Parent',gcf),'Value');
     per = 10^per
     ib = find_ib(nbt,periods,per);
     ivec = UD.EVEC.ivec;
     rho_ref = UD.EVEC.rho_ref;
     if UD.EVEC.l_PWG
        S = squeeze(Sdms.S(:,:,ib));
        var = Sdms.var(:,ib);
        PWGopt = struct('nEvec',UD.EVEC.PWGnum,...
		'evalWt',UD.EVEC.PWGevalWt,...
		'theta',0,'TeTm',0)
        Vpwg = pwg(S,Sdms.Hd,var,PWGopt)
        %  for now l_impNorm is ignored for PWG case
        evecPlt(ib,ivec,Sdms,SdmHd,UD.EVEC,Vpwg);
     else
        if UD.EVEC.l_impNorm
        %  normalize all E components by impedances
           U = impNorm(Sdms,UD.EVEC.Zall,ivec,ib,rho_ref);
           evecPlt(ib,ivec,Sdms,SdmHd,UD.EVEC,U);
        else
           evecPlt(ib,ivec,Sdms,SdmHd,UD.EVEC);
        end
     end
  case 'PlotHz'
     EV = get(gcf,'UserData');
     EV.l_Hz = get(gcbo,'Value');
     set(gcf,'UserData',EV); 
  case 'PWG'
     EV = get(gcf,'UserData');
     EV.l_PWG = get(gcbo,'Value');
     set(gcf,'UserData',EV); 
  case 'PWGnum'
     EV = get(gcf,'UserData');
     EV.PWGnum = str2num(get(gcbo,'string'));
     set(gcf,'UserData',EV); 
  case 'PWGevalWt'
     EV = get(gcf,'UserData');
     EV.PWGevalWt = str2num(get(gcbo,'string'));
     set(gcf,'UserData',EV); 
  case 'Impedance Scaling'
     EV = get(gcf,'UserData');
     EV.l_impNorm = get(gcbo,'Value');
     if EV.l_impNorm 
        %  set up impedances, add to EV
        EV.Zall = setImp(Sdms);
     end
     set(gcf,'UserData',EV); 
  case 'eigSM'
     eigSM; 
     if(exist('hfig_evec_menu')) 
        chk_clr(hfig_evec_menu);
        clear hfig_evec_menu;
     end
  case 'PLTSM'
     l_smthsep=1-l_smthsep; set(gcbo,'Value',l_smthsep);
     if(l_smthsep)
        Uplt = zeros(size(Uplt));
        Uplt(:,1:2,:) = Usm;
        for ib = 1:nb
           if(Neig(ib) > 0 )
              Uplt(:,3:2+Neig(ib),:) = Vsm(:,1:Neig(ib),:);
           end
        end
     else
        Uplt = Sdms.U;
        for ib = 1:nbt
          N = sqrt(Sdms.var(:,ib));
          Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
        end
     end
  case 'PERP'
     l_perp=l_smthsep*(1-l_perp); set(gcbo,'Value',l_perp);
     if(l_perp)
        for ib = 1:nb
           if(Neig(ib) > 0 )
              Uplt(:,3:2+Neig(ib),:) = Uperp(:,1:Neig(ib),:);
           end
        end
     else
        for ib = 1:nb
           if(Neig(ib) > 0 )
              Uplt(:,3:2+Neig(ib),:) = Vsm(:,1:Neig(ib),:);
           end
        end
     end        
     
  case 'MTTFs'
     l_MTTF=l_smthsep*(1-l_MTTF); set(gcbo,'Value',l_MTTF);
     if(l_MTTF)
        Uplt(:,1:2,:) = TFsm;
     else
        Uplt(:,1:2,:) = Usm;
     end        
  case 'Cancel'
    close(gcf);
    l_perp = l_perp_old;
    l_MTTF = l_MTTF_old;
    if(l_smthsep)
       if(l_MTTF)
         Uplt(:,1:2,:) = TFsm;
       else
         Uplt(:,1:2,:) = Usm;
       end        
       if(l_perp)
         for ib = 1:nb
           if(Neig(ib) > 0 )
              Uplt(:,3:2+Neig(ib),:) = Uperp(:,1:Neig(ib),:);
           end
         end
       else
         for ib = 1:nb
           if(Neig(ib) > 0 )
              Uplt(:,3:2+Neig(ib),:) = Vsm(:,1:Neig(ib),:);
           end
         end
       end
     end 
  case '# Evec'
     EV = get(gcf,'UserData');
     EV.n_evec = str2num(get(gcbo,'string'));
     set(gcf,'UserData',EV);
  case 'rho ref'
     EV = get(gcf,'UserData');
     EV.rho_ref=str2num(get(gcbo,'string'));
     set(gcf,'UserData',EV);
  case 'SNR units'
     EV = get(gcf,'UserData');
     EV.snr_units=get(gcbo,'Value');
     set(gcf,'UserData',EV);
  case 'Ellipse/Vectors'
     EV = get(gcf,'UserData')
     EV.l_ellipse=get(gcbo,'Value');
     set(gcf,'UserData',EV);
  case 'Apply'
%    save EVEC_OPTIONS from popup menu to main figure UserData
     EV = get(gcf,'UserData')
     EV.ivec = [1:EV.n_evec];
     UD = get(EV.parentFig,'UserData')
     UD.EVEC = EV;
     set(EV.parentFig,'UserData',UD);
     close(gcf);
     clear EV;
  case 'Clr Evec Figs'
    delete(findobj('Tag','evec'));
  otherwise
     fprintf(1,'%s \n','Case Not Coded in evecCbk')
end
