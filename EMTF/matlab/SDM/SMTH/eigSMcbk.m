cbk = get(gcbo,'String');

plot_ch = ones(nt,1);
plot_ch(1:2) = 0;
switch cbk
   case 'Recalc'
     set(gcf,'pointer','watch')
     delete(hax);
     clear hax;
     nu(1) = str2num(get(h_DAMP(1),'string')) 
     nu(2) = str2num(get(h_DAMP(2),'string')) 
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end
     if weights_changed
        eigSMset
        weights_changed = 0;
     end
     [alpha,res,misfit,datsz] = eigSMfit(nu,nu_comp);
     eigSMchk;
     eigSMplt;
     temp = squeeze(Smth_dot(1:neig_plt,polarization,:));
     temp = cumsum(temp);
     axes(h_fit)
     semilogx(10.^t,temp(2:neig_plt,:));
     set(gca,'Ylim',fit_lim,'FontWeight','bold','Xlim',Tlim);
     fatlines(gca,2);
     title('Eigenvector Fit')
     legend('1-2','1-3','1-4',4)
     set(gcf,'pointer','arrow')

   case 'Replot'
     set(gcf,'pointer','watch')
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end 

     delete(hax);
     clear hax;
     eigSMplt;
     temp = squeeze(Smth_dot(1:neig_plt,polarization,:));
     temp = cumsum(temp);
     axes(h_fit);
     semilogx(10.^t,temp(2:neig_plt,:));
     set(gca,'Ylim',fit_lim,'FontWeight','bold','Xlim',Tlim);
     fatlines(gca,2);
     title('Eigenvector Fit')
     legend('1-2','1-3','1-4',4)
     set(gcf,'pointer','arrow')

   case 'X'
     set(gcf,'pointer','watch')
%    set polarization at reference site to Hx
     val = get(h_X,'value');
     on_off = {'on' 'off'};
     set(h_Y,'value',1-val); 
     polarization = 2-val;
     set(h_DAMP(1),'Enable',char(on_off(2-val)));
     set(h_DAMP(2),'Enable',char(on_off(val+1)));
     for k = 3:nt
        set(h_damp(1,k),'Enable',char(on_off(2-val)));
        set(h_damp(2,k),'Enable',char(on_off(val+1)));
     end
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end 

     delete(hax);
     clear hax;
     eigSMplt;
     temp = squeeze(Smth_dot(1:neig_plt,polarization,:));
     temp = cumsum(temp);
     axes(h_fit);
     semilogx(10.^t,temp(2:neig_plt,:));
     set(gca,'Ylim',fit_lim,'FontWeight','bold','Xlim',Tlim);
     fatlines(gca,2);
     title('Eigenvector Fit');
     legend('1-2','1-3','1-4',4);
     set(gcf,'pointer','arrow')

   case 'Y'
     set(gcf,'pointer','watch')
     val = get(h_Y,'value');
     on_off = {'on' 'off'};
     set(h_X,'value',1-val); 
     polarization = val+1;
     set(h_DAMP(2),'Enable',char(on_off(2-val)));
     set(h_DAMP(1),'Enable',char(on_off(val+1)));
     for k = 3:nt
        set(h_damp(2,k),'Enable',char(on_off(2-val)));
        set(h_damp(1,k),'Enable',char(on_off(val+1)));
     end
     temp = squeeze(Smth_dot(1:neig_plt,polarization,:));
     temp = cumsum(temp);
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end 

     delete(hax);
     clear hax;
     eigSMplt;
     temp = squeeze(Smth_dot(1:neig_plt,polarization,:));
     temp = cumsum(temp);
     axes(h_fit);
     semilogx(10.^t,temp(2:neig_plt,:));
     set(gca,'Ylim',fit_lim,'FontWeight','bold','Xlim',Tlim);
     fatlines(gca,2);
     title('Eigenvector Fit');
     legend('1-2','1-3','1-4',0);
     set(gcf,'pointer','arrow')

   case 'Re'
     set(gcf,'pointer','watch')
     ri = {'real' 'imag' };
     val = get(h_REAL,'value');
     set(h_IMAG,'value',1-val); 
     reim = 2-val;
     ReIm = char(ri(reim)) ;
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end 
     delete(hax);
     clear hax;
     eigSMplt;
     set(gcf,'pointer','arrow')

   case 'Im'
     set(gcf,'pointer','watch')
     ri = {'real' 'imag' };
     val = get(h_IMAG,'value');
     set(h_REAL,'value',1-val); 
     reim = val+1;
     ReIm = char(ri(reim)) ;
     for k = 3:nt
        ctag = [ 'C' num2str(k)];
        plot_ch(k) = get(findobj('Tag',ctag),'value');
        ctag = [ 'D' num2str(k)];
        nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
        ctag = [ 'd' num2str(k)];
        nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
     end 
     delete(hax);
     clear hax;
     eigSMplt;
     set(gcf,'pointer','arrow')

   case 'ChngFac:'
     axes(h_flambda);
     chngWTS_on = 1;
     weights_changed = 0;
     while chngWTS_on
        [ii,jj,button] = index(Xcell,Ycell) ;
        if(ii > 0 & ii <= nb & jj > 0 & jj <= neig_plt)
           if(button(1:3) == 'nor')
              flambda_scales(jj,ii) = fac*flambda_scales(jj,ii);        
              weights_changed = 1;
           elseif(button(1:3) == undobutton)
              flambda_scales(jj,ii) = flambda_scales(jj,ii)/fac;
              weights_changed = 1;
           else
              chngWTS_on = 0;
           end
        else
           chngWTS_on = 0;
        end
        if chngWTS_on & weights_changed
          %   replot
          neig_plt = min(nt-1,neig_plt);
          temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
          temp = log10(temp);
          temp = [ temp temp(:,nb) ];
          handles = get(h_flambda,'Children');
          nh = length(handles);
          for kk = 1:nh
             htemp = get(handles(kk),'type');
             if htemp(1:4) =='surf'
                set(handles(kk),'Cdata',temp);
                break
             end
          end
          delete(findobj('Tag','TMW_COLORBAR'));
          colorbar;
        end
     end
   case 'FBlk'
     axes(h_flambda);
     boxInd
     if button(1:3) == 'nor'
        fac1 = fac;
     else
        fac1 = 1./fac;
     end
     if(j1 >= 1 & j2 <= neig_plt & i1 >= 1 & i2 <= nbt)
       flambda_scales(j1:j2,i1:i2) = fac1*flambda_scales(j1:j2,i1:i2);        
       weights_changed = 1;
    else
       weights_changed = 0;
    end
    if weights_changed
      %   replot
      neig_plt = min(nt-1,neig_plt)
      temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
      temp = log10(temp);
      temp = [ temp temp(:,nb) ];
      handles = get(h_flambda,'Children');
      nh = length(handles);
      for kk = 1:nh
        htemp = get(handles(kk),'type');
        if htemp(1:4) =='surf'
          set(handles(kk),'Cdata',temp);
          break
        end
      end 
      delete(findobj('Tag','TMW_COLORBAR'));
      colorbar;
    end
  case 'SBlk'
    axes(h_flambda);
    boxInd
    if(button(1:3) == 'nor' & ...
          j1 >= 1 & j2 <= neig_plt & i1 >= 1 & i2 <= nbt )
      flambda_scales(j1:j2,i1:i2) = flambda_set*(1./flambda(j1:j2,i1:i2));
      weights_changed = 1;
    else
      weights_changed = 0;
    end
    if weights_changed
      %   replot
      neig_plt = min(nt-1,neig_plt);
      temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
      temp = log10(temp);
      temp = [ temp temp(:,nb) ];
      handles = get(h_flambda,'Children');
      nh = length(handles);
      for kk = 1:nh
        htemp = get(handles(kk),'type');
        if htemp(1:4) =='surf'
          set(handles(kk),'Cdata',temp);
          break
        end
      end
      delete(findobj('Tag','TMW_COLORBAR'));
      colorbar;
    end
 
   case 'SetWts:'
     axes(h_flambda);
     chngWTS_on = 1;
     while chngWTS_on
        [ii,jj,button] = index(Xcell,Ycell) 
        if(ii > 0 & ii <= nb & jj > 0 & jj <= neig_plt)
           if(button(1:3) == 'nor')
              flambda_scales(jj,ii) = (flambda_set/flambda(jj,ii));        
              weights_changed = 1;
           else
              chngWTS_on = 0;
           end
        else
           chngWTS_on = 0;
        end
        if chngWTS_on & weights_changed
          %   replot
          neig_plt = min(nt-1,neig_plt);
          temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
          temp = log10(temp);
          temp = [ temp temp(:,nb) ];
          handles = get(h_flambda,'Children');
          nh = length(handles);
          for kk = 1:nh
             htemp = get(handles(kk),'type');
             if htemp(1:4) =='surf'
                set(handles(kk),'Cdata',temp);
                break
             end
          end
          delete(findobj('Tag','TMW_COLORBAR'));
          colorbar;
       end
     end

   case 'Reset Wts'
      weights_changed = 1;
      flambda_scales = ones(size(flambda_scales));
      %   replot
      neig_plt = min(nt-1,neig_plt);
      temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
      temp = log10(temp);
      temp = [ temp temp(:,nb) ];
      handles = get(h_flambda,'Children');
      nh = length(handles);
      for kk = 1:nh
         htemp = get(handles(kk),'type');
         if htemp(1:4) =='surf'
            set(handles(kk),'Cdata',temp);
            break
         end
      end
      delete(findobj('Tag','TMW_COLORBAR'));
      axes(h_flambda);
      colorbar;
   case 'Load Params'
      [paramfile,parampath] = ...
              uigetfile('*.mat','Load Paramters from') 
      eval(['load ' parampath paramfile ])
      set(h_DAMP(1),'string',num2str(nu(1)))
      set(h_DAMP(2),'string',num2str(nu(2)))
      for k = 3:nt
         ctag = [ 'D' num2str(k)];
         set(findobj('Tag',ctag),'string',num2str(nu_comp(1,k-2)));
         ctag = [ 'd' num2str(k)];
         set(findobj('Tag',ctag),'string',num2str(nu_comp(2,k-2)));
      end
    %   replot flambda_wts
      neig_plt = min(nt-1,neig_plt);
      temp = flambda(1:neig_plt,:).*flambda_scales(1:neig_plt,:);
      temp = log10(temp);
      temp = [ temp temp(:,nb) ];
      handles = get(h_flambda,'Children');
      nh = length(handles);
      for kk = 1:nh
         htemp = get(handles(kk),'type');
         if htemp(1:4) =='surf'
            set(handles(kk),'Cdata',temp);
            break
         end
      end
      weights_changed = 1;

   case 'Save Params'
%     first make sure paramters are up-to-date
      nu(1) = str2num(get(h_DAMP(1),'string'))
      nu(2) = str2num(get(h_DAMP(2),'string'))
      for k = 3:nt
         ctag = [ 'C' num2str(k)];
         plot_ch(k) = get(findobj('Tag',ctag),'value');
         ctag = [ 'D' num2str(k)];
         nu_comp(1,k-2) = str2num(get(findobj('Tag',ctag),'string'));
         ctag = [ 'd' num2str(k)];
         nu_comp(2,k-2) = str2num(get(findobj('Tag',ctag),'string'));
      end 

      [paramfile,parampath] = ...
              uiputfile('*.mat','Save Paramters in') 
      eval(['save ' parampath paramfile ' nu nu_comp flambda_scales' ])

   case 'Save Curves'
      [paramfile,parampath] = ...
              uiputfile('*.smth','Save Smoothed Curves in') 
      curveFile = [ parampath paramfile ];
      fid = fopen( curveFile,'w');
      fprintf(fid,'%d %d \n',nt,nb)
      for k  = 1:nb
        fprintf(fid,'%e \n',Sdms.T(k));
        for j = 1:nt
           fprintf(fid,'%12.4e %12.4e %12.4e %12.4e \n',...
                       real(TF_smth(j,1,k)),imag(TF_smth(j,1,k)),...
                       real(TF_smth(j,2,k)),imag(TF_smth(j,2,k)));
        end
      end
      fclose(fid);
   case 'Separate'
      hfig_sep_menu = figure('Position',[200,200,200,100],...
     'Name','Separate Smoothed/Residual Evecs',...
     'NumberTitle','off',...
     'Tag','Separate Menu');
     %  Max # of eivenvectors : text  
     uicontrol(gcf,'Style','text',...
        'String','Max # of evecs :  ',...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position',[.05,.75,.60,.15],...
        'BackgroundColor',[0.7 0.7 0.7]);
     uicontrol(gcf,'Style','edit',...
        'String',num2str(maxCN),...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position',[.70,.75,.25,.20],...
        'CallBack','maxCN=str2num(get(gcbo,''string''));'); 
     %  Min value for significant eivenvectors : text  
     uicontrol(gcf,'Style','text',...
        'String','Min signif eval:',...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position',[.05,.50,.60,.20],...
        'BackgroundColor',[0.7 0.7 0.7]);
     uicontrol(gcf,'Style','edit',...
        'String',num2str(minCNsig),...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position',[.70,.50,.25,.20],...
        'CallBack','minCNsig=str2num(get(gcbo,''string''));');  
     uicontrol(gcf,'Style','pushbutton',...
        'String','Apply',...
        'FontWeight','bold',...
        'Units','normalized',...
        'Position',[.25,.10,.5,.25],...
        'Callback','eigSMcbk');
     case 'Apply'
        UVall;
        enable_pltsmth = 2;
        delete(findobj('Tag','Separate Menu'));
        figure('Name','POWER IN MT AND (INDEPENDENT) COHERENT NOISE MODES')
%           'NumberTitle','off');
        nfmax = nbt;
        loglog(Sdms.T(1:nfmax),real(squeeze(Sig_MT(1,1,1:nfmax))),'r')
        hold on
        loglog(Sdms.T(1:nfmax),real(squeeze(Sig_MT(2,2,1:nfmax))),'m')
        loglog(Sdms.T(1:nfmax),real(squeeze(Sig_CN(1,1,1:nfmax))),'b')
        loglog(Sdms.T(1:nfmax),real(squeeze(Sig_CN(2,2,1:nfmax))),'c') 
        fatlines(gca,2)
        legend('MT #1','MT #2','CN #1', 'CN #2')

   otherwise
     fprintf(1,'%s \n','Case Not Coded in eigSMcbk')
end
