UD = struct('Sdms',Sdms,'SdmHd',SdmHd',...
	'EVEC',EVEC_OPTIONS,'arrayid',arrayid);
ib = 1;
ismooth = 0;
E = [Sdms.T';Sdms.lambda];
hfig_eval = plt_eval(E,ndf,ismooth,arrayid);
set(hfig_eval,'Tag','SDM Plot','UserData',UD)
h_menu = uimenu('Parent',hfig_eval,'Label','Signal/Noise');
         uimenu(h_menu,'Label','Signal Power',...
                'Tag','Signal Power','Callback','sdmPlotCB')
         uimenu(h_menu,'Label','Noise Power',...
                'Tag','Noise Power','Callback', 'sdmPlotCB')
         uimenu(h_menu,'Label','SNR',...
                'Tag','SNR','Callback','sdmPlotCB')
         uimenu(h_menu,'Label','Correlation',...
                'Tag','Correlation','CallBack','sdmPlotCB')
         uimenu(h_menu,'Label','Can. Coh',...
                'Tag','Canonical Correlation','Callback','sdmPlotCB')
         h_menu_ev_opt = uimenu('Parent',hfig_eval,...
           'Label','Eigenvectors','Tag','Eigenvectors',...
           'Callback','sdmPlotCB');
         uimenu(h_menu,'Label','Print',...
                'Tag','Print','Callback','sdmPlotCB')
         uimenu(h_menu,'Label','Quit',...
                'Tag','Quit','Callback','sdmPlotCB')
h_NoiseEst = uimenu('Parent',hfig_eval,'Label','Ch. Group.');
         uimenu(h_NoiseEst,'Label','Sites',...
                'Tag','Group: Sites','Callback','sdmPlotCB')
         uimenu(h_NoiseEst,'Label','Channels',...
                'Tag','Group: Channels','Callback','sdmPlotCB')
         uimenu(h_NoiseEst,'Label','Sites+E/H',...
                'Tag','Group: Sites+E/H','Callback','sdmPlotCB')
         uimenu(h_NoiseEst,'Label','Custom',...
                'Tag','Group: Custom','Callback','sdmPlotCB')

slide_min = log10(periods(1));
slide_max = log10(periods(nbt));
slide_val = ceil(slide_min);
h_slider = uicontrol('Parent',hfig_eval,...
   'Style','slider',...
   'Units','normalized',...
   'Position',[.13,.30,.84,.03], ...
   'Max',slide_max,...
   'Min',slide_min,...
   'Value',slide_val, ...
   'Tag','GETIB',...
   'CallBack','evecCbk');
h_evec_plot = uicontrol('Parent',hfig_eval,...
   'Style','pushbutton', ...
   'Units','normalized',...
   'Position',[.025,.29,.10,.04],...
   'string','PLOT',...
   'Tag','PLOT',...
   'CallBack','evecCbk');
