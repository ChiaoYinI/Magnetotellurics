function [] = rpltmtF(hfig);
%  mtrplt  replots MT figure ... recoded as a function,
%  which gets all needed plotting info from a structure stored
%  in "UserData" of figure hfig

   FigPos = get(hfig,'Position');
   FigData = get(hfig,'UserData');

   if 1-min(FigData.lims_old == FigData.lims)
      %  set up figure with new limits
      cfile = get(findobj('Tag','Print','Parent',hfig),'UserData');
      delete(hfig);
      [hfig] = set_fig(FigData.lims);
      set(hfig,'Position',FigPos);
      uimenu('Parent',hfig,'Label','Plot Options',...
         'Callback','mt_menu')
      uimenu('Parent',hfig,'Label','Print','Tag','Print',...
         'Callback','printCB','UserData',cfile)
   else
      delete(FigData.rho_axes);
      delete(FigData.ph_axes);
   end   

   % transform/rotate
   [ryx,rxy,pyx,pxy,ryx_se,rxy_se,pyx_se,pxy_se] = ...
                  xform(FigData.Z,FigData.Sig_s,...
		FigData.Sig_e,FigData.periods,...
		FigData.orient,FigData.xypairs,FigData.theta);
   rho = [ rxy   ryx ];
   ph = [ pxy   pyx ];
   rho_se = 2*[ rxy_se   ryx_se ];
   ph_se = 2*[ pxy_se   pyx_se ];
         
%  replot in any case
   NBT = length(FigData.periods);
   [rho_axes,ph_axes] = pltrhom(NBT,FigData.pltind,...
		FigData.periods,rho,rho_se,ph,ph_se,FigData.lims,...
                FigData.c_title,hfig);

%  copy new data into FigData
   FigData.rho_axes = rho_axes;
   FigData.ph_axes = ph_axes;
   FigData.old_lims = FigData.lims;
   FigData.old_theta = FigData.theta;
   FigData.pos = FigPos;
   set(hfig,'UserData',FigData);
