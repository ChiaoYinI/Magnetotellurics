action = get(gcbo,'Tag');
hfig = get(gcbo,'UserData');
FigData = get(hfig,'UserData');
switch action
   case 'xMin'
      xmin = str2num(get(gcbo,'String'));
      FigData.lims(1) = xmin;
      set(hfig,'UserData',FigData);
   case 'xMax'
      xmax = str2num(get(gcbo,'String'));
      FigData.lims(2) = xmax;
      set(hfig,'UserData',FigData);
   case 'rhoMin'
      rmin = str2num(get(gcbo,'String'));
      FigData.lims(3) = rmin;
      set(hfig,'UserData',FigData);
   case 'rhoMax'
      rmax = str2num(get(gcbo,'String'));
      FigData.lims(4) = rmax;
      set(hfig,'UserData',FigData);
   case 'phiMin'
      pmin = str2num(get(gcbo,'String'));
      FigData.lims(5) = pmin;
      set(hfig,'UserData',FigData);
   case 'phiMax'
      pmax = str2num(get(gcbo,'String'));
      FigData.lims(6) = pmax;
      set(hfig,'UserData',FigData);
   case 'theta'
      theta=str2num(get(gcbo,'String'))
      FigData.theta = theta;
      MeasurementCoordinates = ( theta == FigData.th_x);
      FigData.MC = MeasurementCoordinates;
      set(hfig,'UserData',FigData);
   case 'Plot'
      rpltmtF(hfig);
   case 'Cancel'
      FigData.lims = FigData.lims_old;
      FigData.theta = FigData.theta_old;
      delete(gcf)
   case 'Quit'
      delete(gcf)
      delete(hfig)
   case 'New Plot'
      delete(gcf)
      [cfile,dir] = uigetfile(zid);
      z_mtem;
   case 'Add Band'
      delete(gcf)
      addband;
   case 'Print'
      cfile = get(gcbo,UserData);
      ind = find(c_file == '.');
      eval(['print -djpeg ' cfile(1:ind) '.jpg']);
end
