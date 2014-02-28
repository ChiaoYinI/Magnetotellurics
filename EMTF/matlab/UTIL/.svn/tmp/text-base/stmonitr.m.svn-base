% include file to set values of paramters appropriate to
%  running on specific computers (mostly all positions/sizes
%  of figure windows, to accomodate different computer screens
%   FOR SDM_PLOT

l_scale = [ 1, 1, .95,.9,.85,.8,.75,.75,.75,.70,.70,.70];
ctemp = computer;
if(ctemp(1:4)  == 'PCWI' )
%  ASSUME THIS IS THE DELL NOTEBOOK

   undobutton = 'alt';
   %  eval_plt (main window)
   size_fac2 = 50;
   size_fac1 = 60;
   rect_win =  [20,20,400,520];
   
   %  plt_snr  (Signal-to-noise plot)
   scale_snr = 1;
   rect_snr = [100,20,450,450];
   
   %  canonical coherences
   rect_cc =  [200,20,450,450];
   width_cc = 125;
   height_cc_1 = 25;
   loc_cc  = [400,50];
   
   %  evec_plt  (eigenvector plots)
   lower_left = [50,50];
   pix_x = 180;
   extra = 100;
   space = 25;
   sfac = .6;

   line_thick = 1.25;
   
   %  eigSM
   rect_main = [10,50,780,500];
   leg_scale = .75;
   
else
%  ASSUME A SUN workstation with a bigger screen

   undobutton = 'ext';
   %  eval_plt (main window)
   size_fac2 = 80;
   size_fac1 = 90;
   rect_win =  [5,200,450,600];
   
   %  plt_snr (Signal-to-noise)
   scale_snr = .9;
   rect_snr = [550,100,550,800];
   
   %  canonical coherences
   rect_cc =  [50,50,600,600];
   width_cc = 150;
   height_cc_1 = 25;
   loc_cc  = [660,50];
   
   %  evec_plt  (eigenvector plots)
   lower_left = [100,10];
   pix_x = 200;
   extra = 100;
   space = 25;
   sfac = .6;

   line_thick = 2;

   %  eigSM
   rect_main = [10,10,900,800];

   leg_scale = 1;
end
