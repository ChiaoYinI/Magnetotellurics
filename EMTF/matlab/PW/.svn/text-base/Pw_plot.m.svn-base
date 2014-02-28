%  This script reads in reuslts output in the array TF
%  file Pw**** output by multmtrn, and plots any chosen interstation and/or
%  intercomponent transfer functions desired, in any chosen coordinate
%  system

%  GLOBAL VARIABLES
%  channels chosen for plotting, rotation are global variables
%  ref_ch  =  reference channels (currently has to be an actual channel #)
%  rot_ch  =  channel numbers for rotation pairs
%  sing_ch =  channel numbers for single (unpaired for rotation) channels
%  plot_ch =  channels to plot TFs for
%  theta_rot = rotation angle
%  poluse = polarization to plot (1 = Hx , 2 = Hy)
%  chplot = channel number in rotated channel list
global ref_ch sing_ch rot_ch plot_ch theta_rot poluse chplot

%  all of these are associated with dialogue box which controls
%  channel rotations etc.
global chnum chuse chname chhand name_list nlines plt_type

%  default plotting types
ptyps = ['AMPH';'LAMP';'REIM'];
plt_type = ptyps(1,:);
comp=computer;

if ~exist('PWfilt') PWfilt = []; end
if (isempty(PWfilt)) PWfilt = '*.Pw'; end
%   get pathname for Pw file to plot
if comp(1:4)=='PCWI',
 PWfilt=strrep(PWfilt,'/','\');
end
[cfile,cpath] = uigetfile(PWfilt,'Select MMT PW File');
cfile = [cpath cfile]; 
if(cfile ~= 0)

%  load in Pw file and header into Pw and Pwhd data structures
[pw,pwhd] = pwstruct(cfile);

% open dialog box for choosing rotation pairs, polarizations,
[h_pwdialog] = pwdialog(pwhd);

h_plots = [];
end
