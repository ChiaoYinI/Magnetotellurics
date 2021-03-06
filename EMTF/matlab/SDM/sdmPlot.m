%   MODIFIED Jan 2003
%   MODIFIED Jan 2005
%  matlab script (not function) for looking at S0**** files
%  Run this script from the directory containing the sdm files (S0***)
%   that you want to look at
%  After choosing a file from the menu, a plot of sdm eigenvalues
%  vs. period will appear.  There are two additional sorts of plots
%  so far: (1) there is a slider bar beneath the period axis on the
%   eigenvector plot ... with the slider, and the "PLOT" pushbutton
%   you can plot the eigenvectors.  These appear as a plot of
%  magnetic and electric vectors on a map of station locations 
%  (note that this assumes station locations are in the sdm files!)
%  (2) you can plot eigenvalues and canonical coherences for subgroups
%  of data channels.  You get to this via the pull down menu called
%   "Plot Options".  A set of pushbuttons listing all components in 
%  the array pops up ... you use this to pick the components in "group 1";
%  the remainder are in group 2.  The program then plots eigenvalues
%  for each group separately + canonical coherences and covariance.
%   If there is coherent noise present which only occurs at some sites,
%   or which only occurs in E components (e.g.), then this should
%   help you figure out which sites/components are not contaminated.

global SYST
hfig_evec = [];
clear h_ib;

%   Get name of SDM file ...
%    allows standard *.S0 files output by multmtrn and matfiles
%        (contianing appropriate structures!)
[cfile, cpath, filterindex] =     ...
           uigetfile({'*.S0','S0 binary files'; ...
        '*.mat','SDM mat file'},'SDM File to Plot');
eval(['cd ' cpath])

switch filterindex
  case 1
     arrayid = cfile(1:length(cfile)-3);
     s0_file = [ cpath cfile ];
     if (isempty(SYST))
        [fid_sdm,irecl,nbt,nt,nsta,nsig,nch, ...
		ih,stcor,decl,sta,chid,csta, ...
                         orient,periods,ndf] = sdm_init(s0_file);
     else
        [fid_sdm,irecl,nbt,nt,nsta,nsig,nch, ...
		ih,stcor,decl,sta,chid,csta, ...
                         orient,periods,ndf] = sdm_init(s0_file,SYST);
     end
     ch_name = chid(1:2,:);
     SdmHd = struct('nbt',nbt,'nt',nt,'nsta',nsta,'nch',nch,'ih',ih,...
     'stcor',stcor,'decl',decl,'chid',chid,'sta',sta,'orient',orient,...
     'ch_name',ch_name);
     [Sdms,err] = loadsdms(fid_sdm,irecl,nbt,nt);
     Sdms.Hd = SdmHd;
  case 2
     arrayid = cfile(1:length(cfile)-3);
     savg_file = [ cpath cfile ];
     eval(['load ' savg_file]);
     if exist('Sdmhd') == 1
        SdmHd = Sdmhd;
     end
     if isfield(Sdms,'Hd')
         SdmHd = Sdms.Hd;
     else
        Sdms.Hd = SdmHd;
     end
     ndf = Sdms.nf;
end

periods = Sdms.T';
nbt = length(periods);
nch = SdmHd.nch;
nsta = SdmHd.nsta;
nt = SdmHd.nt;
ih = SdmHd.ih;
stcor = SdmHd.stcor;
decl = SdmHd.decl;
chid = SdmHd.chid;
sta = SdmHd.sta;
orient = SdmHd.orient;
ch_name = SdmHd.ch_name;
csta = [];
for k=1:nsta
  for l=1:nch(k)
    csta = [csta [sta(1:3,k)]];
  end
end
SdmHd.csta = csta;

Uplt = Sdms.U;
for ib = 1:nbt
  N = sqrt(Sdms.var(:,ib));
  Uplt(:,:,ib) = diag(N)*Uplt(:,:,ib);
end

stn = [];
for ista = 1:nsta
   stn = [ stn ista*ones(1,nch(ista)) ];
end

%  rho_ref is the apparent resistivity assumed for scaling E-field vectors
%  into nT ... if actual apparent resistivity = rho_ref, then H and E vectors
%  will have the same length on eigenvector plots 
%  ivec is an array which tells which eigenvectors of the SDM should
%  be plotted for each band; 1 = eigenvalue associated with largest eigenvector, etc.
%           ( both of these should ultimately be changeable from menu ...)
rho_ref = 100;
n_evec = 4; ivec = [1:n_evec];
snr_units = 0; l_ellipse = 0;l_Hz = 0;
l_PWG = 0; l_impNorm = 0;

title_sig = [arrayid ': Signal Power'];
title_noise = [arrayid ': Noise Power'];
l_smthsep = 0;
l_perp = 0;
l_MTTF = 0;
enable_pltsmth = 1;
maxCN  = 2;
minCNsig = 1.5;

[SdmHd.Hp,SdmHd.Ep,SdmHd.Hz] = ch_pair(nch,chid);

% some calculations for scaling plot sizes ...
lat_max =max(max(stcor(1,:)));
lat_min =min(min(stcor(1,:)));
lon_max =max(max(stcor(2,:)));
lon_min =min(min(stcor(2,:)));
lat_av = (lat_max+lat_min)/2.;
lat_range = (lat_max-lat_min);
lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);

marg = 1.5/(sqrt(nsta))*max([lat_range,lon_range]);
 if(marg == 0)
%&*&^(*(7%$%$#  user forgot to define site locations in sp file
%              ...  yet again !!!
%   so define dummy site locations ...
   nn = size(stcor);nsta = nn(2);
   stcor = [ 1:nsta; 1:nsta ];
   lat_max =max(max(stcor(1,:)));
   lat_min =min(min(stcor(1,:)));
   lon_max =max(max(stcor(2,:)));
   lon_min =min(min(stcor(2,:)));
   lat_av = (lat_max+lat_min)/2.;
   lat_range = (lat_max-lat_min);
   lon_range = cos(lat_av*pi/180)*(lon_max-lon_min);
   marg = .6*max([lat_range,lon_range]);
 end
if (marg == 0), marg =1 ; end

asp = (2*marg+lon_range)/(2*marg+lat_range);
lat_max = lat_max+marg;
lon_max = lon_max+marg;
lat_min = lat_min-marg;
lon_min = lon_min-marg;
ll_lim = [ lat_min,lat_max,lon_min,lon_max];
EVEC_OPTIONS = struct('asp',asp,'rho_ref',rho_ref,'n_evec',n_evec,...
	'ivec',ivec,'snr_units',snr_units,'l_ellipse',l_ellipse,...
	'll_lim',ll_lim,'l_Hz',l_Hz,...
	'l_PWG',l_PWG,'PWGnum',5,'PWGevalWt',1,'l_impNorm',l_impNorm);
sdmPlotSet
