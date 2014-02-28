function [hfig_evec] = evecPlt(ib,ivec,Sdms,SdmHd,OPTIONS,u)
%  Usage: 
%	 [hfig_evec] = evecPlt(ib,ivec,Sdms,SdmHd,OPTIONS,U)
%   plots eigenvectors numbers specified in "ivec"
%   for a single specified frequency band "ib"
%   ib = frequency band # to plot
%   ivec = list of eigenvector #s
%   Sdms = standard SDM data structure
%   SdmHd = corresponding SDM meta-data structure
%   OPTIONS = data structure with plotting options
%   U = optional data structure giving vectors to plot
%       (the eigenvectors are already in the Sdms data structure;
%       this allows various transformed or derived vectors to
%       replace the eigenvectors in Sdms.U; if this input
%       argument is present, U is plotted instead of Sdms.U

stmonitr

nvec = length(ivec);
period = Sdms.T(ib);
var = Sdms.var(:,ib);
nf = Sdms.nf(ib);
if nargin < 6
   u = Sdms.U(:,:,ib);
   if(~OPTIONS.snr_units) 
      u = diag(sqrt(var))*u;
   end
else
   if(OPTIONS.snr_units) 
      u = diag(1./(sqrt(var)))*u;
   end
end

nn = size(SdmHd.Hp);
if length(SdmHd.Hp) > 0
   Hind = reshape(SdmHd.Hp(:,1:2),nn(1)*2,1);
end

% calculate window sizes, figure scalings
pix_y = pix_x/OPTIONS.asp;
width = (pix_x+space)*nvec+extra;
height = pix_y + extra;
width_paper = 9;
height_paper = 9*height/width;

norm_width = pix_x/width;
x_step = (pix_x+space)/width;
x0 = extra/width;
y0 = .5*extra/height;
norm_height = pix_y/height;
rect_fig = [ lower_left [width height ]]; 
rect_paper = [ 1 1 width_paper height_paper ];
if(OPTIONS.snr_units) 
  figname = [ 'Band =  ' num2str(ib) ...
             '     ::   Period =  ' num2str(period)  ' sec.      '...
                          'SNR Units'];
else
  figname = [ 'Band = ' num2str(ib) ...
             '     ::   Period =  ' num2str(period)  ' sec.      '...
             ' ::   Ref rho =  ' num2str(OPTIONS.rho_ref) ];
end
hfig_evec=figure('Position',rect_fig,...
	'Name',figname,...
	'PaperPosition',rect_paper,...
	'PaperOrientation','landscape',...
	'NumberTitle','off',...
	'Tag','evec');
rect_plt = [ x0 , y0 , norm_width, norm_height ];

% loop over desired eigenvectors
axlab = [1,1];
l_label = 1;
l_Hz = 1;
ll_lim = OPTIONS.ll_lim
for k=1:nvec
   ctit = ['Eigenvector #', num2str(ivec(k))];
   if(length(SdmHd.Hp) > 0 )
      u(:,k) = chngph(u(:,k),Hind);
   end
   [uH,uE,uZ,H_sta,E_sta] = u_pair(u(:,k),SdmHd.Hp,SdmHd.Ep,...
	SdmHd.Hz,SdmHd.orient,SdmHd.decl,SdmHd.stcor,...
	SdmHd.csta,period,OPTIONS.rho_ref,OPTIONS.snr_units);
   if( k > 1 ) l_label = 0; end
   hfig = evplt_HE(rect_plt,ll_lim,sfac,uH,uE,uZ,H_sta,...
       axlab,ctit,OPTIONS.l_ellipse,l_label,OPTIONS.l_Hz);
   if(k == 1 ) 
       xtxt = ll_lim(1) + .2*(ll_lim(2)-ll_lim(1));
       ytxt = ll_lim(3) - .15*(ll_lim(4)-ll_lim(3));
       text('Position',[xtxt,ytxt],'string','Green: Magnetics',...
           'Color',[0,.7,0],'FontSize',12,'FontWeight','bold');
   elseif (k == nvec)
       xtxt = ll_lim(1) + .2*(ll_lim(2)-ll_lim(1));
       ytxt = ll_lim(3) - .15*(ll_lim(4)-ll_lim(3));
       text('Position',[xtxt,ytxt],'string','Red: Electrics',...
            'Color','r','FontSize',12,'FontWeight','bold');
   end
   rect_plt(1)  = rect_plt(1) + x_step;
   axlab = [1,0];
end
