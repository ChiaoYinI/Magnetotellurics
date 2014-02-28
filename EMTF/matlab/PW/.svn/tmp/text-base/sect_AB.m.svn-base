%***********************************************************************
% 	sect_AB : PLOTS a pseudo-section of two matrices A,B
%	(e. g. amplitude/phase, real/imag.)  
%   	on log or linear z-axes, y-axes is log10(period), xaxis is 
%	time/distance (linear)
%  	lims = [xmax,xmin,ymax,ymin] are plot limnits; 
%	empty array => automatic
%  	USAGE [hfig] = sect_amp_ph(T,amp,amp_se,ph,ph_se,lims,...
%				ll,c_title,pos_rect,labels,ifref)
%   returns figure handle in hfig
%   Inputs : T(nbt,1) = periods
%            A(nbt,nplot) = amplitudes, real
%            B(nbt,nplot) = phases, imag
%            lims(6) = plotting limts for T, amp, ph
%            c_title = plot title (character string)
%            pos_rect = position of figure on screen
%            lables(nplot+nref,6) = strings to label each curve with
%					(including reference channels)
%			    ifref = 1/0 array = 1 for reference labels

function hfig = sect_AB(T,A,B,lims,...
   plt_type,c_title,pos_rect,labels,ifref)

global time_fact time_offs time_str
pwstmon;
%  need this for label color
line_styles = ['b-';'r-';'g-';'c-';'m-';'k-';'y-';'b-';'r-';'g-';'c-';'m-'];

nplot = length(A(1,:));

%	Set up some figure parameters

paper_width = 6 * log10(nplot);
pos_rect(3) = pos_rect(3) * log10(nplot);
paper_height = 7.5;
width_n = 0.75;
height_B_n = .35;
height_A_n = .35;
x0_n = .1;
y0_B_n = .1;
y0_A_n = .55;

x0_cbr = x0_n + width_n + 0.07;
rect_cbr_A = [x0_cbr,y0_A_n,0.05,height_A_n];
rect_cbr_B = [x0_cbr,y0_B_n,0.05,height_A_n];

rect_paper = [.5,.5,paper_width,paper_height];
rect_A = [x0_n,y0_A_n,width_n,height_A_n];
rect_B = [x0_n,y0_B_n,width_n,height_B_n];
hfig = figure('Position',pos_rect,'PaperPosition',rect_paper);

%  some specials to plot the full matrix
TT = log10(T(1)) + (0:1:length(T))* ...
(log10(T(length(T)))-log10(T(1)))/length(T);
%TT = (TT);
%  set time_fact in Pw_sect.m
XX = (0:1:nplot-1) * time_fact + time_offs;

%       plot upper section (amplitudes/real)
    A_axes = axes('position',rect_A) ;

    pcolor(XX,TT,A);
    set(A_axes,'CLIM',lims(1:2));

    pltlabl;

    hold off ;
    set(gca,'FontSize',12,'FontWeight','bold');
    y_pos = TT(length(TT)) + (TT(length(TT))-TT(1))/5;
    text(XX(1),y_pos,char(c_title(3,1)), ...
	'FontSize',12,'FontWeight','demi');
    ylabel('log (Period) [s]');
    set(get(gco,'ylabel'),...
	'FontSize',14,...
	'FontWeight','bold',...
	'Interpreter', 'none');
    title(char(c_title(1,1)));
    set(get(gco,'Title'),...
	'FontSize',14,...
	'FontWeight','bold',...
	'Interpreter', 'none');

%      plot lower section (phases/imag)
    B_axes = axes('position',rect_B);
    pcolor(XX,TT, B);
    set(gca,'CLIM',lims(3:4));
    set(gca,'FontSize',12,'FontWeight','bold')
    title(char(c_title(2,1)));
    set(get(gco,'Title'),...
	'FontSize',14,...
	'FontWeight','bold',...
	'Interpreter', 'none');

    ylabel('log (Period) [s]');
    xlabel(time_str);
    hold off;
%
%  Make the plot nice
%
    set(A_axes,'FontWeight','bold');
    set(B_axes,'FontWeight','bold');

%  colorbars: get the range from the limits, make a linear color scale

    div = (lims(2)-lims(1))/32;
    vcbr_lin(rect_cbr_A,lims(1:2),div,'default',' ',10,'demi');
    div = (lims(4)-lims(3))/32;
    vcbr_lin(rect_cbr_B,lims(3:4),div,'default',' ',10,'demi');


