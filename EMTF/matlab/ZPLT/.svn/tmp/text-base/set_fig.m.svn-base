%  set_fig:
%  sets up figure, given plotting limits for rho
%  returns figure handle
% Usage:  [hfig] = set_fig(lims)

function [hfig] = set_fig(lims,pltnum)

%  load in size factors appropriate for this system
pltfrm

%  call with two arguments to make multiple plots
if nargin ==2
  size_fac=size_fac2;
else
  size_fac = size_fac1;
end

globinc
%global sym_line_width line_styles symbol_styles rectZ rect_paper rect_rho rect_ph

paperSizeFac = .65;
sym_line_width = 1.5;
one_dec = 1.6;
xdecs = log10(lims(2)) - log10(lims(1));
one_dec = one_dec*4/xdecs;
ydecs = log10(lims(4)) - log10(lims(3));
paper_width = xdecs*one_dec;
paper_height = ( ydecs + 3 ) * one_dec;
paper_height = min([paper_height,9]);
rectZ = [0.5,0.5,paper_width,paper_height] * size_fac;
if nargin > 1
   rectZ(1) = rectZ(1) + ((0.2+paper_width)*(pltnum-1))*size_fac;
end
rect_paper = [1.,1.,paper_width*paperSizeFac,...
		paper_height*paperSizeFac];

rect_rho = [.15,.15+2.3/(ydecs+3),.8,ydecs/(ydecs+3)*.8];
rect_ph = [.15,.15,.8,2/(ydecs+3)*.8];
hfig = figure('Position',rectZ,'PaperPosition',rect_paper,...
	'Tag','MT Plot');
marg = 1.25;

%  First half are
%%%   new colors ... better for matlab 5
line_styles = ['b-';'r-';'g-';'m-';'c-';'k-';'y-';...
               'b-';'r-';'g-';'m-';'c-';'k-';'y-'];
symbol_styles = ['bo';'rx';'go';'mx';'co';'kx';'yo';...
                 'b+';'r*';'g+';'m*';'c+';'k*';'y+'];
symbol_colors = [ 0 0 .7; .7 0 0; 0 .7 0; .35 0 .35; 0 .35 .35 ; 0 0 0 ; .3 .3 0; ... 
                  0 0 .7; .7 0 0; 0 .7 0; .35 0 .35; 0 .35 .35 ; 0 0 0 ; .3 .3 0];
Nsym = 7;
