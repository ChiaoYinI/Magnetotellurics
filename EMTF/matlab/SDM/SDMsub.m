%*********************************************************************
%  SDMsub  new canonical coherence plotting function 
function [hCC] = SDMsub(Sdms,ind1,ind2)
%       given indices ind1, ind2 for two groups of channels compute
%	eigenvalues for each diagonal submarix
%	plus singular values for cross-product matrix (i.e.,
%	canonical covariances) and canonical coherences

CC = compCC(Sdms,ind1,ind2);

line_sty1 = 'r-';
line_sty2 = 'b--';
ismooth = 0;

n1 = length(ind1);
n2 = length(ind2);
nc = min(n1,n2);
nsta = Sdms.Hd.nsta;
nch = Sdms.Hd.nch;
sta = Sdms.Hd.sta;
nbt = Sdms.Hd.nbt;
csta = [];
for ista = 1:nsta
   for ich = 1:nch(ista)
       csta = [csta  sta(:,ista)];
   end
end 
chid = Sdms.Hd.chid;

TITLE1 = ['GROUP 1 =  '];
for k=1:n1-1
   TITLE1 = [ TITLE1 csta(:,ind1(k))' ' ' chid(:,ind1(k))' ' : ' ];
end
TITLE1 = [ TITLE1 csta(:,ind1(n1))' ' ' chid(:,ind1(n1))' ];
TITLE2 = [ 'GROUP 2 =  '];
for k=1:n2-1
   TITLE2 = [ TITLE2 csta(:,ind2(k))' ' ' chid(:,ind2(k))' ' : ' ];
end
TITLE2 = [ TITLE2 csta(:,ind2(n2))' ' ' chid(:,ind2(n2))' ];

period = CC.T;
stmonitr
rect_win = rect_cc;
rect_paper = [1.,1.,6.5,6.5];
rect_11 = [.1,.52,.35,.35];
rect_22 = [.55,.1,.35,.35];
rect_21 = [.55,.52,.35,.35];
rect_12 = [.1,.1,.35,.35];
rect_tit = [ .2 .93 .6 .05 ];
hCC = figure('Position',rect_win,'PaperPosition',rect_paper,...
	'Tag','CC Plot');
xmin = period(1)/1.1 ;
xmax = period(nbt)*1.1 ;
ev1 = 10*log10(CC.ev1);
ev2 = 10*log10(CC.ev2);
emax = max(max(ev1));
emax2= max(max(ev2));
emax = max([emax emax2]);
ccov = 10*log10(CC.ccov);
ymin = -10
ymax = 10*(ceil(max(emax)/10)+1);
ymax = max(ymax,10)
x_lab = (log10(xmax)-log10(xmin))*.1+log10(xmin);
x_lab = 10^x_lab;
y_lab = .85*ymax;

h11 = axes('Position',[rect_11]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ev1,3);
   for c = 1:n1 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:n1 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,ev1(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end
tit_all = 'Eigenvalues of S11';
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)

lims = [ xmin,xmax,ymin,ymax];
xt1 = ceil(log10(xmin)); xt2 = floor(log10(xmax));
xt2 = max(xt1,xt2);
xt = 10.^[xt1:1:xt2];

axis(lims);
set(gca,'FontWeight','bold','FontSize',13,'XTick',xt,...
	'TickLength',[.07,.07])
ylabel('SNR : dB');
text(-.1,1.25,TITLE1,'Units','Normalized','FontSize',8,...
            'Color','red');
text(-.1,1.15,TITLE2,'Units','Normalized','FontSize',8,...
             'Color','blue');

h22 = axes('Position',[rect_22]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ev2,3);
   for c = 1:n2 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:n2 ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,ev2(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

tit_all = 'Eigenvalues of S22';
text(x_lab,y_lab,tit_all,'FontSize',11,'FontWeight','bold');
lims = [ xmin,xmax,ymin,ymax];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13,'Xtick',xt,...
	'TickLength',[.07,.07])
xlabel('PERIOD (s)');
ylabel('SNR : dB');

h12 = axes('Position',[rect_21]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(ccov,3);
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1; 
      else
          linestyle = line_sty2;
      end
      semilogx(period,ccov(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

%   last canonical correlations
tit_all = 'Canonical Covariances' ;
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)
lims = [ xmin,xmax,ymin,ymax];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13,'Xtick',xt,...
	'TickLength',[.07,.07])
ylabel('SNR : dB');

h12 = axes('Position',[rect_12]);
if(ismooth > 0)
   period_sm = smooth(period,3);
   E_sm = smooth(CC.ccor,3);
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period_sm,E_sm(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
else
   for c = 1:nc ;
      if( c <= 2 )
          linestyle = line_sty1;
      else
          linestyle = line_sty2;
      end
      semilogx(period,CC.ccor(c,:),linestyle);
      line_handles = get(gca,'Children');
      if(c == 1) 
         set(line_handles(length(line_handles)),'LineWidth',2);
      else
         set(line_handles(length(1)),'LineWidth',2);
      end
      hold on ;
   end
end

y_lab = .85;
tit_all = 'Canonical Coherence';
text(x_lab,y_lab,tit_all,'FontWeight','bold','FontSize',11)
lims = [ xmin,xmax,0,1];
axis(lims);
set(gca,'FontWeight','bold','FontSize',13,'Xtick',xt,...
	'TickLength',[.07,.07])
xlabel('PERIOD (s)');
ylabel('Coherence ');
