%  rotates each array pw of a single Pw file, and sets up 
%  large array V, SIG_V of all TFs and covariances (all components,
%  all periods and all files
  

    if (NEW_ROT)
	V = []; SIG_V = []; 
	for k = 1:nfiles
	    [pwRot] = rotatePw(PW(k,1),PWHD(k,1),rot_ch,sing_ch,theta_rot);
	    [v,sig_v] = pwrfsite(pwRot,[1 2]);
	    V = cat(3,V,v);
	    SIG_V = cat(3,SIG_V,sig_v);
	end
    end
%  now V is a 2 x ntf x size(T)*nfiles array (e. g. 2 x 10 x 280)


    nt_rot = 2*length(rot_ch(:,1))+length(sing_ch);
    XY = ['xy'];
    l = poluse; 
    labels = [ labels ; [chname(l,1,1) XY(l) squeeze(chname(l,1,7:10))']];
    ifref = [ ifref 1 ];
    kk = 0;
    nrot = length(rot_ch(:,1));   
    nsing = length(sing_ch);   

%  there can be only one component selected for plotting so:

    chplot = [];
   for k = 1:nrot 
      for l = 1:2 
         kk = kk + 1;
         if(chuse(l,k))  
            chplot = [ kk ]; 
            labels  = [ labels ; ...
                  [ chname(l,k,1)  XY(l) squeeze(chname(l,k,7:10))' ]];
            ifref = [ ifref 0 ];
         end
      end
   end
   ifref = [ ifref 0 ];
   
   for k = 1:nsing
      kk = kk + 1;
      	if(chuse(1,k+nrot))  
         chplot = [ kk ];
         labels  = [labels; squeeze(chname(1,nrot+k,[1:2 7:10]))'];
         ifref = [ ifref 0 ];
      end
   end
   
%  set nplot to the number of files to be plotted + 1 because 
%  v1, sig_v1 are padded below

    nplot = nfiles + 1;
    labels = char(labels);

    v1 = squeeze(V(poluse,chplot,:));
    sig_v1 = squeeze(SIG_V(poluse,chplot,:));

%   v1, sig_v1 are now column vectors of size nbt*nfiles
%   reshape them into a nbt x nfiles matrix

    v1 = reshape(v1,length(T),nfiles);
    sig_v1 = reshape(sig_v1,length(T),nfiles);

%  for plotting pseudo-sections pad the the last column and row to v1
%  and sig_v1 so that the plots show all values
nm = size(v1)
    v1 = [v1,v1(:,nm(2))];
    v1 = [v1;v1(nm(1),:)];
    sig_v1 = [sig_v1,sig_v1(:,nm(2))];
    sig_v1 = [sig_v1;sig_v1(nm(1),:)];


    if(nplot == 0)
	fprintf(2,'%s','No channels selected for plotting');
    else
    plt_type
	c_title = cell(3,1);
	c_title(3,1) = cellstr( ...
	    ['Rotation Angle = ' num2str(fix(theta_rot)) ]);
	if plt_type == 'AMPH' 
% plot selected channels as amplitudes and phases
	    c_title(1,1) = cellstr('Amplitude');
	    c_title(2,1) = cellstr('Phase [deg]');
	    amp = abs(v1);
	    ph = (180/pi)*atan2(imag(v1),real(v1));
	    lims = sect_lims(amp,ph);
	    hfig_plt = sect_AB(T,amp,ph,lims,...
		    plt_type,c_title,fig_pos,labels,ifref);
	    h_plots = [ h_plots hfig_plt];
	elseif plt_type == 'DADP'
% plot the varances of amplitude and phase
	    c_title(1,1) = cellstr('Var (Amplitude)');
	    c_title(2,1) = cellstr('Var (Phase) [deg]');
	    amp = abs(v1);
	    amp_se = 2*real(sig_v1)/sqrt(2);
	    ph_se = 2*(180/pi)*real(sig_v1)./(amp*sqrt(2));
	    lims = sect_lims(amp_se,ph_se);
	    hfig_plt = sect_AB(T,amp_se,ph_se,lims,...
		    plt_type,c_title,fig_pos,labels,ifref);
	    h_plots = [ h_plots hfig_plt];
	elseif plt_type == 'REIM'
% plot real and imaginary parts of selected channels 
	    c_title(1,1) = cellstr('Real (TF)');
	    c_title(2,1) = cellstr('Imag (TF)');
	    tr = real(v1) ;
	    ti = imag(v1);
	    lims = sect_lims(tr,ti);
	    hfig_plt = sect_AB(T,tr,ti,lims,plt_type,c_title, ...
		    fig_pos,labels,ifrefs);
	    h_plots = [ h_plots hfig_plt];
	elseif plt_type == 'DRDI'
% plot the variances of real and imaginary
	    c_title(1,1) = cellstr('Var (Real(TF))');
	    c_title(2,1) = cellstr('Var (Imag(TF))');
	    tr_se = 2*real(sig_v1)/sqrt(2);
	    lims = sect_lims(tr_se, tr_se);
	    hfig_plt = sect_AB(T,tr_se,tr_se,lims,c_title, ...
		    fig_pos,labels,ifref);
	    h_plots = [ h_plots hfig_plt];
	elseif plt_type == 'LAMP'	    
	    c_title(1,1) = cellstr('Log10 (Amplitude)');
	    c_title(2,1) = cellstr('Phase [deg]')
	    amp = log10(abs(v1));
	    ph = (180/pi)*atan2(imag(v1),real(v1));
	    lims = sect_lims(amp,ph);
	    hfig_plt = sect_AB(T,amp,ph,lims,...
		    plt_type,c_title,fig_pos,labels,ifref);
	    h_plots = [ h_plots hfig_plt];
	end
    end
