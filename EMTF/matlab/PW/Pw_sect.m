%	Pw_sect		MATLAB script
%
%	Pw_sect reads a set of Pw files and plots selected TFs
%	as a function of time/location. Pw_sect uses pwstruct to
%	read in the TFs of the specified files. When choosing a 
%	site as a reference it is required that this site is at
%	the same position in the sequence of sites in the Pw files.
%	It is also required that the number of channels is the same
%	in all Pw files.
%
%	Used scripts/functions:
%
%	pw_files: 	reads a *.piw file with list of Pw files 
%			<- cfile(nfiles,:),nfiles
%	pwstruct:	reads Pw files <- struct pw,pwhd
%	pwsect_dlg:	opens dialog window for selecting TFs
%			<- h_dialog (handle to figure)
%	new_sect:	callback for PLOT in dialog
%	plot_sect:	calculates selected TFs and plots
%	rotatePw:	rotates array PW of a single Pw file
%			(called by plot_sect)

global ref_ch sing_ch rot_ch plot_ch theta_rot poluse chplot
global chnum chuse chname chhand name_list nlines plt_type
global time_fact time_offs time_str
    time_fact = 10;
    time_offs = 50;
    time_str = 'Julian Day (1997)';
    h_plots = [];
    ptyps = ['AMPH';'LAMP';'REIM';'DADP';'DRDI'];
    plt_type = ptyps(1,:);
    NEW_ROT = 1
    theta = 0.0;	% rotate into geographic north

    [cfile,nfiles] = pw_files;


    %  make one large struct array for both pw and pwhd
    [PW,PWHD] = pwstruct(cfile(1,:));

    for k = 2: nfiles
		[pw,pwhd] = pwstruct(cfile(k,:));
		PW = [PW;pw];
   	PWHD = [PWHD;pwhd];
    end
%	[rot_ch, sing_ch] = rot_setup(pwhd);
%	[pwrot] = rotatePw(pw, pwhd, rot_ch, sing_ch, theta);
% now the sequence of channels has changed. It is now in the sequence
% given by the rot_ch, sing_ch arrays.
%	TF = cat(1,TF,pwrot.tf);
%	XXINV = cat(1,XXINV,pwrot.xxinv);
%	COV = cat(1,COV,pwrot.cov);

%  call a dialog box similar to the PW_PLOT dialog
%  Only the callbacks for plotting and adding are 
%  different

	[h_pwdialog] = pwsect_dlg(PWHD(1,1));
    
%    TF1 = [];
%    ref_ch = [1 2];
%    [nfiles,nch,nfreq] = size(TF);
%    A = zeros(nfiles,nfiles);
%    nfiles = nfiles/2;
%    for k = 1: nfreq
%	for l = 1:nfiles
%	    i1 = (l-1)*2+1;
%	    i2 = i1 + 1;
%	    A(i1:i2,i1:i2) = inv(TF(i1:i2,ref_ch,k));
%	end
%	TF1 = cat(1,TF1,A*TF(:,:,k));
%    end

%    figure;
%    semilogx(pw.T,abs(TF1(1:24:672,5)),'rx')
