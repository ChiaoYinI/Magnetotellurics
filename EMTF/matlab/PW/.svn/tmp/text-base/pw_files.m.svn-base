%	pw_files		MATLAB function
%
%	pw_files prompts the user for a file with a list of Pw files to
%	plot with Pw_sect

    function [pwfiles,nfiles] = pw_files;
    pwfiles = [];
    nfiles = 0;
    [cfile, cdir] = uigetfile('*.pwi');
    cfile = [cdir cfile]
    fl_id = fopen(cfile,'r');
    cline = fgets(fl_id);
    while (cline ~= -1)
	pwfiles = [pwfiles; cline(1:length(cline)-1)];
	nfiles = nfiles + 1;
	cline = fgets(fl_id);
    end
    fclose(fl_id);
%
%	assuming that all file names are relativ to the directory
%	of the pwi file, cd in to that directory
%    cd cdir;

