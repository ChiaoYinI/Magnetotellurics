cfile = get(gcbo,'UserData');
ind = find(cfile == '.');
eval(['print -djpeg ' cfile(1:ind) 'jpg']);
