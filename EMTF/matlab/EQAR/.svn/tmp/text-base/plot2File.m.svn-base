function [] = plot2File(plottype, plotFile)
%usage plot2File(plottype, plotFile)
%save the current figure as a type plottype ('jpg', 'eps', etc.)
%in the file plotFile

switch plottype
    case 'eps'
        eval(['print -depsc ' plotFile]); 
    case 'jpg'
        eval(['print -djpeg90 ' plotFile]);
end;