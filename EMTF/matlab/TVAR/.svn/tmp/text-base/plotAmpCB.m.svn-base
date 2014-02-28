action = get(gcbo,'Tag');
hFig = gcf;
OPTIONS = get(hFig,'UserData');
nSeg = length(OPTIONS.segments);
k = OPTIONS.kSeg;
if(OPTIONS.plotAmp ==1)
    ResSig = 'Sig'; 
else
    ResSig = 'Res';
end

switch action
    case 'Back'
        if k > 1
            delete(hFig);
            k = k-1;
            OPTIONS.kSeg = k;
            plotAmp(ResSig,Sig{k},Res{k},OPTIONS);
        end
    case 'Forward'
        if k < nSeg
            delete(hFig);
            k = k+1;
            OPTIONS.kSeg = k;
            plotAmp(ResSig,Sig{k},Res{k},OPTIONS);
        end   
    case 'Select Segment'
        k = get(gcbo,'Value');
        delete(hFig);
        OPTIONS.kSeg = k;
        plotAmp(ResSig,Sig{k},Res{k},OPTIONS);
    case 'Residual'
        delete(hFig);
        ResSig = 'Res';
        OPTIONS.plotAmp = 0;
        OPTIONS.plotRes = 1;
        plotAmp(ResSig,Sig{k},Res{k},OPTIONS);
    case 'Signal'
        delete(hFig);
        ResSig = 'Sig';
        OPTIONS.plotAmp = 1;
        OPTIONS.plotRes = 0;
        plotAmp(ResSig,Sig{k},Res{k},OPTIONS);
    case 'Print'
end
