function [RD] = resStat(RES,Pval,TimeDiv)
%   Usage : [RD] = resStat(Res,Pval,TimeDiv);
%   Computes quantiles of distribution of residual amplitudes
%
%   INPUTS:  
%       RES{nSeg} = cell array of FT data structures containing
%          (time/frequency band averaged) residual amplitudes for a
%          series of processed segments, each resulting from a run
%          of one of the "resplt" functions (e.g., resplt_logf.m)
%           NOTE:  all FT structures are assumed to have the same channels
%            ferequencies, etc.
%       Pval(nP) = array of p-values (between 0 and 1) that quantiles
%          of residual distributions will be computed for 
%       TimeDiv = structure used to define time division, consisting of
%          two fields: nDiv = number of divisions, and func = name of a
%          matlab function that maps the array of residual times Res{k}.t
%          to an array of time groupings  (i.e., an integer in the range 1:nDiv)
%          A typical case would be to sort the data by time of day, in 2 hour
%          bins, but other divisions could be used by changing the function named
%          in TimeDiv.func.
%          NOTE: this argument is optional; if missing nDiv=1 is assumed and 
%           distibutions are computed over all time indicies
%
%   OUTPUTS:
%       RD(nDiv,nCh,nFreq,nP) = quantiles of residual distributions for each
%          time division, channel and frequency.

nSeg = length(RES);
nP = length(Pval);
nCh = RES{1}.nch;
nFreq = length(RES{1}.f)
if nargin == 2
    nDiv = 1;
    tDiv = 0;
else
    nDiv = TimeDiv.nDiv;
    tDiv = 1;
end

RD = zeros(nDiv,nCh,nFreq,nP);

%   might reorder loops here to get faster performance (with greater memory requirements)
for iDiv = 1:nDiv
    qInd = [];
    for iCh = 1:nCh
        for iFreq = 1:nFreq
            %  accumlate residuals for channel ich, time division iDiv, iFreq
            if tDiv
                temp = [];
                for k = 1:nSeg
                    eval(['ind = ' TimeDiv.func '(RES{k}.t,iDiv);']);
                    temp = [temp squeeze(RES{k}.data(iCh,ind,iFreq))];
                end
            else
                temp = squeeze(Res{k}.data(iCh,:,iFreq));
            end
                
            % comute quantiles ... first sort
            temp = sort(temp);
            if(isempty(qInd))
                %  (only need to compute quantile indices once for each segment)
                qInd = round(Pval*length(temp));
            end
            RD(iDiv,iCh,iFreq,:) = temp(qInd);
        end
    end
end
