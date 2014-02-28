%  TEST SCRIPT
%   nch is total number of channels in the array
nch = 15;
%   inds1, inds2 are the indicies for channel groups  1 and 2
inds1 = [1:5];
inds2 = [6:15];
%   Pvals are the significance values to compute quantiles for
Pvals = [.5 .75 .90 .95 ];
%   Nsamp is the size of the monte-carlo sample ... increase for
%    better accuracy of the significance curves ... experiment
Nsamp = 500;
%   Nsegs is an array giving the sample sizes to do the
%   monte-carlo calculations for
Nsegs = round(10.^[1:.125:3.5]);
%   V is the "signal vector" ... here idealized to correspond to
%    a 1D Earth with equal SNR in all channels (everything is
%    non-dimensional!)  NOte there are two signal vectors here
V = [1 0;
    0 1;
    0 0;
    0 -1;
    1 0];
V = [V;V;V];
%   lambda is SNR for the two signal vectors
lambda = [20,20];
%    END OF SETUP 
%    pass arguments to CCsignif ...
[CC] = CCsignif(nch,inds1,inds2,Pvals,Nsegs,V,lambda,Nsamp);
%       ... results are returned in 
%          CC(length(Nsegs),length(Pvals),nchMin-2)
%      where nchMin is min(length(inds1),length(inds2)) = number
%        of canonical coherences which should be zero (i.e., omitting
%        the first two, which should be >> 0).

%   Plot some results
figure
semilogx(Nsegs,CC(:,1,1),'r-');
hold on
semilogx(Nsegs,CC(:,1,2),'b-');
semilogx(Nsegs,CC(:,1,3),'g-');
fatlines(gca,2)
set(gca,'FontWeight','demi','FontSize',14,'Ylim',[0,1])
xlabel('Sample Size');
ylabel('Canonical Coherence');
legend({'CC #3','CC #4','CC#5'})
title('50% Significance Level')
figure
semilogx(Nsegs,CC(:,4,1),'r-');
hold on
semilogx(Nsegs,CC(:,4,2),'b-');
semilogx(Nsegs,CC(:,4,3),'g-');
fatlines(gca,2)
set(gca,'FontWeight','demi','FontSize',14,'Ylim',[0,1])
xlabel('Sample Size');
ylabel('Canonical Coherence');
legend({'CC #3','CC #4','CC#5'})
title('95% Significance Level')
