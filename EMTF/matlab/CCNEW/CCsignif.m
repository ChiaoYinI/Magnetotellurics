function [CC] = CCsignif(nch,inds1,inds2,Pvals,Nsegs,V,lambda,Nsamp)
% Usage: [CC] = CCsignif(nch,inds1,inds2,Pvals,Nsegs,V,lambda,Nsamp);
%   CCsignif computes significance levels for canonical coherences
%    for Gaussian random vectors.   To be precise the statistical
%    model for the K dimensional data vector X is:
%     X = sum_i=1^K a_i V_i + e
%     where K is the dimension of the signal space, which is
%     defined by the columns of the nch x K matrix V
%     and e is an nch-dimensional vector of iid unit variance (zero
%     mean) Gaussian errors.   The coefficients a_i are also
%     random, and Gaussian, with variances given by the input vector
%     lambda.  All random vectors (a, e) are complex (and the variances
%     given correspond to complex variances, the expected value of the
%     sum of real and imaginary parts of the random vectors).
%    Input arguments lambda, V, thus define the model
%    Other inputs:
%       inds1, inds2  define the two groups of channels for the
%           canonical coherence analysis
%       Pvals = list of one or more probability values (e.g., .95)
%                to compute significance levels for
%       Nsegs = list of sample sizes (i.e., number of complex vectors
%                  X used to form the cross-product matrix used for
%                  CC analysis
%       Nsamp = number of times to generate a random SDM and do CC
%                 analysis for the monte-carlo simulation used to
%                 estimate significance levels
%      CC =  estimated significance levels, a 3-D array.
%      CC(i,j,k) gives the significance level for:
%                    i:  sample size Nseg(i)
%                    j:  probability value Pval(j)                      
%                    k:  CC # j+K  (i.e., if K = size of assumed signal 
%                           space = 2, j = 1 corresponds to CC#3 ...
%                            the largest of those that should be zero
%                            for infinite sample size.
%     The significance level gives the CC value that will be exceeded
%      a fraction Pval of the time, if the statistical model
%       (Gaussian random vectors, K signal vectors) holds
%     The idea is to use this to make a table of significance levels
%       for commonly used CC configurations (e.g., 10  channels from
%        two sites in one group, 5 channels from a third site in the
%        other).  The table can be saved, and used (with interpolation)
%        to compute significance levels for actual sample sizes, so
%        that a significance line (e.g., 95%) can be added to CC plots
%     I hypothesize (but have not tested extensively) that the tables
%       are insensitive to the actual vectors V in the model above.

N = length(Nsegs);
M = length(Pvals);
K = length(lambda);
n1 = length(inds1);
n2 = length(inds2);
nCC = min(n1,n2);
L = nCC-K;
CC = zeros(N,M,L);
ccor = zeros(nCC,Nsamp);
%   loop over sample sizes
for n = 1:N
    %   generate random SDMs
    [S] = mkRandSDM(nch,Nsegs(n),Nsamp,V,lambda);
    %  compute CCs for each sample
    for ib = 1:Nsamp
       S11 = squeeze(S(inds1,inds1,ib));
       S22 = squeeze(S(inds2,inds2,ib)); 
       S12 = squeeze(S(inds1,inds2,ib));
       temp = real(eig((S12/S22)*S12',S11));
       temp = sort(temp);
       ccor(:,ib) = temp(end:-1:1);
    end
    %  for each canonical correlation (among those supposed to be zero)
    %   compute sample quantiles for requested P-values
    for l = 1:L
        temp = sort(ccor(l+K,:));
        pInd = round(Pvals*Nsamp);
        CC(n,:,l) = temp(pInd);
    end
end
