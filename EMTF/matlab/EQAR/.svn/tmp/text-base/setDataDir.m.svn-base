function [dataDir] = setDataDir(DataEnvirVbl);
%  sets data directory from data directory environment variable 
%   (i.e., PKDSAOdata as we have been doing this)
%  unix version ... could generalize, somehow
unixCmd = ['echo $' DataEnvirVbl];
[s,dataDir] = unix(unixCmd);
%  strip off EOL
dataDir = [dataDir(1:end-1) '/'];
%dataDir = '/home/gauss/scratch/PKDSAOg/data/';
