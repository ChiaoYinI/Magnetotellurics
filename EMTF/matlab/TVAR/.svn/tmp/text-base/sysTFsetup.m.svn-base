function [sysTF] = sysTFsetup(sysTFfiles)
%   Usage:   [sysTF] = sysTFsetup(sysTFfiles);
%   could update to check for Timing error file, and modify
%    system response accordinglya : see rdTS.m
nFiles = length(sysTFfiles);
TF = [];
for k = 1:nFiles
   [sysTF] = rdSysTF(sysTFfiles{k});
   if ~isstruct(sysTF)
      TF = [];
      return
   else
      TF = [TF ; sysTF.TF];
   end
end
sysTF.TF = TF;
sysTF.PhysU = 1;
