TFTYPE.TFtype = 'CC1';
TFTYPE.UseAllCh = 1;
TFTYPE.CohUseIn = [1 1 1];
TFTYPE.CohUseOut = [1,1,1];
TFTYPE.SigLevFunc{1} = 'CcovSig';
TFTYPE.SigLevFunc{2} = 'SDMSig';
TFTYPE.SigLevFunc{3} = 'SDMSig';
%   set up for PKD ...
if site == 'PKD'
   KinV = [8:12];
   KoutV = [1:7];
   KallV = [KinV KoutV];
else
   KinV = [1:7];
   KoutV = [8:12];
   KallV = [KinV KoutV];
end

k = 1;
ir = max(find(char(FCfiles{k})=='.'));
root = FCfiles{k}(1:ir-1);
ir = min(find(root=='_'));
root = root(ir+1:end);
FileNameSDM = [OPTIONS.dir '/SDM/' root '.S0'];

[TFres] = TFsdmSetX(FileNameSDM,KoutV,KinV,TFTYPE);

K = [];
for k = 1:length(Kout)
   ind = find(KallV == Kout(k));
   K = [K ind];
end
TFres.Kout = K;
