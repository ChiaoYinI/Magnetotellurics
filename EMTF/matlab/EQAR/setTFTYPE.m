TFTYPE.TFtype = 'CC1';
TFTYPE.UseAllCh = 1;
TFTYPE.CohUseIn = [1 1 1];
TFTYPE.CohUseOut = [1,1,1];
TFTYPE.SigLevFunc{1} = 'CcovSig';
TFTYPE.SigLevFunc{2} = 'SDMSig';
TFTYPE.SigLevFunc{3} = 'SDMSig';
TFTYPE.Kin = [1:7];
TFTYPE.Kout = [8:12];
TFTYPE.title = 'All(111,111)';