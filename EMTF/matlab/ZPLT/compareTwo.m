cfile1 =  'ELT14ab_T16.zrr';
cfile2 = 'ELT16ab_T14.zrr';
dir = '/home/gauss/egbert/project/MTproc/ELECTRODE/MT/';

Z1 = readZfile(cfile1);
Z2 = readZfile(cfile2);

ElectrodeType = {'Hungarian','Park','Booker'};
Zyx{1} = squeeze(Z1.TF(1,3,:));
Sig{1} = squeeze(Z1.SIG_S(1,1,:));
Err{1} = squeeze(Z1.SIG_E(3,3,:));
Zyx{2} = squeeze(Z2.TF(1,2,:));
Sig{2} = squeeze(Z2.SIG_S(1,1,:));
Err{2} = squeeze(Z2.SIG_E(2,2,:));
Zyx{3} = squeeze(Z2.TF(1,3,:));
Sig{3} = squeeze(Z2.SIG_S(1,1,:));
Err{3} = squeeze(Z2.SIG_E(3,3,:));

ryx = {};
pyx = {};
ryx_se = {};
pyx_se = {};
for k = 1:3
   ryx{k} = abs(Zyx{k}).^2;
   pyx{k} = (180/pi)*atan(imag(Zyx{k})./real(Zyx{k}));
   ryx_se{k} = real(Sig{k}.*Err{k});
   pyx_se{k} = (180/pi)*sqrt(ryx_se{k}./ryx{k}); 
   ryx{k} = ryx{k}.*Z1.T'/5;
   ryx_se{k} = sqrt(ryx_se{k}.*ryx{k}.*Z1.T'*4/5);
end

rho = [];
phi = [];
rho_se = [];
phi_se = [];
for k = 1:3
  rho = [rho ryx{k}];
  phi = [phi pyx{k}];
  rho_se = [rho_se ryx_se{k}];
  phi_se = [phi_se pyx_se{k}];
end

% reset plotting limits now that rho is known
[lims] = set_lims(dir,periods,rho);
lims = [10 100000 10 10000 0 90];
[hfig] = set_fig(lims);
c_title = ''


periods = Z1.T;
pltall = ones(length(periods),1);
NBT = Z1.nbt;
[rho_axes,ph_axes] = pltrhom(NBT,pltall,periods,rho,rho_se,phi,phi_se,lims,...
           c_title,hfig);


