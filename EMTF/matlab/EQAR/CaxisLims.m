%   This script is used to set color axis limits for apparent resistivity
%    and phase plots for the PKD and SAO sites.  Different limits are allowd
%    for 1 hz and 40 hz, for xy and yx modes, and for the 100/200 meter dipoles
%    at PKD.  Need to set BAND = 1 or 40 before calling
%   Adjust these number to make plots clearer

if BAND == 40
   % limts for xy, yx modes, for site 1 (PKD 100 m)
   rhoLims{1} = [.8,2.;.1,1.3];
   % limts for xy, yx modes, for site 2 (SAO)
   rhoLims{2} = [2,3.5;2,3.5];
   % limts for xy, yx modes, for site 1 (PKD 200 m)
   rhoLims{3} = [.8,2.;1.2,2.8];
   %   phase limits for same three cases
   phiLims{1} = [20,50;10,55];
   phiLims{2} = [10,55;10,55];
   phiLims{3} = [10,55;10,55];
else
   % limts for xy, yx modes, for site 1 (PKD 100 m)
   rhoLims{1} = [1.2,2.2;.6,1.6];
   % limts for xy, yx modes, for site 2 (SAO)
   rhoLims{2} = [2,3.5;2,3.5];
   % limts for xy, yx modes, for site 1 (PKD 200 m)
   rhoLims{3} = [1.2,2.2;1.7,3.0];
   %   phase limits for same three cases
   phiLims{1} = [10,55;10,55];
   phiLims{2} = [10,55;10,55];
   phiLims{2} = [10,55;10,55];
end
if subMed == 1
   % reset limits for deviations from time avg. ... just specify
   %   fractional error in impedance; the rest is computed
   fracErr = .05;
   %  following are computed from fracErr for rho, impedance phase
   %   in degrees
   rhoLims{1} = [-fracErr,fracErr;-fracErr,fracErr];
   phiLims{1} = rhoLims{1}*180/pi;
   rhoLims{1} = rhoLims{1}*2;
   rhoLims{2} = rhoLims{1};
   rhoLims{3} = rhoLims{1};
   phiLims{2} = phiLims{1};
   phiLims{3} = phiLims{1};
end
