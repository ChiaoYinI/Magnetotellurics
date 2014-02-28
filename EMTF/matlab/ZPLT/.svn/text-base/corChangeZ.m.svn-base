function [Z] = corChangeZ(Z,theta0)
%
%  Usage: [Z] = corChangeZ(Z,theta0);
%
%  change coordinates of transfer function structure
%   to fixed coordinate system, with x-axis oriented
%   theta degrees E of N.
%
%  limited generality, assumes 4 or 5 channel standard MT 
%  with two horizontal H channels first.  If nch = 5 assumes
%  channel # 3 is Hz
%

if nargin ==1
   theta0 = 0;
end

orient = Z.orient;
theta1 = Z.orient(1,1)+Z.theta0-theta0;
theta2 = Z.orient(1,2)+Z.theta0-theta0;
orient(1,1) = theta0;
orient(1,2) = theta0+90;
theta1 = pi*theta1/180;
theta2 = pi*theta2/180;
Uh =  [cos(theta1) cos(theta2); ...
       sin(theta1) sin(theta2)];
Uh = inv(Uh);

switch Z.Nch
  case 5
     theta1 = Z.orient(1,4)+Z.theta0-theta0;
     theta2 = Z.orient(1,5)+Z.theta0-theta0;
     orient(1,4) = theta0;
     orient(1,5) = theta0+90;
     theta1 = pi*theta1/180;
     theta2 = pi*theta2/180;
     Ue = [1        0           0 ;
	   0  cos(theta1) cos(theta2); ...
           0  sin(theta1) sin(theta2)];

  case 4
     theta1 = Z.orient(1,3)+Z.theta0-theta0;
     theta2 = Z.orient(1,4)+Z.theta0-theta0;
     orient(1,3) = theta1;
     orient(1,4) = theta2;
     theta1 = pi*theta1/180;
     theta2 = pi*theta2/180;
     Ue =  [cos(theta1) cos(theta2); ...
       sin(theta1) sin(theta2)];
  otherwise
     Z = -1;
     return
end

TF = zeros(size(Z.TF));
SIG_E = zeros(size(Z.SIG_E));
SIG_S = zeros(size(Z.SIG_S));
for ib = 1:Z.nbt
   TF(:,:,ib) = Uh'*squeeze(Z.TF(:,:,ib))*Ue';
   SIG_E(:,:,ib) = Ue*squeeze(Z.SIG_E(:,:,ib))*Ue';
   SIG_S(:,:,ib) = Uh'*squeeze(Z.SIG_S(:,:,ib))*Uh;
end
Z.TF = TF;
Z.SIG_E = SIG_E;
Z.SIG_S = SIG_S;
Z.orient = orient;
Z.theta0 =  theta0;
