%----- DriveTubeDesign -----%
%                           %
%       Todd Reichert       %
%        May 29, 2012       %
%                           %
%---------------------------%

% DriveTubeDesign computes the torsional material failure and torsional
% buckling failure of the drive tube based on the diameter, number of
% layers and wrap angle

% Power delivery
P = 3000; % Watts
Omega = 0.2*2*pi; % rad/sec
Torque = P/Omega/4; % Nm

% % Tube specifications
% flags.CFRPType = 1;
% d = 0.5*2.54/100; % spar diameter, m
% nTube = 2;
% theta = 25*pi/180; % wrap angle, rad
% L = 0.755; % unsupported length

% Tube specifications
flags.CFRPType = 1; % HS40
d = 1.625*.0254; % spar diameter, m
nTube = 4;
theta = 45*pi/180; % wrap angle, rad
L = 0.8; % unsupported length

% Tube properties
[EIx, EIz, EA, GJ, mSpar] = SparProperties([0 L], d, theta, nTube, 0, L, flags);

% Displacement 
k = zeros(6,6);
k(1,1) = 12*EIx/(L*L*L);
k(1,6) = -6*EIx/(L*L);              
k(6,1) = k(1,6);
k(2,2) = EA/L;
k(3,3) = 12*EIz/(L*L*L);
k(3,4) = 6*EIz/(L*L);             
k(4,3) = k(3,4);
k(4,4) = 4*EIz/L;
k(5,5) = GJ/L;
k(6,6) = 4*EIx/L;
f = [0 0 0 0 Torque 0]';
qq = k\f;
qq = [zeros(6,1); qq];
%qq = [0 0 0 0 0 0 0 0 0 0 Torque*L/GJ 0]';


% Strains
x_hat = d/2;
z_hat = d/2;
r_hat = d/2;
strain.bending_x = -[(-(6*x_hat)/(L^2)) ((4*x_hat)/L) ((6*x_hat)/(L^2)) ((2*x_hat)/L)]*[qq(1) qq(6) qq(7) qq(12)]'; 
strain.bending_z = -[(-(6*z_hat)/(L^2)) ((-4*z_hat)/L) ((6*z_hat)/(L^2)) ((-2*z_hat)/L)]*[qq(3) qq(4) qq(9) qq(10)]'; 
strain.axial_y = [(-1/L) (1/L)]*[qq(2) qq(8)]'; 
strain.torsion_y = r_hat*[(-1/L) (1/L)]*[qq(5) qq(11)]';
    
strain.top(:,1) = [strain.bending_z+strain.axial_y 0 strain.torsion_y]';
strain.bottom(:,1) = [-strain.bending_z+strain.axial_y 0 strain.torsion_y]';
strain.back(:,1) = [strain.bending_x+strain.axial_y 0 strain.torsion_y]';
strain.front(:,1) = [-strain.bending_x+strain.axial_y 0 strain.torsion_y]';

% Failure calculation
Finternal(:,1) = [0 0 0 0 Torque 0]';
flags.WireType = 1;
fail = FailureCalc([0 L], Finternal, strain, d, theta, nTube, 0, 1, 1, 1, 1, L, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, flags);
         

% Output
disp('------------------------------')
disp(sprintf('Power:   %10.2f W',P))
disp(sprintf('Torque:  %10.2f Nm\n',Torque))

disp(sprintf('Diameter:%10.2f mm',d*1000))
disp(sprintf('EIx:     %10.2f ',EIx))
disp(sprintf('EIz:     %10.2f ',EIz))
disp(sprintf('EA:      %10.2f ',EA))
disp(sprintf('GJ:      %10.2f ',GJ))
disp(sprintf('Mass:    %10.2f kg\n',mSpar))

disp(sprintf('Fibre failure:  %5.2f ',fail.top.minus(1,1)));
disp(sprintf('Matrix failure: %5.2f ',fail.top.minus(2,1)));
disp(sprintf('Shear failure:  %5.2f ',fail.top.minus(3,1)));
disp(sprintf('Tor buckling:   %5.2f ',fail.buckling.torsion(1)));
disp(sprintf('Tor buckling (corr):   %5.2f ',0.742*fail.buckling.torsion(1)));



% Tube specifications
flags.CFRPType = 3; % MTM-28
d = 0.05; % spar diameter, m
nTube = 4;
theta = 45*pi/180; % wrap angle, rad
L = 0.8; % unsupported length
