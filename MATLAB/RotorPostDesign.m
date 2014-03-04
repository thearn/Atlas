%----- RotorPostDesign -----%
%                           %
%       Todd Reichert       %
%        July 30, 2012       %
%                           %
%---------------------------%

% RotorPostDesign computes the material failure of the rotor post
% at the point of maximum bending moment based on the diameter, number of
% layers and wrap angle

% Max Bending moment
P = 1100; %Watts
RPM = 11; %RPM
v = RPM/60*pi*1.4;
F = P/v/4;
L = 1.111; %m
BM = F*L;

% Tube specifications
flags.CFRPType = 1;
t = 0.0055*2.54*10; %mm
d = (23+28)/2/1000; % spar diameter, mm
nTube = 23;
theta = 20*pi/180; % wrap angle, rad

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
f = [0 0 0 BM 0 0]';
qq = k\f;
qq = [zeros(6,1); qq];


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
Finternal(:,1) = [0 0 0 BM 0 0]';
flags.WireType = 1;
fail = FailureCalc([0 L], Finternal, strain, d, theta, nTube, 0, 1, 1, 1, 1, L, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, flags);
         

% Output
disp('------------------------------')
disp(sprintf('Force:   %10.2f N',F))
%disp(sprintf('Length:  %10.2f m',L))
disp(sprintf('Bending Moment:  %10.2f Nm',BM))

disp(sprintf('Diameter:%10.2f mm',d*1000))
disp(sprintf('Thickness:%10.2f mm',t*nTube))
disp(sprintf('EIx:     %10.2f ',EIx))
disp(sprintf('EIz:     %10.2f ',EIz))
disp(sprintf('EA:      %10.2f ',EA))
disp(sprintf('GJ:      %10.2f ',GJ))
disp(sprintf('Mass:    %10.2f kg\n',mSpar))

disp(sprintf('Fibre failure:  %5.2f ',fail.top.minus(1,1)));
disp(sprintf('Matrix failure: %5.2f ',fail.top.minus(2,1)));
disp(sprintf('Shear failure:  %5.2f ',fail.top.minus(3,1)));
